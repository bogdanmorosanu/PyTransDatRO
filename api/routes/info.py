"""Grid information and telemetry statistics endpoints.
"""
import sqlite3
import datetime
import math
from pathlib import Path
from typing import List, Optional
from fastapi import APIRouter, Depends, Query, Request
from api.dependencies import get_trans_engine
from api.schemas.info import (
    BoundingBox,
    GridInfoResponse,
    DailyStatItem,
    TelemetrySummary,
    TelemetryStatsResponse,
    TimelineResponse,
    TimelineSummary,
    TimelineSeriesItem,
    TimelineSeriesMetrics,
    SourceItem,
)
from pytransdatro.trans_ro import TransRO

router = APIRouter(prefix="/api/v1", tags=["Metadata & Telemetry"])


@router.get("/grid/info", response_model=GridInfoResponse)
async def get_grid_info(trans: TransRO = Depends(get_trans_engine)) -> GridInfoResponse:
    """Returns active .spg grid metadata, coverage boundaries, and CRS parameters."""
    grid = trans._t_gr2d
    grid_filename = Path(grid.file_name).name

    return GridInfoResponse(
        grid_file=grid_filename,
        bounds=BoundingBox(
            min_northing=float(grid.n_min),
            min_easting=float(grid.e_min),
            max_northing=float(grid.n_max),
            max_easting=float(grid.e_max),
        ),
        crs={
            "projected": "EPSG:31700 (Pulkovo 1942(58) / Stereo70)",
            "geographic": "EPSG:4258 (ETRS89)"
        },
        quasigeoid_model="Modern Romanian Quasigeoid (1D nearest-neighbor)"
    )


@router.get("/telemetry/stats", response_model=TelemetryStatsResponse)
async def get_telemetry_stats(request: Request) -> TelemetryStatsResponse:
    """Returns high-level, anonymized usage telemetry from SQLite WAL database."""
    log_dir = Path("./logs")
    db_path = log_dir / "usage.db"

    total_calls = 0
    total_points = 0
    daily_items: List[DailyStatItem] = []

    if db_path.exists():
        try:
            # Connect in read-only mode to prevent locking
            conn = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True, timeout=2.0)
            cursor = conn.cursor()

            # Query computation_usage (grouped by period_start and source)
            cursor.execute(
                """
                SELECT c.period_start, COALESCE(s.name, 'unknown'), SUM(c.count_calls), SUM(c.count_points) 
                FROM computation_usage c
                LEFT JOIN sources s ON c.source = s.id
                GROUP BY c.period_start, c.source
                ORDER BY c.period_start DESC, c.source ASC
                LIMIT 30
                """
            )
            rows = cursor.fetchall()

            for date_str, src_name, calls, points in rows:
                calls_int = int(calls)
                points_int = int(points)
                daily_items.append(
                    DailyStatItem(
                        date=date_str,
                        source=src_name,
                        calls=calls_int,
                        points=points_int
                    )
                )

            # Query overall totals
            cursor.execute("SELECT COALESCE(SUM(count_calls), 0), COALESCE(SUM(count_points), 0) FROM computation_usage")
            tot_row = cursor.fetchone()
            if tot_row:
                total_calls = int(tot_row[0])
                total_points = int(tot_row[1])

            conn.close()
        except Exception:
            # If reading telemetry database fails for any reason, return default
            pass

    return TelemetryStatsResponse(
        summary=TelemetrySummary(total_calls=total_calls, total_points=total_points),
        daily_trend=daily_items
    )


@router.get("/telemetry/timeline", response_model=TimelineResponse)
async def get_telemetry_timeline(
    from_date: Optional[str] = Query(None, description="Start date filter (YYYY-MM-DD)"),
    to_date: Optional[str] = Query(None, description="End date filter (YYYY-MM-DD)"),
    granularity: str = Query("daily", pattern="^(daily|weekly|monthly)$", description="Aggregation bucket: daily, weekly, or monthly"),
    source: Optional[int] = Query(None, description="Filter by client source ID (1=rest_api, 2=web_map, 3=desktop_cli)"),
    operation: Optional[str] = Query(None, description="Filter by operation: 1, 2, Stereo70ToETRS89, or ETRS89ToStereo70"),
    is_3d: Optional[int] = Query(None, description="Filter by 2D (0) or 3D (1)"),
) -> TimelineResponse:
    """Returns time-series telemetry data for the Telemetry Visualization Dashboard.
    
    Supports dynamic date filtering, daily, weekly, vs. monthly rollups, and breakdown by client source.
    """
    log_dir = Path("./logs")
    db_path = log_dir / "usage.db"

    if not db_path.exists():
        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=0,
                total_points=0,
                avg_points_per_call=0.0,
                from_date=from_date,
                to_date=to_date,
                days_count=0
            ),
            series=[]
        )

    try:
        conn = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True, timeout=2.0)
        cursor = conn.cursor()

        # Query all available registered sources
        cursor.execute("SELECT id, name, description FROM sources ORDER BY id ASC")
        source_rows = cursor.fetchall()
        available_sources = [
            SourceItem(id=s_id, name=s_name, description=s_desc)
            for s_id, s_name, s_desc in source_rows
        ]

        # Build dynamic query filters
        filters = []
        params = []

        if from_date:
            filters.append("c.period_start >= ?")
            params.append(from_date)
        if to_date:
            filters.append("c.period_start <= ?")
            params.append(to_date)
        if source is not None:
            filters.append("c.source = ?")
            params.append(source)
        if is_3d is not None:
            filters.append("c.is_3d = ?")
            params.append(1 if is_3d else 0)

        # Normalize operation filter
        if operation is not None:
            op_str = str(operation).strip()
            if op_str in ("1", "st70_to_etrs89", "Stereo70ToETRS89"):
                filters.append("c.operation = 1")
            elif op_str in ("2", "etrs89_to_st70", "ETRS89ToStereo70"):
                filters.append("c.operation = 2")

        where_clause = f"WHERE {' AND '.join(filters)}" if filters else ""

        # Time grouping: 'daily' groups by period_start; 'weekly' groups by date(period_start, 'weekday 0', '-6 days'); 'monthly' groups by substr(period_start, 1, 7)
        if granularity == "monthly":
            time_expr = "substr(c.period_start, 1, 7)"
        elif granularity == "weekly":
            time_expr = "date(c.period_start, 'weekday 0', '-6 days')"
        else:
            time_expr = "c.period_start"

        sql = f"""
            SELECT 
                {time_expr} AS time_bucket,
                COALESCE(s.name, 'unknown') AS source_name,
                c.operation,
                c.is_3d,
                SUM(c.count_calls) AS total_calls,
                SUM(c.count_points) AS total_points
            FROM computation_usage c
            LEFT JOIN sources s ON c.source = s.id
            {where_clause}
            GROUP BY time_bucket, c.source, c.operation, c.is_3d
            ORDER BY time_bucket ASC
        """

        cursor.execute(sql, params)
        rows = cursor.fetchall()

        def create_empty_metrics():
            return {
                "web_map": 0, "rest_api": 0, "desktop_cli": 0, "unknown": 0,
                "st70_to_etrs": 0, "etrs_to_st70": 0,
                "dim_2d": 0, "dim_3d": 0,
                "total": 0,
                "sources": {}
            }

        # Aggregate series buckets
        series_map = {}
        total_calls_all = 0
        total_points_all = 0
        all_dates = set()

        for time_bucket, src_name, op_val, is_3d_val, calls, points in rows:
            c_int = int(calls)
            p_int = int(points)
            total_calls_all += c_int
            total_points_all += p_int
            all_dates.add(time_bucket)

            if time_bucket not in series_map:
                series_map[time_bucket] = {
                    "calls": create_empty_metrics(),
                    "points": create_empty_metrics(),
                }

            b_calls = series_map[time_bucket]["calls"]
            b_points = series_map[time_bucket]["points"]

            # Map source dynamically into sources dict
            src_key = str(src_name)
            b_calls["sources"][src_key] = b_calls["sources"].get(src_key, 0) + c_int
            b_points["sources"][src_key] = b_points["sources"].get(src_key, 0) + p_int

            # Also map into named attribute for backward compatibility
            if src_key in b_calls:
                b_calls[src_key] += c_int
            else:
                b_calls[src_key] = c_int

            if src_key in b_points:
                b_points[src_key] += p_int
            else:
                b_points[src_key] = p_int

            # Map operation (1=st70_to_etrs, 2=etrs_to_st70)
            if op_val == 1:
                b_calls["st70_to_etrs"] += c_int
                b_points["st70_to_etrs"] += p_int
            elif op_val == 2:
                b_calls["etrs_to_st70"] += c_int
                b_points["etrs_to_st70"] += p_int

            # Map 2D / 3D
            if is_3d_val == 1:
                b_calls["dim_3d"] += c_int
                b_points["dim_3d"] += p_int
            else:
                b_calls["dim_2d"] += c_int
                b_points["dim_2d"] += p_int

            # Total for the bucket
            b_calls["total"] += c_int
            b_points["total"] += p_int

        conn.close()

        # Format series items sorted chronologically
        series_items: List[TimelineSeriesItem] = []
        for bucket in sorted(series_map.keys()):
            b_data = series_map[bucket]
            series_items.append(
                TimelineSeriesItem(
                    date=bucket,
                    calls=TimelineSeriesMetrics(**b_data["calls"]),
                    points=TimelineSeriesMetrics(**b_data["points"]),
                )
            )

        # Calculate days_count, period_count, and date boundaries
        resolved_from = from_date
        resolved_to = to_date
        days_count = None
        period_count = None

        if resolved_from and resolved_to:
            try:
                d1 = datetime.date.fromisoformat(resolved_from)
                d2 = datetime.date.fromisoformat(resolved_to)
                days_diff = max(0, (d2 - d1).days + 1)
                days_count = days_diff
                if granularity == "monthly":
                    period_count = max(1, (d2.year - d1.year) * 12 + (d2.month - d1.month) + 1)
                elif granularity == "weekly":
                    period_count = max(1, math.ceil(days_diff / 7.0))
                else:
                    period_count = days_diff
            except ValueError:
                pass
        elif series_items:
            resolved_from = series_items[0].date
            resolved_to = series_items[-1].date
            days_count = len(all_dates)
            period_count = len(series_items)

        avg_ratio = round(total_points_all / total_calls_all, 2) if total_calls_all > 0 else 0.0

        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=total_calls_all,
                total_points=total_points_all,
                avg_points_per_call=avg_ratio,
                from_date=resolved_from,
                to_date=resolved_to,
                days_count=days_count,
                period_count=period_count
            ),
            series=series_items,
            available_sources=available_sources
        )

    except Exception:
        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=0,
                total_points=0,
                avg_points_per_call=0.0,
                from_date=from_date,
                to_date=to_date,
                days_count=0,
                period_count=0
            ),
            series=[]
        )
