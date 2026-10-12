"""Grid information and telemetry statistics endpoints.
"""
import sqlite3
import datetime
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
    granularity: str = Query("daily", pattern="^(daily|monthly)$", description="Aggregation bucket: daily or monthly"),
    source: Optional[int] = Query(None, description="Filter by client source ID (1=rest_api, 2=web_map, 3=desktop_cli)"),
    operation: Optional[str] = Query(None, description="Filter by operation: 1, 2, Stereo70ToETRS89, or ETRS89ToStereo70"),
    is_3d: Optional[int] = Query(None, description="Filter by 2D (0) or 3D (1)"),
) -> TimelineResponse:
    """Returns time-series telemetry data for the Telemetry Visualization Dashboard.
    
    Supports dynamic date filtering, daily vs. monthly rollups, and breakdown by client source.
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

        # Time grouping: 'daily' groups by period_start; 'monthly' groups by substr(period_start, 1, 7)
        time_expr = "substr(c.period_start, 1, 7)" if granularity == "monthly" else "c.period_start"

        sql = f"""
            SELECT 
                {time_expr} AS time_bucket,
                COALESCE(s.name, 'unknown') AS source_name,
                SUM(c.count_calls) AS total_calls,
                SUM(c.count_points) AS total_points
            FROM computation_usage c
            LEFT JOIN sources s ON c.source = s.id
            {where_clause}
            GROUP BY time_bucket, c.source
            ORDER BY time_bucket ASC, c.source ASC
        """

        cursor.execute(sql, params)
        rows = cursor.fetchall()

        # Aggregate series buckets
        series_map = {}
        total_calls_all = 0
        total_points_all = 0
        all_dates = set()

        for time_bucket, src_name, calls, points in rows:
            c_int = int(calls)
            p_int = int(points)
            total_calls_all += c_int
            total_points_all += p_int
            all_dates.add(time_bucket)

            if time_bucket not in series_map:
                series_map[time_bucket] = {
                    "calls": {"web_map": 0, "rest_api": 0, "desktop_cli": 0, "unknown": 0, "total": 0},
                    "points": {"web_map": 0, "rest_api": 0, "desktop_cli": 0, "unknown": 0, "total": 0},
                }

            # Map known source names
            clean_src = src_name if src_name in ("web_map", "rest_api", "desktop_cli") else "unknown"
            series_map[time_bucket]["calls"][clean_src] += c_int
            series_map[time_bucket]["calls"]["total"] += c_int
            series_map[time_bucket]["points"][clean_src] += p_int
            series_map[time_bucket]["points"]["total"] += p_int

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

        # Calculate days_count / date boundaries
        resolved_from = from_date
        resolved_to = to_date
        days_count = None

        if resolved_from and resolved_to:
            try:
                d1 = datetime.date.fromisoformat(resolved_from)
                d2 = datetime.date.fromisoformat(resolved_to)
                days_count = max(0, (d2 - d1).days + 1)
            except ValueError:
                pass
        elif series_items:
            resolved_from = series_items[0].date
            resolved_to = series_items[-1].date
            days_count = len(all_dates)

        avg_ratio = round(total_points_all / total_calls_all, 2) if total_calls_all > 0 else 0.0

        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=total_calls_all,
                total_points=total_points_all,
                avg_points_per_call=avg_ratio,
                from_date=resolved_from,
                to_date=resolved_to,
                days_count=days_count
            ),
            series=series_items
        )

    except Exception:
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
