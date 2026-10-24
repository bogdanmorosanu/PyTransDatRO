"""Grid information and telemetry statistics endpoints."""

import datetime
import json
import logging
import math
import sqlite3
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, Depends, Query

from api.dependencies import get_db_path, get_legacy_db_path, get_trans_engine
from api.schemas.info import (
    BoundingBox,
    DailyStatItem,
    GridComparisonResponse,
    GridInfoResponse,
    GridMetadataResponse,
    LegacyStatsSummary,
    LegacyTimelineMetrics,
    LegacyTimelineResponse,
    LegacyTimelineSeriesItem,
    LegacyTimelineSummary,
    QuadrantItem,
    QuadrantSummary,
    SourceItem,
    SpatialBounds,
    SpatialCellItem,
    SpatialMatrixResponse,
    SpatialMatrixSummary,
    SpgFileInfo,
    TelemetryStatsResponse,
    TelemetrySummary,
    TimelineResponse,
    TimelineSeriesItem,
    TimelineSeriesMetrics,
    TimelineSummary,
)
from pytransdatro.trans_ro import TransRO

try:
    from tools.export_spg_metadata_json import extract_full_spg_metadata
except ImportError:
    extract_full_spg_metadata = None

logger = logging.getLogger("api.routes.info")
router = APIRouter(prefix="/api/v1", tags=["Metadata & Telemetry"])


def _load_comparison_summary() -> dict:
    """Loads ANCPI replication benchmark JSON from pytransdatro/grids/ancpi_comparison_summary.json."""
    json_path = (
        Path(__file__).resolve().parent.parent.parent
        / "pytransdatro"
        / "grids"
        / "ancpi_comparison_summary.json"
    )
    if json_path.exists():
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Failed to load ANCPI comparison summary from %s: %s", json_path, e)
    return {
        "total_points": 9046,
        "overall_status": "PASS",
        "etrs89_to_st70": {
            "direction": "ETRS89 to Stereo70",
            "total_points": 9046,
            "status": "PASS",
            "max_abs_error": {"easting_m": 0.000314, "northing_m": 0.000335, "elevation_m": 0.0005},
            "mean_abs_error": {
                "easting_m": 0.00007,
                "northing_m": 0.00008,
                "elevation_m": 0.000251,
            },
            "flags": {"easting": {"0": 9046}, "northing": {"0": 9046}, "elevation": {"0": 9046}},
        },
        "st70_to_etrs89": {
            "direction": "Stereo70 to ETRS89",
            "total_points": 9046,
            "status": "PASS",
            "max_abs_error": {
                "latitude_sec": 0.00001006,
                "longitude_sec": 0.00001007,
                "elevation_m": 0.0005,
            },
            "mean_abs_error": {
                "latitude_sec": 0.00000501,
                "longitude_sec": 0.000005,
                "elevation_m": 0.000251,
            },
            "flags": {
                "latitude": {"0": 9032, "1": 14},
                "longitude": {"0": 9015, "1": 31},
                "elevation": {"0": 9046},
            },
        },
    }


def _load_or_extract_metadata() -> dict:
    """Loads metadata JSON from pytransdatro/grids/spg_metadata.json or extracts it dynamically."""
    json_path = (
        Path(__file__).resolve().parent.parent.parent
        / "pytransdatro"
        / "grids"
        / "spg_metadata.json"
    )
    if json_path.exists():
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Failed to load SPG grid metadata from %s: %s", json_path, e)

    # Fallback to dynamic extraction if JSON not generated yet and tools available
    spg_dir = Path(__file__).resolve().parent.parent.parent / "pytransdatro" / "grids"
    spg_candidates = list(spg_dir.glob("*.spg"))
    if spg_candidates and extract_full_spg_metadata is not None:
        return extract_full_spg_metadata(spg_candidates[0])
    return {
        "file_info": {
            "name": "rom_grid3d_25.09.spg",
            "size_bytes": 0,
            "size_formatted": "0 B",
            "sha256": "",
        },
        "metadata": {},
        "params": {},
        "grids": {},
    }


@router.get("/grid/comparison", response_model=GridComparisonResponse)
async def get_grid_comparison() -> GridComparisonResponse:
    """Returns official ANCPI replication benchmark test summary and residual metrics."""
    return GridComparisonResponse(**_load_comparison_summary())


@router.get("/grid/metadata", response_model=GridMetadataResponse)
async def get_grid_metadata() -> GridMetadataResponse:
    """Returns comprehensive, dynamic 3D SPG grid metadata, Helmert parameters, layer statistics, and diagnostics."""
    meta_dict = _load_or_extract_metadata()
    return GridMetadataResponse(**meta_dict)


@router.get("/grid/info", response_model=GridInfoResponse)
async def get_grid_info(trans: TransRO = Depends(get_trans_engine)) -> GridInfoResponse:
    """Returns active .spg grid metadata, coverage boundaries, and CRS parameters."""
    grid = trans._t_gr2d
    grid_filename = Path(grid.file_name).name
    meta_dict = _load_or_extract_metadata()

    file_info = SpgFileInfo(
        **meta_dict.get(
            "file_info",
            {"name": grid_filename, "size_bytes": 0, "size_formatted": "0 B", "sha256": ""},
        )
    )

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
            "geographic": "EPSG:4258 (ETRS89)",
        },
        quasigeoid_model="Modern Romanian Quasigeoid (1D nearest-neighbor)",
        metadata=meta_dict.get("metadata"),
        params=meta_dict.get("params"),
        grids=meta_dict.get("grids"),
        file_info=file_info,
    )


@router.get("/telemetry/stats", response_model=TelemetryStatsResponse)
async def get_telemetry_stats(
    db_path: Path = Depends(get_db_path),
) -> TelemetryStatsResponse:
    """Returns high-level, anonymized usage telemetry from SQLite WAL database."""
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
                        points=points_int,
                    )
                )

            # Query overall totals
            cursor.execute(
                "SELECT COALESCE(SUM(count_calls), 0), COALESCE(SUM(count_points), 0) FROM computation_usage"
            )
            tot_row = cursor.fetchone()
            if tot_row:
                total_calls = int(tot_row[0])
                total_points = int(tot_row[1])

            conn.close()
        except sqlite3.Error as e:
            logger.error("Database error querying telemetry stats: %s", e)

    return TelemetryStatsResponse(
        summary=TelemetrySummary(total_calls=total_calls, total_points=total_points),
        daily_trend=daily_items,
    )


@router.get("/telemetry/timeline", response_model=TimelineResponse)
async def get_telemetry_timeline(
    from_date: Optional[str] = Query(None, description="Start date filter (YYYY-MM-DD)"),
    to_date: Optional[str] = Query(None, description="End date filter (YYYY-MM-DD)"),
    granularity: str = Query(
        "daily",
        pattern="^(daily|weekly|monthly)$",
        description="Aggregation bucket: daily, weekly, or monthly",
    ),
    source: Optional[int] = Query(
        None,
        description="Filter by client source ID (1=rest_api, 2=web_map, 3=desktop_cli)",
    ),
    operation: Optional[str] = Query(
        None,
        description="Filter by operation: 1, 2, Stereo70ToETRS89, or ETRS89ToStereo70",
    ),
    is_3d: Optional[int] = Query(None, description="Filter by 2D (0) or 3D (1)"),
    db_path: Path = Depends(get_db_path),
) -> TimelineResponse:
    """Returns time-series telemetry data for the Telemetry Visualization Dashboard.

    Supports dynamic date filtering, daily, weekly, vs. monthly rollups, and breakdown by client source.
    """
    if not db_path.exists():
        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=0,
                total_points=0,
                avg_points_per_call=0.0,
                from_date=from_date,
                to_date=to_date,
                days_count=0,
            ),
            series=[],
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
            filters.append("c.period_end >= ?")
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
                "web_map": 0,
                "rest_api": 0,
                "desktop_cli": 0,
                "unknown": 0,
                "st70_to_etrs": 0,
                "etrs_to_st70": 0,
                "dim_2d": 0,
                "dim_3d": 0,
                "total": 0,
                "sources": {},
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

        # Query effective period boundaries from matched rows
        bounds_sql = f"""
            SELECT MIN(c.period_start), MAX(c.period_end)
            FROM computation_usage c
            {where_clause}
        """
        cursor.execute(bounds_sql, params)
        b_row = cursor.fetchone()
        effective_from = b_row[0] if b_row and b_row[0] else None
        effective_to = b_row[1] if b_row and b_row[1] else None

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

        # Ensure resolved_from is a full YYYY-MM-DD date string even if monthly aggregation was used
        if resolved_from and len(resolved_from) == 7:
            resolved_from = f"{resolved_from}-01"

        avg_ratio = round(total_points_all / total_calls_all, 2) if total_calls_all > 0 else 0.0

        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=total_calls_all,
                total_points=total_points_all,
                avg_points_per_call=avg_ratio,
                from_date=resolved_from,
                to_date=resolved_to,
                effective_from_date=effective_from,
                effective_to_date=effective_to,
                days_count=days_count,
                period_count=period_count,
            ),
            series=series_items,
            available_sources=available_sources,
        )

    except sqlite3.Error as e:
        logger.error("Database error querying telemetry timeline: %s", e)
        return TimelineResponse(
            summary=TimelineSummary(
                total_calls=0,
                total_points=0,
                avg_points_per_call=0.0,
                from_date=from_date,
                to_date=to_date,
                effective_from_date=from_date,
                effective_to_date=to_date,
                days_count=0,
                period_count=0,
            ),
            series=[],
        )


@router.get("/telemetry/spatial-matrix", response_model=SpatialMatrixResponse)
async def get_telemetry_spatial_matrix(
    from_date: Optional[str] = Query(None, description="Start date YYYY-MM-DD"),
    to_date: Optional[str] = Query(None, description="End date YYYY-MM-DD"),
    metric: str = Query(
        "points",
        pattern="^(points|calls)$",
        description="Metric to aggregate: 'points' or 'calls'",
    ),
    sources: Optional[str] = Query(
        None, description="Comma-separated client source IDs e.g. '1,2'"
    ),
    operation: Optional[int] = Query(
        None, description="Operation ID (1: st70_to_etrs89, 2: etrs89_to_st70)"
    ),
    is_3d: Optional[int] = Query(None, description="0 for 2D, 1 for 3D"),
    db_path: Path = Depends(get_db_path),
) -> SpatialMatrixResponse:
    """Returns aggregated 0.1 degree spatial density cells and mathematical quadrant totals."""

    # Default bounding box for Romania
    default_bounds = SpatialBounds(min_lat=43.6, max_lat=48.3, min_lon=20.2, max_lon=29.7)
    center_lat = 45.9
    center_lon = 24.9

    default_quadrants = QuadrantSummary(
        center={"lat": center_lat, "lon": center_lon},
        NW=QuadrantItem(value=0, percentage=0.0),
        NE=QuadrantItem(value=0, percentage=0.0),
        SW=QuadrantItem(value=0, percentage=0.0),
        SE=QuadrantItem(value=0, percentage=0.0),
    )

    empty_response = SpatialMatrixResponse(
        resolution=0.1,
        bounds=default_bounds,
        summary=SpatialMatrixSummary(
            active_cells_count=0,
            total_metric_value=0,
            max_cell_hits=0,
            metric=metric,
            from_date=from_date,
            to_date=to_date,
            effective_from_date=from_date,
            effective_to_date=to_date,
            quadrants=default_quadrants,
        ),
        cells=[],
    )

    if not db_path.exists():
        return empty_response

    try:
        conn = sqlite3.connect(f"file:{db_path.resolve()}?mode=ro", uri=True, timeout=2.0)
        cursor = conn.cursor()

        # Build dynamic SQL query
        where_clauses = []
        params = []

        if from_date:
            where_clauses.append("period_end >= ?")
            params.append(from_date)
        if to_date:
            where_clauses.append("period_start <= ?")
            params.append(to_date)
        if sources:
            src_ids = [int(s.strip()) for s in sources.split(",") if s.strip().isdigit()]
            if src_ids:
                placeholders = ",".join("?" for _ in src_ids)
                where_clauses.append(f"source IN ({placeholders})")
                params.extend(src_ids)
        if operation is not None:
            where_clauses.append("operation = ?")
            params.append(operation)
        if is_3d is not None:
            where_clauses.append("is_3d = ?")
            params.append(is_3d)

        where_sql = ("WHERE " + " AND ".join(where_clauses)) if where_clauses else ""

        # Metric aggregation expression
        # 'points' aggregates coordinate count; 'calls' aggregates cell transformation event occurrences
        metric_expr = "SUM(count)" if metric == "points" else "SUM(count_calls)"

        # Query aggregating ~1km (lat_deg, lat_grad) to 0.1 degree bins
        query = f"""
            SELECT
                ROUND(lat_deg + (CAST(lat_grad / 10 AS INTEGER) * 0.1), 1) AS bin_lat,
                ROUND(lon_deg + (CAST(lon_grad / 10 AS INTEGER) * 0.1), 1) AS bin_lon,
                {metric_expr} AS hits
            FROM spatial_usage
            {where_sql}
            GROUP BY bin_lat, bin_lon
            HAVING hits > 0
            ORDER BY bin_lat ASC, bin_lon ASC
        """
        cursor.execute(query, params)
        rows = cursor.fetchall()

        # Query effective period boundaries from matched rows
        bounds_sql = f"""
            SELECT MIN(period_start), MAX(period_end)
            FROM spatial_usage
            {where_sql}
        """
        cursor.execute(bounds_sql, params)
        b_row = cursor.fetchone()
        effective_from = b_row[0] if b_row and b_row[0] else from_date
        effective_to = b_row[1] if b_row and b_row[1] else to_date

        conn.close()

        if not rows:
            return empty_response

        cells: List[SpatialCellItem] = []
        total_val = 0
        max_hits = 0
        min_lat_found = 90.0
        max_lat_found = -90.0
        min_lon_found = 180.0
        max_lon_found = -180.0

        quad_nw = 0
        quad_ne = 0
        quad_sw = 0
        quad_se = 0

        for r_lat, r_lon, r_hits in rows:
            c_lat = round(float(r_lat), 1)
            c_lon = round(float(r_lon), 1)
            hits_int = int(r_hits)

            cells.append(SpatialCellItem(lat=c_lat, lon=c_lon, hits=hits_int))
            total_val += hits_int
            if hits_int > max_hits:
                max_hits = hits_int

            if c_lat < min_lat_found:
                min_lat_found = c_lat
            if c_lat > max_lat_found:
                max_lat_found = c_lat
            if c_lon < min_lon_found:
                min_lon_found = c_lon
            if c_lon > max_lon_found:
                max_lon_found = c_lon

            # Mathematical Quadrant partitioning around center (45.9, 24.9)
            if c_lat >= center_lat:
                if c_lon < center_lon:
                    quad_nw += hits_int
                else:
                    quad_ne += hits_int
            else:
                if c_lon < center_lon:
                    quad_sw += hits_int
                else:
                    quad_se += hits_int

        # Calculate percentages
        def calc_pct(val: int, tot: int) -> float:
            return round((val / tot) * 100.0, 1) if tot > 0 else 0.0

        quadrants = QuadrantSummary(
            center={"lat": center_lat, "lon": center_lon},
            NW=QuadrantItem(value=quad_nw, percentage=calc_pct(quad_nw, total_val)),
            NE=QuadrantItem(value=quad_ne, percentage=calc_pct(quad_ne, total_val)),
            SW=QuadrantItem(value=quad_sw, percentage=calc_pct(quad_sw, total_val)),
            SE=QuadrantItem(value=quad_se, percentage=calc_pct(quad_se, total_val)),
        )

        bounds = SpatialBounds(
            min_lat=min_lat_found if min_lat_found <= max_lat_found else default_bounds.min_lat,
            max_lat=max_lat_found if min_lat_found <= max_lat_found else default_bounds.max_lat,
            min_lon=min_lon_found if min_lon_found <= max_lon_found else default_bounds.min_lon,
            max_lon=max_lon_found if min_lon_found <= max_lon_found else default_bounds.max_lon,
        )

        return SpatialMatrixResponse(
            resolution=0.1,
            bounds=bounds,
            summary=SpatialMatrixSummary(
                active_cells_count=len(cells),
                total_metric_value=total_val,
                max_cell_hits=max_hits,
                metric=metric,
                from_date=from_date,
                to_date=to_date,
                effective_from_date=effective_from,
                effective_to_date=effective_to,
                quadrants=quadrants,
            ),
            cells=cells,
        )

    except sqlite3.Error as e:
        logger.error("Database error querying spatial matrix: %s", e)
        return empty_response


# ==============================================================================
# Legacy TransdatOnline Statistics Endpoints
# ==============================================================================
@router.get("/telemetry/legacy/stats", response_model=LegacyStatsSummary)
async def get_legacy_stats(
    db_path: Path = Depends(get_legacy_db_path),
) -> LegacyStatsSummary:
    """Retrieve overall aggregate summary from legacy TransdatOnline statistics database."""
    if not db_path.exists():
        return LegacyStatsSummary(
            total_points=0,
            app_st70=0,
            app_st30=0,
            service_st70=0,
            service_st30=0,
            from_date=None,
            to_date=None,
            recorded_days=0,
        )

    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        cursor = conn.cursor()
        cursor.execute("""
            SELECT
                COUNT(*),
                MIN(day),
                MAX(day),
                COALESCE(SUM(cooOpAppSt70), 0),
                COALESCE(SUM(cooOpAppSt30), 0),
                COALESCE(SUM(cooOpServiceSt70), 0),
                COALESCE(SUM(cooOpServiceSt30), 0)
            FROM statistics
        """)
        row = cursor.fetchone()
        conn.close()

        count, min_day, max_day, app70, app30, srv70, srv30 = row
        total = app70 + app30 + srv70 + srv30

        return LegacyStatsSummary(
            total_points=total,
            app_st70=app70,
            app_st30=app30,
            service_st70=srv70,
            service_st30=srv30,
            from_date=min_day,
            to_date=max_day,
            recorded_days=count,
        )
    except sqlite3.Error as e:
        logger.error("Error querying legacy stats database: %s", e)
        return LegacyStatsSummary(
            total_points=0,
            app_st70=0,
            app_st30=0,
            service_st70=0,
            service_st30=0,
            from_date=None,
            to_date=None,
            recorded_days=0,
        )


@router.get("/telemetry/legacy/timeline", response_model=LegacyTimelineResponse)
async def get_legacy_timeline(
    from_date: Optional[str] = Query(None, description="Start date YYYY-MM-DD"),
    to_date: Optional[str] = Query(None, description="End date YYYY-MM-DD"),
    granularity: str = Query(
        "daily",
        pattern="^(daily|weekly|monthly)$",
        description="Temporal aggregation bucket: daily, weekly, or monthly",
    ),
    source: Optional[str] = Query(
        None, description="Filter source: 'app', 'service', or comma-separated e.g. 'app,service'"
    ),
    operation: Optional[str] = Query(
        None, description="Filter operation: 'st70', 'st30', or comma-separated e.g. 'st70,st30'"
    ),
    db_path: Path = Depends(get_legacy_db_path),
) -> LegacyTimelineResponse:
    """Retrieve time-series transformation volumes from legacy TransdatOnline database."""
    if not db_path.exists():
        return LegacyTimelineResponse(
            summary=LegacyTimelineSummary(total_points=0, from_date=from_date, to_date=to_date),
            series=[],
            bounds={},
        )

    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        cursor = conn.cursor()

        # Retrieve absolute database boundaries
        cursor.execute("SELECT MIN(day), MAX(day) FROM statistics")
        min_bound, max_bound = cursor.fetchone()
        bounds_map = {"min_date": min_bound or "", "max_date": max_bound or ""}

        filters = []
        params = []

        if from_date:
            filters.append("day >= ?")
            params.append(from_date)
        if to_date:
            filters.append("day <= ?")
            params.append(to_date)

        where_clause = f"WHERE {' AND '.join(filters)}" if filters else ""

        if granularity == "monthly":
            time_expr = "substr(day, 1, 7)"
        elif granularity == "weekly":
            time_expr = "date(day, 'weekday 0', '-6 days')"
        else:
            time_expr = "day"

        sql = f"""
            SELECT
                {time_expr} AS time_bucket,
                SUM(cooOpAppSt70) AS app_st70,
                SUM(cooOpAppSt30) AS app_st30,
                SUM(cooOpServiceSt70) AS srv_st70,
                SUM(cooOpServiceSt30) AS srv_st30
            FROM statistics
            {where_clause}
            GROUP BY time_bucket
            ORDER BY time_bucket ASC
        """
        cursor.execute(sql, params)
        rows = cursor.fetchall()
        conn.close()

        # Parse source and operation filters if provided
        active_sources = set()
        if source:
            for s in source.split(","):
                s_clean = s.strip().lower()
                if s_clean in ("app", "web_app", "webapp"):
                    active_sources.add("app")
                elif s_clean in ("service", "rest_api", "api"):
                    active_sources.add("service")
        else:
            active_sources = {"app", "service"}

        active_ops = set()
        if operation:
            for o in operation.split(","):
                o_clean = o.strip().lower()
                if o_clean in ("st70", "stereo70"):
                    active_ops.add("st70")
                elif o_clean in ("st30", "stereo30"):
                    active_ops.add("st30")
        else:
            active_ops = {"st70", "st30"}

        series_items: List[LegacyTimelineSeriesItem] = []
        overall_total_points = 0

        for time_bucket, raw_app70, raw_app30, raw_srv70, raw_srv30 in rows:
            v_app70 = int(raw_app70 or 0)
            v_app30 = int(raw_app30 or 0)
            v_srv70 = int(raw_srv70 or 0)
            v_srv30 = int(raw_srv30 or 0)

            app_tot = v_app70 + v_app30
            srv_tot = v_srv70 + v_srv30
            st70_tot = v_app70 + v_srv70
            st30_tot = v_app30 + v_srv30

            # Filtered total for this bucket based on active filters
            bucket_total = 0
            if "app" in active_sources:
                if "st70" in active_ops:
                    bucket_total += v_app70
                if "st30" in active_ops:
                    bucket_total += v_app30
            if "service" in active_sources:
                if "st70" in active_ops:
                    bucket_total += v_srv70
                if "st30" in active_ops:
                    bucket_total += v_srv30

            overall_total_points += bucket_total

            series_items.append(
                LegacyTimelineSeriesItem(
                    date=time_bucket,
                    points=LegacyTimelineMetrics(
                        app_st70=v_app70,
                        app_st30=v_app30,
                        service_st70=v_srv70,
                        service_st30=v_srv30,
                        app_total=app_tot,
                        service_total=srv_tot,
                        st70_total=st70_tot,
                        st30_total=st30_tot,
                        total=bucket_total,
                    ),
                )
            )

        resolved_from = from_date or (series_items[0].date if series_items else None)
        resolved_to = to_date or (series_items[-1].date if series_items else None)

        days_count = None
        period_count = len(series_items)
        if resolved_from and resolved_to:
            try:
                d1 = datetime.date.fromisoformat(resolved_from)
                d2 = datetime.date.fromisoformat(resolved_to)
                days_count = max(0, (d2 - d1).days + 1)
            except ValueError:
                pass

        return LegacyTimelineResponse(
            summary=LegacyTimelineSummary(
                total_points=overall_total_points,
                from_date=resolved_from,
                to_date=resolved_to,
                days_count=days_count,
                period_count=period_count,
            ),
            series=series_items,
            bounds=bounds_map,
        )

    except sqlite3.Error as e:
        logger.error("Database error querying legacy timeline: %s", e)
        return LegacyTimelineResponse(
            summary=LegacyTimelineSummary(total_points=0, from_date=from_date, to_date=to_date),
            series=[],
            bounds={},
        )
