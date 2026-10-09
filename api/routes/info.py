"""Grid information and telemetry statistics endpoints.
"""
import sqlite3
from pathlib import Path
from typing import List
from fastapi import APIRouter, Depends, Request
from api.dependencies import get_trans_engine
from api.schemas.info import (
    BoundingBox,
    GridInfoResponse,
    DailyStatItem,
    TelemetrySummary,
    TelemetryStatsResponse,
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

            # Query daily stats
            cursor.execute(
                """
                SELECT d.date, s.name, d.count_calls, d.count_points 
                FROM daily_stats d
                JOIN sources s ON d.source = s.id
                ORDER BY d.date DESC, d.source ASC
                LIMIT 30
                """
            )
            rows = cursor.fetchall()

            for date_str, src_name, calls, points in rows:
                daily_items.append(
                    DailyStatItem(
                        date=date_str,
                        source=src_name,
                        calls=int(calls),
                        points=int(points)
                    )
                )
                total_calls += int(calls)
                total_points += int(points)

            conn.close()
        except Exception:
            # If reading telemetry database fails for any reason, return default
            pass

    return TelemetryStatsResponse(
        summary=TelemetrySummary(total_calls=total_calls, total_points=total_points),
        daily_trend=daily_items
    )
