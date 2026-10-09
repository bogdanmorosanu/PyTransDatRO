"""Pydantic schemas for grid information and telemetry stats.
"""
from typing import Dict, List, Optional
from pydantic import BaseModel, Field


class BoundingBox(BaseModel):
    min_northing: float
    min_easting: float
    max_northing: float
    max_easting: float


class GridInfoResponse(BaseModel):
    grid_file: str
    bounds: BoundingBox
    crs: Dict[str, str]
    quasigeoid_model: str


class DailyStatItem(BaseModel):
    date: str
    source: str
    calls: int
    points: int


class TelemetrySummary(BaseModel):
    total_calls: int
    total_points: int


class TelemetryStatsResponse(BaseModel):
    summary: TelemetrySummary
    daily_trend: List[DailyStatItem]
