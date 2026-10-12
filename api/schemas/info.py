"""Pydantic schemas for grid information and telemetry stats/timeline.
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


# Timeline schemas adhering to specs/telemetry_visualization.md
class TimelineSummary(BaseModel):
    total_calls: int
    total_points: int
    avg_points_per_call: float
    from_date: Optional[str] = None
    to_date: Optional[str] = None
    days_count: Optional[int] = None
    period_count: Optional[int] = None


class SourceItem(BaseModel):
    id: int
    name: str
    description: Optional[str] = None


class TimelineSeriesMetrics(BaseModel):
    # Allow dynamic source names as attributes if needed
    model_config = {"extra": "allow"}

    # Sources (well-known defaults)
    web_map: int = 0
    rest_api: int = 0
    desktop_cli: int = 0
    unknown: int = 0
    # Operations
    st70_to_etrs: int = 0
    etrs_to_st70: int = 0
    # Dimensionality
    dim_2d: int = 0
    dim_3d: int = 0
    # Total
    total: int = 0
    # Dynamic dictionary of all sources for this bucket {source_name: count}
    sources: Dict[str, int] = Field(default_factory=dict)


class TimelineSeriesItem(BaseModel):
    date: str
    calls: TimelineSeriesMetrics
    points: TimelineSeriesMetrics


class TimelineResponse(BaseModel):
    summary: TimelineSummary
    series: List[TimelineSeriesItem]
    available_sources: List[SourceItem] = Field(default_factory=list)

