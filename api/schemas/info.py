"""Pydantic schemas for grid information and telemetry stats/timeline."""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


class BoundingBox(BaseModel):
    min_northing: float
    min_easting: float
    max_northing: float
    max_easting: float


class SpgFileInfo(BaseModel):
    name: str
    size_bytes: int
    size_formatted: str
    sha256: str
    model_config = {"extra": "allow"}


class GridMetadataResponse(BaseModel):
    """Exhaustive, dynamic SPG grid metadata model for public REST consumers and web app."""

    file_info: SpgFileInfo
    metadata: Dict[str, Any]
    params: Dict[str, Any]
    grids: Dict[str, Any]
    model_config = {"extra": "allow"}


class GridInfoResponse(BaseModel):
    grid_file: str
    bounds: BoundingBox
    crs: Dict[str, str]
    quasigeoid_model: str
    metadata: Optional[Dict[str, Any]] = None
    params: Optional[Dict[str, Any]] = None
    grids: Optional[Dict[str, Any]] = None
    file_info: Optional[SpgFileInfo] = None
    model_config = {"extra": "allow"}


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
    effective_from_date: Optional[str] = None
    effective_to_date: Optional[str] = None
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


# Spatial Matrix schemas adhering to specs/telemetry_visualization.md
class SpatialCellItem(BaseModel):
    lat: float
    lon: float
    hits: int


class QuadrantItem(BaseModel):
    value: int
    percentage: float


class QuadrantSummary(BaseModel):
    center: Dict[str, float] = Field(default_factory=lambda: {"lat": 45.9, "lon": 24.9})
    NW: QuadrantItem
    NE: QuadrantItem
    SW: QuadrantItem
    SE: QuadrantItem


class SpatialBounds(BaseModel):
    min_lat: float
    max_lat: float
    min_lon: float
    max_lon: float


class SpatialMatrixSummary(BaseModel):
    active_cells_count: int
    total_metric_value: int
    max_cell_hits: int
    metric: str
    from_date: Optional[str] = None
    to_date: Optional[str] = None
    effective_from_date: Optional[str] = None
    effective_to_date: Optional[str] = None
    quadrants: QuadrantSummary


class SpatialMatrixResponse(BaseModel):
    resolution: float = 0.1
    bounds: SpatialBounds
    summary: SpatialMatrixSummary
    cells: List[SpatialCellItem]


# Benchmark comparison schemas for the Web App Help section
class DirectionComparisonSummary(BaseModel):
    direction: str
    total_points: int
    status: str
    max_abs_error: Dict[str, float]
    mean_abs_error: Dict[str, float]
    flags: Dict[str, Dict[str, int]]
    model_config = {"extra": "allow"}


class GridComparisonResponse(BaseModel):
    total_points: int
    overall_status: str
    etrs89_to_st70: DirectionComparisonSummary
    st70_to_etrs89: DirectionComparisonSummary
    model_config = {"extra": "allow"}


# Legacy TransdatOnline Statistics Schemas
class LegacyStatsSummary(BaseModel):
    total_points: int
    app_st70: int
    app_st30: int
    service_st70: int
    service_st30: int
    from_date: Optional[str] = None
    to_date: Optional[str] = None
    recorded_days: int
    model_config = {"extra": "allow"}


class LegacyTimelineMetrics(BaseModel):
    app_st70: int = 0
    app_st30: int = 0
    service_st70: int = 0
    service_st30: int = 0
    app_total: int = 0
    service_total: int = 0
    st70_total: int = 0
    st30_total: int = 0
    total: int = 0
    model_config = {"extra": "allow"}


class LegacyTimelineSeriesItem(BaseModel):
    date: str
    points: LegacyTimelineMetrics


class LegacyTimelineSummary(BaseModel):
    total_points: int
    from_date: Optional[str] = None
    to_date: Optional[str] = None
    days_count: Optional[int] = None
    period_count: Optional[int] = None


class LegacyTimelineResponse(BaseModel):
    summary: LegacyTimelineSummary
    series: List[LegacyTimelineSeriesItem]
    bounds: Dict[str, str] = Field(default_factory=dict)
