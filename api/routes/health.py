"""Health and liveness probe endpoint."""

from typing import Optional

from fastapi import APIRouter, Request
from pydantic import BaseModel, Field

from api.config import get_settings
from pytransdatro import __version__ as pytransdatro_version

router = APIRouter(tags=["Health & Monitoring"])


class RuntimeLimitsModel(BaseModel):
    max_batch_points: int = Field(..., description="Maximum points in an in-memory batch")
    max_file_points: int = Field(..., description="Maximum points in a streamed file")
    max_file_upload_bytes: int = Field(..., description="Maximum upload payload bytes")


class HealthResponse(BaseModel):
    status: str = Field(..., description="Overall health status: 'healthy' or 'degraded'")
    grid_loaded: bool = Field(..., description="Whether SPG binary grid is loaded in memory")
    engine_version: str = Field(..., description="PyTransDatRO geodetic engine semantic version")
    telemetry_active: bool = Field(..., description="Whether telemetry worker thread is active")
    telemetry_status: str = Field(
        ..., description="Telemetry status: 'healthy', 'degraded', or 'disabled'"
    )
    telemetry_error: Optional[str] = Field(
        default=None, description="Diagnostics error message if degraded"
    )
    limits: RuntimeLimitsModel = Field(
        ..., description="Active server transformation volume and size limits"
    )


@router.get("/health", response_model=HealthResponse)
async def health_check(request: Request) -> HealthResponse:
    """Returns application health, grid status, active limits, and telemetry diagnostics."""
    has_trans = hasattr(request.app.state, "trans") and request.app.state.trans is not None
    telemetry = getattr(request.app.state, "telemetry", None)
    settings = getattr(request.app.state, "settings", None) or get_settings()

    status = "healthy" if has_trans else "degraded"
    telemetry_active = False
    telemetry_status = "disabled"
    telemetry_error: Optional[str] = None

    if telemetry is not None:
        try:
            diag = telemetry.get_diagnostics()
            is_alive = diag.get("is_alive", False)
            last_err = diag.get("last_error")

            if not is_alive or last_err is not None:
                status = "degraded"
                telemetry_status = "degraded"
                telemetry_active = is_alive
                telemetry_error = last_err or "Worker thread is not running"
            else:
                telemetry_status = "healthy"
                telemetry_active = True
        except Exception as exc:
            status = "degraded"
            telemetry_status = "degraded"
            telemetry_error = str(exc)

    return HealthResponse(
        status=status,
        grid_loaded=has_trans,
        engine_version=pytransdatro_version,
        telemetry_active=telemetry_active,
        telemetry_status=telemetry_status,
        telemetry_error=telemetry_error,
        limits=RuntimeLimitsModel(
            max_batch_points=settings.max_batch_points,
            max_file_points=settings.max_file_points,
            max_file_upload_bytes=settings.max_file_upload_bytes,
        ),
    )
