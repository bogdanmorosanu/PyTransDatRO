"""Health and liveness probe endpoint."""

from typing import Any, Dict

from fastapi import APIRouter, Request

from pytransdatro import __version__ as pytransdatro_version

router = APIRouter(tags=["Health & Monitoring"])


@router.get("/health")
async def health_check(request: Request) -> Dict[str, Any]:
    """Returns application health, grid status, and telemetry diagnostics."""
    has_trans = hasattr(request.app.state, "trans") and request.app.state.trans is not None
    telemetry = getattr(request.app.state, "telemetry", None)

    response: Dict[str, Any] = {
        "status": "healthy" if has_trans else "degraded",
        "grid_loaded": has_trans,
        "engine_version": pytransdatro_version,
        "telemetry_active": False,
        "telemetry_status": "disabled",
    }

    if telemetry is not None:
        try:
            diag = telemetry.get_diagnostics()
            is_alive = diag.get("is_alive", False)
            last_err = diag.get("last_error")

            if not is_alive or last_err is not None:
                response["status"] = "degraded"
                response["telemetry_status"] = "degraded"
                response["telemetry_active"] = is_alive
                response["telemetry_error"] = last_err or "Worker thread is not running"
            else:
                response["telemetry_status"] = "healthy"
                response["telemetry_active"] = True
        except Exception as exc:
            response["status"] = "degraded"
            response["telemetry_status"] = "degraded"
            response["telemetry_error"] = str(exc)

    return response
