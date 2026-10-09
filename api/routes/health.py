"""Health and liveness probe endpoint.
"""
from typing import Dict, Any
from fastapi import APIRouter, Request

router = APIRouter(tags=["Health & Monitoring"])


@router.get("/health")
async def health_check(request: Request) -> Dict[str, Any]:
    """Returns application health and grid status."""
    has_trans = hasattr(request.app.state, "trans") and request.app.state.trans is not None
    has_telemetry = hasattr(request.app.state, "telemetry") and request.app.state.telemetry is not None

    return {
        "status": "healthy" if has_trans else "degraded",
        "grid_loaded": has_trans,
        "engine_version": "1.0.0",
        "telemetry_active": has_telemetry
    }
