"""Modern in-memory batch coordinate transformation endpoint."""

from typing import List

from fastapi import APIRouter, Depends, HTTPException

from api.config import Settings, get_settings
from api.dependencies import get_client_source, get_trans_engine
from api.schemas.batch import BatchRequest, BatchResponse
from api.schemas.point import PointResponse
from api.services.transformer import execute_batch_transformation, format_modern_response
from pytransdatro.trans_ro import TransRO

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])


@router.post("/batch", response_model=BatchResponse)
async def transform_batch(
    req: BatchRequest,
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
    settings: Settings = Depends(get_settings),
) -> BatchResponse:
    """Batch processes a list of coordinates in memory (up to configured max points).

    Provides per-point failure isolation (e.g., points out of grid are marked with warnings
    without failing valid coordinates in the same batch).
    """
    if len(req.points) > settings.max_batch_points:
        raise HTTPException(
            status_code=422,
            detail=f"Batch size exceeds maximum limit of {settings.max_batch_points} points.",
        )

    batch_tuples = execute_batch_transformation(
        trans, req.op.value, req.points, unit=req.unit, source=source
    )
    results: List[PointResponse] = [
        format_modern_response(out_coos, warn) for out_coos, warn in batch_tuples
    ]

    return BatchResponse(results=results, count=len(results))
