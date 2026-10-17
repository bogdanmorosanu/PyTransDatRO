"""Modern single point coordinate transformation endpoint."""

from fastapi import APIRouter, Depends

from api.dependencies import get_client_source, get_trans_engine
from api.schemas.point import PointRequest, PointResponse
from api.services.transformer import execute_point_transformation, format_modern_response
from pytransdatro.trans_ro import TransRO

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])


@router.post("/point", response_model=PointResponse)
async def transform_point(
    req: PointRequest,
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
) -> PointResponse:
    """Transforms a single coordinate between Stereo 70 and ETRS89.

    Supports both radians (default geodetic standard) and decimal degrees.
    """
    out_coos, warn = execute_point_transformation(
        trans, req.op.value, req.coos, unit=req.unit, source=source
    )
    return format_modern_response(out_coos, warn)
