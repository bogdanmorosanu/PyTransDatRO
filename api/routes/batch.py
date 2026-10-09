"""Modern in-memory batch coordinate transformation endpoint.
"""
import math
from typing import List
from fastapi import APIRouter, Depends
from api.dependencies import get_trans_engine
from api.schemas.batch import BatchRequest, BatchResponse
from api.schemas.point import OperationEnum, AngleUnitEnum, PointResponse
from pytransdatro.trans_ro import TransRO
import pytransdatro.exceptions as exc

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])


@router.post("/batch", response_model=BatchResponse)
async def transform_batch(
    req: BatchRequest,
    trans: TransRO = Depends(get_trans_engine)
) -> BatchResponse:
    """Batch processes a list of coordinates in memory (up to 20,000 points).
    
    Provides per-point failure isolation (e.g., points out of grid are marked with warnings
    without failing valid coordinates in the same batch).
    """
    results: List[PointResponse] = []
    is_degrees = (req.unit == AngleUnitEnum.DEGREES)

    for pt in req.points:
        dim = len(pt)
        if dim < 2:
            results.append(PointResponse(coos=pt, warning="Wrong dimension"))
            continue

        try:
            if req.op == OperationEnum.STEREO70_TO_ETRS89:
                n, e = pt[0], pt[1]
                z = pt[2] if dim >= 3 else None
                res = trans.st70_to_etrs89(n, e, z=z, source=2)

                if is_degrees:
                    lat_deg = math.degrees(res[0])
                    lon_deg = math.degrees(res[1])
                    out_coords = [lat_deg, lon_deg]
                    if dim >= 3:
                        out_coords.append(res[2])
                    results.append(PointResponse(coos=out_coords, warning=None))
                else:
                    results.append(PointResponse(coos=list(res), warning=None))

            elif req.op == OperationEnum.ETRS89_TO_STEREO70:
                lat, lon = pt[0], pt[1]
                if is_degrees:
                    lat = math.radians(lat)
                    lon = math.radians(lon)

                h = pt[2] if dim >= 3 else None
                res = trans.etrs89_to_st70(lat, lon, h=h, source=2)
                results.append(PointResponse(coos=list(res), warning=None))

        except exc.OutOfGridErr:
            results.append(PointResponse(coos=pt, warning="Out of grid"))
        except exc.NoDataGridErr:
            results.append(PointResponse(coos=pt, warning="No data on grid"))

    return BatchResponse(results=results, count=len(results))
