"""Modern single point coordinate transformation endpoint.
"""
import math
from fastapi import APIRouter, Depends
from api.dependencies import get_trans_engine
from api.schemas.point import PointRequest, PointResponse, OperationEnum, AngleUnitEnum
from pytransdatro.trans_ro import TransRO
import pytransdatro.exceptions as exc

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])


@router.post("/point", response_model=PointResponse)
async def transform_point(
    req: PointRequest,
    trans: TransRO = Depends(get_trans_engine)
) -> PointResponse:
    """Transforms a single coordinate between Stereo 70 and ETRS89.
    
    Supports both radians (default geodetic standard) and decimal degrees.
    """
    raw_coords = req.coos
    dim = len(raw_coords)

    try:
        if req.op == OperationEnum.STEREO70_TO_ETRS89:
            n, e = raw_coords[0], raw_coords[1]
            z = raw_coords[2] if dim >= 3 else None
            res = trans.st70_to_etrs89(n, e, z=z, source=2)

            if req.unit == AngleUnitEnum.DEGREES:
                lat_deg = math.degrees(res[0])
                lon_deg = math.degrees(res[1])
                out_coords = [lat_deg, lon_deg]
                if dim >= 3:
                    out_coords.append(res[2])
                return PointResponse(coos=out_coords, warning=None)

            return PointResponse(coos=list(res), warning=None)

        elif req.op == OperationEnum.ETRS89_TO_STEREO70:
            lat, lon = raw_coords[0], raw_coords[1]
            if req.unit == AngleUnitEnum.DEGREES:
                lat = math.radians(lat)
                lon = math.radians(lon)

            h = raw_coords[2] if dim >= 3 else None
            res = trans.etrs89_to_st70(lat, lon, h=h, source=2)
            return PointResponse(coos=list(res), warning=None)

    except exc.OutOfGridErr:
        return PointResponse(coos=raw_coords, warning="Out of grid")
    except exc.NoDataGridErr:
        return PointResponse(coos=raw_coords, warning="No data on grid")
