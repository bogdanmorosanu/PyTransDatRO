"""Legacy TransDatOnline compatibility router.
Reproduces the 10+ year-old Java servlet contract at /transdatonline/cooOpService and /cooOpService.
"""
import json
from typing import Any, Dict, List, Optional
from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import JSONResponse
import pytransdatro.exceptions as exc

router = APIRouter(tags=["Legacy TransDatOnline Compatibility"])


def _transform_legacy_point(trans, coo_op: str, raw_coords: List[float]) -> Dict[str, Any]:
    """Transforms a single coordinate according to legacy rules and warning strings."""
    dim = len(raw_coords)
    if dim < 2:
        return {"coos": raw_coords, "warning": "Wrong dimension"}

    if coo_op == "Stereo70ToETRS89":
        n, e = raw_coords[0], raw_coords[1]
        z = raw_coords[2] if dim >= 3 else None
        try:
            res = trans.st70_to_etrs89(n, e, z=z, source=2)
            return {"coos": list(res)}
        except exc.OutOfGridErr:
            return {"coos": raw_coords, "warning": "Out of grid"}
        except exc.NoDataGridErr:
            return {"coos": raw_coords, "warning": "No data on grid"}

    elif coo_op == "ETRS89ToStereo70":
        lat, lon = raw_coords[0], raw_coords[1]
        h = raw_coords[2] if dim >= 3 else None
        try:
            res = trans.etrs89_to_st70(lat, lon, h=h, source=2)
            return {"coos": list(res)}
        except exc.OutOfGridErr:
            return {"coos": raw_coords, "warning": "Out of grid"}
        except exc.NoDataGridErr:
            return {"coos": raw_coords, "warning": "No data on grid"}

    elif coo_op in ("Stereo30ToETRS89", "ETRS89ToStereo30"):
        return {
            "coos": raw_coords,
            "warning": "Stereo30 transformation is obsolete and unsupported"
        }

    else:
        raise HTTPException(
            status_code=400,
            detail=f"There is no defined transformation with name {coo_op}"
        )


@router.get("/transdatonline/cooOpService")
@router.get("/cooOpService")
async def legacy_coo_op_get(
    request: Request,
    cooOp: Optional[str] = Query(None),
    coos: Optional[str] = Query(None)
):
    """Legacy GET coordinate transformation endpoint."""
    if not cooOp or not coos:
        raise HTTPException(
            status_code=400,
            detail="Invalid request! cooOp or coos variable is missing!"
        )

    # Semicolon-delimited coordinate string
    tokens = coos.split(";")
    parsed_coords: List[float] = []
    for token in tokens:
        stripped = token.strip()
        if not stripped:
            continue
        try:
            parsed_coords.append(float(stripped))
        except ValueError:
            return JSONResponse(content={"coos": [], "warning": "Invalid coordinate data"})

    if len(parsed_coords) < 2:
        return JSONResponse(content={"coos": parsed_coords, "warning": "Wrong dimension"})

    trans = request.app.state.trans
    result = _transform_legacy_point(trans, cooOp, parsed_coords)
    return JSONResponse(content=result)


@router.post("/transdatonline/cooOpService")
@router.post("/cooOpService")
async def legacy_coo_op_post(request: Request):
    """Legacy POST coordinate transformation endpoint with dual form/JSON content-type support."""
    coo_op: Optional[str] = None
    raw_coos_array: Any = None

    content_type = request.headers.get("content-type", "").lower()

    if "application/json" in content_type:
        try:
            body = await request.json()
            if isinstance(body, dict):
                coo_op = body.get("cooOp")
                raw_coos_array = body.get("coosArray")
        except Exception:
            raise HTTPException(status_code=400, detail="Invalid JSON payload")
    else:
        # Form-urlencoded or multipart form data
        form = await request.form()
        coo_op = form.get("cooOp")
        raw_coos_array = form.get("coosArray")

    if not coo_op or raw_coos_array is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid request! cooOp or coosArray variable is missing!"
        )

    # If coosArray arrived as a string, parse it
    if isinstance(raw_coos_array, str):
        try:
            items = json.loads(raw_coos_array)
        except Exception:
            raise HTTPException(status_code=400, detail="Malformed JSON in coosArray")
    elif isinstance(raw_coos_array, list):
        items = raw_coos_array
    else:
        raise HTTPException(status_code=400, detail="coosArray must be a list or JSON string array")

    trans = request.app.state.trans
    results: List[Dict[str, Any]] = []

    for item in items:
        if not isinstance(item, dict) or "coos" not in item:
            results.append({"coos": [], "warning": "Invalid coordinate data"})
            continue

        raw_coords = item["coos"]
        if not isinstance(raw_coords, (list, tuple)):
            results.append({"coos": [], "warning": "Invalid coordinate data"})
            continue

        try:
            numeric_coords = [float(c) for c in raw_coords]
        except (ValueError, TypeError):
            results.append({"coos": [], "warning": "Invalid coordinate data"})
            continue

        res = _transform_legacy_point(trans, coo_op, numeric_coords)
        results.append(res)

    return JSONResponse(content=results)
