"""Legacy TransDatOnline compatibility router.

Reproduces the 10+ year-old Java servlet contract at /transdatonline/cooOpService.
"""

import json
from typing import Any, Dict, List

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import JSONResponse

from api.config import Settings, get_settings
from api.dependencies import get_client_source, get_trans_engine
from api.schemas.point import AngleUnitEnum
from api.services.transformer import (
    execute_batch_transformation,
    execute_point_transformation,
    format_legacy_dict,
    parse_legacy_semicolon_coos,
    validate_legacy_operation,
)
from pytransdatro.trans_ro import TransRO

router = APIRouter(tags=["Legacy TransDatOnline Compatibility"])


@router.get("/transdatonline/cooOpService")
async def legacy_coo_op_get(
    request: Request,
    cooOp: str = Query(
        ...,
        description="Legacy transformation operation (Stereo70ToETRS89, ETRS89ToStereo70)",
        examples=["Stereo70ToETRS89"],
    ),
    coos: str = Query(
        ...,
        description="Semicolon-delimited coordinate string (e.g. '500000;500000;100')",
        examples=["500000;500000;100"],
    ),
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
):
    """Legacy GET coordinate transformation endpoint."""
    if cooOp is None or coos is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid request! cooOp or coos variable is missing!",
        )

    try:
        validate_legacy_operation(cooOp)
    except ValueError as err:
        raise HTTPException(status_code=400, detail=str(err)) from err

    parsed_coords, parse_warn = parse_legacy_semicolon_coos(coos)
    if parse_warn is not None:
        return JSONResponse(content=format_legacy_dict(parsed_coords, parse_warn))

    out_coos, warn = execute_point_transformation(
        trans, cooOp, parsed_coords, unit=AngleUnitEnum.RADIANS, source=source
    )
    return JSONResponse(content=format_legacy_dict(out_coos, warn))


@router.post("/transdatonline/cooOpService")
async def legacy_coo_op_post(
    request: Request,
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
    settings: Settings = Depends(get_settings),
):
    """Legacy POST coordinate transformation endpoint with dual form/JSON content-type support."""
    coo_op: Any = None
    raw_coos_array: Any = None

    content_type = request.headers.get("content-type", "").lower()

    if "application/json" in content_type:
        try:
            body = await request.json()
            if isinstance(body, dict):
                coo_op = body.get("cooOp")
                raw_coos_array = body.get("coosArray")
        except (json.JSONDecodeError, UnicodeDecodeError) as err:
            raise HTTPException(status_code=400, detail="Invalid JSON payload") from err
    else:
        # Form-urlencoded or multipart form data
        form = await request.form()
        coo_op = form.get("cooOp")
        raw_coos_array = form.get("coosArray")

    if coo_op is None or raw_coos_array is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid request! cooOp or coosArray variable is missing!",
        )

    if not isinstance(coo_op, str):
        raise HTTPException(status_code=400, detail="Invalid cooOp parameter type")

    try:
        validate_legacy_operation(coo_op)
    except ValueError as err:
        raise HTTPException(status_code=400, detail=str(err)) from err

    # If coosArray arrived as a string, parse it
    if isinstance(raw_coos_array, str):
        try:
            items = json.loads(raw_coos_array)
        except (json.JSONDecodeError, TypeError) as err:
            raise HTTPException(status_code=400, detail="Malformed JSON in coosArray") from err
    elif isinstance(raw_coos_array, list):
        items = raw_coos_array
    else:
        raise HTTPException(status_code=400, detail="coosArray must be a list or JSON string array")

    if not isinstance(items, list):
        raise HTTPException(status_code=400, detail="coosArray must decode to a JSON array")

    if len(items) > settings.max_batch_points:
        raise HTTPException(
            status_code=422,
            detail=f"Batch size exceeds maximum limit of {settings.max_batch_points} points.",
        )

    raw_points: List[Any] = [
        item["coos"] if (isinstance(item, dict) and "coos" in item) else None for item in items
    ]
    batch_tuples = execute_batch_transformation(
        trans, coo_op, raw_points, unit=AngleUnitEnum.RADIANS, source=source
    )
    results: List[Dict[str, Any]] = [
        format_legacy_dict(out_coos, warn) for out_coos, warn in batch_tuples
    ]

    return JSONResponse(content=results)
