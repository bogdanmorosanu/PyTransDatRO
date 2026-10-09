# FastAPI Backward Compatibility Adapter Plan

## 1. Goal

Implement a drop-in compatibility router within the upcoming FastAPI application (`api/routes/legacy.py`) that strictly adheres to [`docs/transdatonline/02_legacy_service_api_specification.md`](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/02_legacy_service_api_specification.md).

This guarantees that existing legacy consumers calling `https://www.geo-spatial.org/transdatonline/cooOpService` will experience zero regressions or breaking changes when routed to the new Python backend.

---

## 2. Architecture & File Placement

```text
PyTransDatRO/
├── pytransdatro/              # Pure Python Geodetic Engine (UNMODIFIED)
│   ├── trans_ro.py
│   └── exceptions.py
│
├── api/                       # REST API Application
│   ├── main.py                # FastAPI entrypoint
│   ├── dependencies.py        # Lifespan & TransRO singleton injection
│   └── routes/
│       ├── legacy.py          # <--- THE COMPATIBILITY ADAPTER
│       ├── point.py           # Modern /api/v1/transform/point
│       ├── batch.py           # Modern /api/v1/transform/batch
│       └── file.py            # Modern /api/v1/transform/file
```

---

## 3. Router Implementation Blueprint (`api/routes/legacy.py`)

### A. Routes Definition
Mount the router with aliases for both historical root and application paths:
* `@router.get("/transdatonline/cooOpService")`
* `@router.get("/cooOpService")`
* `@router.post("/transdatonline/cooOpService")`
* `@router.post("/cooOpService")`

### B. Python Adapter Logic (Conceptual Skeleton)

```python
import json
from typing import Optional, List, Dict, Any
from fastapi import APIRouter, Request, Query, Form, HTTPException
from fastapi.responses import JSONResponse
import pytransdatro.exceptions as exc

router = APIRouter(tags=["Legacy TransDatOnline Compatibility"])

def transform_single_point(trans, coo_op: str, raw_coords: List[float]) -> Dict[str, Any]:
    dim = len(raw_coords)
    if dim < 2:
        return {"coos": raw_coords, "warning": "Wrong dimension"}
    
    try:
        if coo_op == "Stereo70ToETRS89":
            n, e = raw_coords[0], raw_coords[1]
            z = raw_coords[2] if dim >= 3 else None
            res = trans.st70_to_etrs89(n, e, z=z, source=2) # source=2 for rest_api
            return {"coos": list(res)}

        elif coo_op == "ETRS89ToStereo70":
            lat, lon = raw_coords[0], raw_coords[1]
            h = raw_coords[2] if dim >= 3 else None
            res = trans.etrs89_to_st70(lat, lon, h=h, source=2)
            return {"coos": list(res)}

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

    except exc.OutOfGridErr:
        return {"coos": raw_coords, "warning": "Out of grid"}
    except exc.NoDataGridErr:
        return {"coos": raw_coords, "warning": "No data on grid"}


@router.get("/transdatonline/cooOpService")
@router.get("/cooOpService")
async def legacy_coo_op_get(
    request: Request,
    cooOp: Optional[str] = Query(None),
    coos: Optional[str] = Query(None)
):
    if not cooOp or not coos:
        raise HTTPException(
            status_code=400, 
            detail="Invalid request! cooOp or coos variable is missing!"
        )

    # Parse semicolon-separated coordinates
    tokens = coos.split(";")
    parsed_coords = []
    for token in tokens:
        stripped = token.strip()
        if not stripped:
            continue
        try:
            parsed_coords.append(float(stripped))
        except ValueError:
            return JSONResponse(content={"coos": [], "warning": "Invalid coordinate data"})

    trans = request.app.state.trans
    result = transform_single_point(trans, cooOp, parsed_coords)
    return JSONResponse(content=result)


@router.post("/transdatonline/cooOpService")
@router.post("/cooOpService")
async def legacy_coo_op_post(
    request: Request,
    cooOp: Optional[str] = Form(None),
    coosArray: Optional[str] = Form(None)
):
    if not cooOp or not coosArray:
        raise HTTPException(
            status_code=400, 
            detail="Invalid request! cooOp or coosArray variable is missing!"
        )

    try:
        items = json.loads(coosArray)
    except Exception:
        raise HTTPException(status_code=400, detail="Malformed JSON in coosArray")

    trans = request.app.state.trans
    results = []

    for item in items:
        raw_coords = item.get("coos", [])
        # Check if coordinates contain non-numeric types
        try:
            numeric_coords = [float(c) for c in raw_coords]
        except (ValueError, TypeError):
            results.append({"coos": [], "warning": "Invalid coordinate data"})
            continue

        res = transform_single_point(trans, cooOp, numeric_coords)
        results.append(res)

    return JSONResponse(content=results)
```

---

## 4. Telemetry Integration

* In `transdatonline/server/JSONCooOpService.java`, the legacy app recorded usage into MySQL table `statistics`:
  * `cooOpServiceSt70`
  * `cooOpServiceSt30`
* In the new adapter, every call passes `source=2` (`rest_api`) into `trans.st70_to_etrs89()` and `trans.etrs89_to_st70()`.
* The background `SqliteUsageLogger` in `pytransdatro.telemetry` will automatically:
  1. Increment `daily_stats` counters (`count_calls`, `count_points`) for `source=2`.
  2. Aggregate non-sensitive ~1 km spatial cells in `usage`.
  3. Require zero disk I/O on the request worker thread.

---

## 5. Verification & Testing Strategy

A dedicated test suite (`tests/test_legacy_api_compatibility.py`) will validate:
1. **GET Standard 3D**: `cooOp=Stereo70ToETRS89&coos=500000;500000;100` matches planar coords to $< 10^{-15}\text{ rad}$.
2. **GET Standard 2D**: `cooOp=Stereo70ToETRS89&coos=500000;500000` returns 2-element array.
3. **GET Out-of-Grid**: `coos=5;1;1` returns `{"coos":[5,1,1],"warning":"Out of grid"}`.
4. **GET Invalid Data**: `coos=a;b;c` returns `{"coos":[],"warning":"Invalid coordinate data"}`.
5. **POST Batch**: Multiple items with mixed valid and out-of-grid coordinates.
6. **Error handling**: Missing `cooOp` or `coos` parameters.
