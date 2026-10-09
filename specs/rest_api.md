# PyTransDatRO: REST API & Web Service Specification

## 1. Overview & Objectives

This document defines the formal architecture, endpoint contracts, Pydantic data schemas, operational constraints, and deployment specifications for the **PyTransDatRO REST API and Web Service**.

### Primary Objectives:
1. **100% Backward Compatibility**: Expose `/transdatonline/cooOpService` matching the 10+ year old TransDatRO/RomGeo Java service contract byte-for-byte, ensuring existing clients (QGIS, CAD plugins, scripts, external services) continue to operate without disruption.
2. **Standardized Geodetic Units & Modern Units Support**:
   - **Legacy Service**: Strictly **radians** for geographic angular coordinates (ETRS89) and **meters** for planimetric coordinates (Stereo 70) and heights.
   - **Modern API (`/api/v1`)**: Defaults to standard **radians**, with optional first-class support for **decimal degrees** (`unit: "degrees" | "radians"`), allowing modern web GIS clients (Leaflet, MapLibre, Turf.js) to interact natively without manual conversion.
3. **High-Performance Stateless Execution**: Execute single-point transformations in microsecond timescales ($< 0.05\text{ ms}$ per point) using an in-memory `TransRO` singleton.
4. **Massive Bulk File Streaming**: Support streaming delimited text files up to **100,000 points** per request with chunked HTTP transfer encoding (benchmarked at $\approx 50,000\text{ points/sec}$ in pure Python).
5. **Lean Transfer Principle**: Process strictly coordinate numbers, matching input to output by row index ($N$ coordinates in $\rightarrow$ $N$ results out) without transferring redundant metadata.
6. **Non-Blocking Telemetry**: Integrate transparently with the SQLite WAL telemetry engine (`SqliteUsageLogger`), tagging calls with `source=2` (`rest_api`) for daily volume and spatial density tracking.

---

## 2. Invariants & System Constraints

- **Decoupled Architecture (Option 1 - Monorepo Extras)**:
  - The core geodetic library `pytransdatro` remains pure Python with standard library dependencies only.
  - Web service dependencies are defined as optional extras (`pytransdatro[api]`):
    ```toml
    [project.optional-dependencies]
    api = [
        "fastapi>=0.110.0",
        "uvicorn[standard]>=0.28.0",
        "pydantic>=2.6.0",
        "slowapi>=0.1.9",
        "python-multipart>=0.0.9",
    ]
    ```
- **RAM Grid Singleton**: The active `.spg` binary grid is loaded once into memory during application lifespan startup. Each worker process retains an independent instance in memory ($\approx 8\text{ MB}$ RAM per process; $\approx 32\text{ MB}$ for 4 Uvicorn workers). Requests share this read-only grid instance with zero I/O per query.
- **Coordinate Conventions**:
  - **Stereo 70 (Projected)**: Order `[Northing, Easting]` (2D) or `[Northing, Easting, Normal_Height]` (3D), units in **meters**.
  - **ETRS89 (Geographic)**: Order `[Latitude, Longitude]` (2D) or `[Latitude, Longitude, Ellipsoidal_Height]` (3D). Latitude and longitude in **radians** by default (or **decimal degrees** when explicitly requested in `/api/v1`), height in **meters**.
- **Dimension Matching Invariant**:
  - 2D input (`[N, E]` or `[Lat, Lon]`) $\rightarrow$ strictly 2D output.
  - 3D input (`[N, E, H]` or `[Lat, Lon, h]`) $\rightarrow$ strictly 3D output.
- **Index-Based Match**: No arbitrary point IDs or attribute codes are ingested by the service. Row $i$ of the output corresponds to Row $i$ of the input.

---

## 3. Endpoints Matrix

| Route | Method | Purpose | Protocol / Payload Format |
| :--- | :--- | :--- | :--- |
| `/transdatonline/cooOpService`<br/>`/cooOpService` | `GET` | **Legacy Single Point** | Query params: `cooOp`, `coos` (semicolon-delimited numbers) |
| `/transdatonline/cooOpService`<br/>`/cooOpService` | `POST` | **Legacy Batch** | Form fields or JSON: `cooOp`, `coosArray` (dual content-type support) |
| `/api/v1/transform/point` | `POST` | **Modern Single Point** | JSON: `{"op": "Stereo70ToETRS89", "coos": [N, E, H], "unit": "degrees"}` |
| `/api/v1/transform/batch` | `POST` | **Modern In-Memory Batch** | JSON: `{"op": "Stereo70ToETRS89", "points": [[...], [...]], "unit": "degrees"}` |
| `/api/v1/transform/file` | `POST` | **Delimited File Streaming** | Multipart form: `file`, `op`, `unit`, `delimiter` (up to 100,000 points) |
| `/api/v1/grid/info` | `GET` | **Grid Metadata & Bounds** | GeoJSON extent, grid boundaries, SPG version, CRS info |
| `/api/v1/telemetry/stats` | `GET` | **Usage Transparency** | Daily point/call counts, country-wide aggregated density heatmap |
| `/health` | `GET` | **Health & Liveness Probe** | Health status, uptime, grid memory status |

---

## 4. Legacy Compatibility: `/transdatonline/cooOpService`

### 4.1. HTTP GET (Single Point)
* **Query Parameters**:
  * `cooOp` (string, required): `Stereo70ToETRS89`, `ETRS89ToStereo70`, `Stereo30ToETRS89`, `ETRS89ToStereo30`.
  * `coos` (string, required): Semicolon-delimited numbers in **radians/meters**, e.g. `500000;500000;100` or `500000;500000`.
* **Responses (Status Code: 200 OK)**:
  * **Success (3D)**:
    ```json
    {
      "coos": [0.8028465450500997, 0.43630521911977527, 139.60844039916992]
    }
    ```
  * **Success (2D)**:
    ```json
    {
      "coos": [0.8028465450500997, 0.43630521911977527]
    }
    ```
  * **Out of Grid (`OutOfGridErr`)**: Original input coordinates echoed back, with warning:
    ```json
    {
      "coos": [5.0, 1.0, 1.0],
      "warning": "Out of grid"
    }
    ```
  * **No Data on Grid (`NoDataGridErr`)**: Original input coordinates echoed back, with warning:
    ```json
    {
      "coos": [500000.0, 500000.0, 100.0],
      "warning": "No data on grid"
    }
    ```
  * **Invalid Input**: Empty array, with warning:
    ```json
    {
      "coos": [],
      "warning": "Invalid coordinate data"
    }
    ```
  * **Obsolete Stereo 30**: Original input coordinates echoed back, with warning:
    ```json
    {
      "coos": [500000.0, 500000.0],
      "warning": "Stereo30 transformation is obsolete and unsupported"
    }
    ```
  * **Missing Parameter Error (Status: 400 Bad Request)**:
    ```json
    {"detail": "Invalid request! cooOp or coos variable is missing!"}
    ```

---

### 4.2. HTTP POST (Batch Points with Dual Content-Type Support)
The legacy adapter accepts requests from both classical form submissions and modern JSON clients:

1. **Form-Encoded (`application/x-www-form-urlencoded` or `multipart/form-data`)**:
   - `cooOp`: string
   - `coosArray`: stringified JSON array, e.g. `[{"coos": [500000, 500000, 100]}, {"coos": [5, 1, 1]}]`
2. **Direct JSON (`application/json`)**:
   ```json
   {
     "cooOp": "Stereo70ToETRS89",
     "coosArray": [
       {"coos": [500000.0, 500000.0, 100.0]},
       {"coos": [5.0, 1.0, 1.0]}
     ]
   }
   ```
* **Response (Status Code: 200 OK)**:
  ```json
  [
    {"coos": [0.8028465450500997, 0.43630521911977527, 139.60844039916992]},
    {"coos": [5.0, 1.0, 1.0], "warning": "Out of grid"}
  ]
  ```

---

## 5. Modern API Specification (`/api/v1`)

### 5.1. Pydantic Domain Schemas (`api/schemas/`)

#### A. Single Point (`api/schemas/point.py`)
```python
from enum import Enum
from typing import List, Optional
from pydantic import BaseModel, Field

class OperationEnum(str, Enum):
    STEREO70_TO_ETRS89 = "Stereo70ToETRS89"
    ETRS89_TO_STEREO70 = "ETRS89ToStereo70"

class AngleUnitEnum(str, Enum):
    RADIANS = "radians"
    DEGREES = "degrees"

class PointRequest(BaseModel):
    op: OperationEnum
    coos: List[float] = Field(..., min_length=2, max_length=3, description="[N, E[, H]] or [Lat, Lon[, h]]")
    unit: AngleUnitEnum = Field(default=AngleUnitEnum.RADIANS, description="Angular unit for ETRS89 coordinates")

class PointResponse(BaseModel):
    coos: List[float]
    warning: Optional[str] = None
```

#### B. Batch Transformation (`api/schemas/batch.py`)
```python
class BatchRequest(BaseModel):
    op: OperationEnum
    points: List[List[float]] = Field(..., min_length=1, max_length=20000, description="List of 2D or 3D coordinate arrays")
    unit: AngleUnitEnum = Field(default=AngleUnitEnum.RADIANS)

class BatchResponse(BaseModel):
    results: List[PointResponse]
    count: int
```

#### C. Grid Info & Telemetry (`api/schemas/info.py`)
```python
class BoundingBox(BaseModel):
    min_northing: float
    min_easting: float
    max_northing: float
    max_easting: float

class GridInfoResponse(BaseModel):
    grid_file: str
    bounds: BoundingBox
    crs: dict
    quasigeoid_model: str

class DailyStatItem(BaseModel):
    date: str
    source: str
    calls: int
    points: int

class TelemetrySummary(BaseModel):
    total_calls: int
    total_points: int

class TelemetryStatsResponse(BaseModel):
    summary: TelemetrySummary
    daily_trend: List[DailyStatItem]
```

---

### 5.2. Modern Endpoints Behavior

#### A. `POST /api/v1/transform/point`
* **Request**:
  ```json
  {
    "op": "Stereo70ToETRS89",
    "coos": [500000.0, 500000.0, 100.0],
    "unit": "degrees"
  }
  ```
* **Response** (HTTP 200):
  ```json
  {
    "coos": [45.99999999999999, 24.99999999999999, 139.60844039916992],
    "warning": null
  }
  ```

#### B. `POST /api/v1/transform/batch`
* Limit: up to **20,000 points** in-memory.
* **Request**:
  ```json
  {
    "op": "Stereo70ToETRS89",
    "points": [
      [500000.0, 500000.0, 100.0],
      [5.0, 1.0, 1.0]
    ],
    "unit": "degrees"
  }
  ```
* **Response** (HTTP 200):
  ```json
  {
    "results": [
      {"coos": [45.99999999999999, 24.99999999999999, 139.60844039916992], "warning": null},
      {"coos": [5.0, 1.0, 1.0], "warning": "Out of grid"}
    ],
    "count": 2
  }
  ```

---

### 5.3. Bulk File Streaming: `POST /api/v1/transform/file`

* **Content-Type**: `multipart/form-data`
* **Form Parameters**:
  * `file`: Uploaded `.csv` or `.txt` file (max **20 MB**, up to **100,000 coordinate rows**).
  * `op`: `Stereo70ToETRS89` or `ETRS89ToStereo70`.
  * `unit`: `radians` or `degrees` (default: `radians`).
  * `delimiter`: `auto`, `comma`, `semicolon`, `tab`, `space` (default: `auto`).
* **Input File Robustness**:
  * Strips empty lines and whitespace.
  * Ignores header comments starting with `#`.
  * Automatically detects delimiter if set to `auto` from the first non-comment line.
* **Output Streaming (`StreamingResponse`)**:
  * `Content-Type`: `text/csv; charset=utf-8`
  * `Content-Disposition`: `attachment; filename="transformed_coordinates.csv"`
  * Returns transformed rows line-by-line using chunked transfer encoding.
  * If a point is out of grid, the input coordinates are echoed with `# Out of grid` appended at the end of the line:
    ```csv
    46.0000000000,25.0000000000,139.6084
    5.0000000000,1.0000000000,1.0000 # Out of grid
    ```

---

## 6. Metadata, Monitoring & Telemetry Endpoints

### 6.1. `GET /health`
```json
{
  "status": "healthy",
  "grid_loaded": true,
  "engine_version": "1.0.0",
  "active_workers": 1,
  "telemetry_active": true
}
```

### 6.2. `GET /api/v1/grid/info`
```json
{
  "grid_file": "rom_grid3d_25.09.spg",
  "bounds": {
    "min_northing": 213634.564,
    "min_easting": 109783.04,
    "max_northing": 785634.564,
    "max_easting": 890783.04
  },
  "crs": {
    "projected": "EPSG:31700 (Pulkovo 1942(58) / Stereo70)",
    "geographic": "EPSG:4258 (ETRS89)"
  },
  "quasigeoid_model": "Modern Romanian Quasigeoid (1D nearest-neighbor)"
}
```

### 6.3. `GET /api/v1/telemetry/stats`
Queries SQLite `logs/usage.db` in read-only WAL mode to return high-level usage metrics:
```json
{
  "summary": {
    "total_calls": 14250,
    "total_points": 892300
  },
  "daily_trend": [
    {"date": "2026-10-08", "source": "rest_api", "calls": 420, "points": 18500},
    {"date": "2026-10-09", "source": "rest_api", "calls": 510, "points": 24200}
  ]
}
```

---

## 7. Security, Rate Limiting & Abuse Prevention

1. **Fair-Use Rate Limiting (`slowapi`)**:
   - Single-point endpoints (`/point`, `/transdatonline/cooOpService` GET): **60 requests / minute / IP**.
   - Batch and file streaming endpoints (`/batch`, `/file`, `/transdatonline/cooOpService` POST): **10 requests / minute / IP**.
2. **Payload Size Clamps**:
   - In-memory JSON batch (`/batch`): Maximum **20,000 points**.
   - Streaming file upload (`/file`): Maximum **100,000 rows** and **20 MB** file size.
3. **CORS Policy**:
   - Open CORS enabled (`Allow-Origin: *`, `Allow-Methods: GET, POST, OPTIONS`, `Allow-Headers: *`) to ensure compatibility with client-side mapping tools, QGIS plugins, and external web applications.

---

## 8. Directory Layout & Package Structure

```text
PyTransDatRO/
├── pytransdatro/              # Pure Python Geodetic Engine (Unmodified)
│   ├── trans_ro.py
│   ├── exceptions.py
│   ├── telemetry.py
│   └── grids/
│
├── api/                       # REST API Application Package
│   ├── __init__.py
│   ├── main.py                # FastAPI factory, lifespan, CORS & rate-limiter setup
│   ├── dependencies.py        # App state accessors (TransRO singleton, rate limiter)
│   ├── schemas/
│   │   ├── __init__.py
│   │   ├── point.py           # PointRequest, PointResponse, enums
│   │   ├── batch.py           # BatchRequest, BatchResponse
│   │   └── info.py            # GridInfoResponse, TelemetryStatsResponse
│   └── routes/
│       ├── __init__.py
│       ├── legacy.py          # /transdatonline/cooOpService adapter (GET & POST)
│       ├── point.py           # /api/v1/transform/point
│       ├── batch.py           # /api/v1/transform/batch
│       ├── file.py            # /api/v1/transform/file (streaming generator)
│       ├── info.py            # /api/v1/grid/info & /api/v1/telemetry/stats
│       └── health.py          # /health
│
├── tests/
│   ├── test_legacy_api_compatibility.py # Full parity tests matching docs/transdatonline
│   ├── test_api_v1_endpoints.py         # Tests for modern point, batch, and file streaming
│   └── test_api_rate_limiting.py        # Slowapi rate limit enforcement tests
│
└── pyproject.toml             # Defines [project.optional-dependencies] api = [...]
```

---

## 9. Launch & Operational Commands

* **Local Development**:
  ```bash
  uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
  ```
* **Production Multi-Worker**:
  ```bash
  gunicorn api.main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000
  ```
* **OpenAPI Documentation**:
  - Swagger UI: `http://localhost:8000/docs`
  - ReDoc: `http://localhost:8000/redoc`
