# Multi-Layered Coordinate Input Validation & Legacy Compatibility Hardening

## Objectives
- [x] Establish a strict multi-layered validation architecture separating **semantic/domain validation** in `pytransdatro/` from **syntactic/wire-format validation** and **exception-to-warning translation** in `api/`.
- [x] Add core domain exceptions (`InvalidCoordinateErr`, `WrongDimensionErr`, `OutOfGrid2DErr`, `OutOfGrid1DErr`) in [`pytransdatro/exceptions.py`](../../pytransdatro/exceptions.py) and enforce finite numeric checks, non-empty sequences, and matching sequence lengths in [`TransRO`](../../pytransdatro/trans_ro.py).
- [x] Implement **6-Point Extended Dual-CRS Bounding Boxes** (`grid2d_bounds_st70`, `grid2d_bounds_etrs89`, `grid1d_bounds_etrs89`, `grid1d_bounds_st70`) in [`TransRO._init_dual_crs_extents`](../../pytransdatro/trans_ro.py) using the 4 grid corners plus the 2 North/South intersections with the Stereo70 projection center ($\varphi_0 = 46^\circ\text{ N}, \lambda_0 = 25^\circ\text{ E} \leftrightarrow E_0 = 500,000\text{ m}$), preventing false-positive rejections on the outward-bowing North parallel while retaining exact grid checks (`covered_by_grid`) in [`trans_grid.py`](../../pytransdatro/trans_grid.py) for points inside curved border wedges.
- [x] Implement a **4-Stage Fault-Isolated Batch Pipeline** (`TransRO.transform_batch`) in [`pytransdatro/trans_ro.py`](../../pytransdatro/trans_ro.py) that isolates invalid/out-of-grid/no-data points per index without aborting valid points and records telemetry in a single batch call.
- [x] Upgrade [`pytransdatro/utils.py`](../../pytransdatro/utils.py) with strict sexagesimal DMS validation ($0 \le m < 60$, $0 \le s < 60$, cardinal sign parsing `N/E/S/W/V`, rejection of conflicting prefix/suffix signs and fractional degrees/minutes when sub-units exist) and centralized DMS formatting (`format_rad_to_dms`).
- [x] Enforce explicit bounds (`OutOfGrid1DErr` / `OutOfGrid2DErr`), pole singularity guards ($\varphi = \pm 90^\circ \to$ `InvalidCoordinateErr`), and No-Data (`NoDataGridErr`) checks in [`trans_grid.py`](../../pytransdatro/trans_grid.py) and [`proj_stereo.py`](../../pytransdatro/proj_stereo.py).
- [x] Centralize single-point and batch API transformation dispatch (`execute_point_transformation`, `execute_batch_transformation`), unit conversion, and warning mapping in [`api/services/transformer.py`](../../api/services/transformer.py).
- [x] Ensure 100% Legacy `TransDatOnline` compatibility in [`api/routes/legacy.py`](../../api/routes/legacy.py) (Java `split(";")` empty token semantics, Java `getParameter()` non-null empty string `is None` semantics, eager `cooOp` validation on empty `coosArray`, strict 2D/3D dimension bounds, and `max_batch_points` enforcement).
- [x] Harden Modern `/api/v1/transform/*` endpoints ([`point.py`](../../api/routes/point.py), [`batch.py`](../../api/routes/batch.py), [`file.py`](../../api/routes/file.py)) with Option B per-point batch isolation (`points: List[Any]`), balanced CSV quote stripping (`_strip_balanced_quotes`), empty intermediate CSV column rejection, strict delimiter validation, file upload size enforcement (`HTTP 413`), and 5+ column rejection.
- [x] Add offline SPG Metadata vs. Actual Grid Array Boundary verification in [`tools/extract_spg_metadata.py`](../../tools/extract_spg_metadata.py) and [`tools/extract_spg_grid_csv.py`](../../tools/extract_spg_grid_csv.py), documenting findings in [`reports/spg_grid_metadata.txt`](../../reports/spg_grid_metadata.txt) and [`reports/spg_grid_extraction_summary.txt`](../../reports/spg_grid_extraction_summary.txt).

---

## Architectural Design & Validation Matrix

### 1. Layer Responsibilities
| Layer | Module(s) | Validation Responsibility | Failure Mechanism |
| :--- | :--- | :--- | :--- |
| **Core Geodetic Engine** | [`pytransdatro/exceptions.py`](../../pytransdatro/exceptions.py)<br>[`pytransdatro/trans_ro.py`](../../pytransdatro/trans_ro.py)<br>[`pytransdatro/trans_grid.py`](../../pytransdatro/trans_grid.py)<br>[`pytransdatro/proj_stereo.py`](../../pytransdatro/proj_stereo.py)<br>[`pytransdatro/utils.py`](../../pytransdatro/utils.py) | Finite float validation (`math.isfinite`), scalar/sequence parity, sequence length equality, pole singularity guard ($\|\varphi\| \ge \pi/2$), strict sexagesimal DMS grammar, 6-point extended dual-CRS bounding box pre-check, exact 2D/1D grid coverage (`OutOfGrid2DErr`, `OutOfGrid1DErr`), No-Data checks (`NoDataGridErr`), and 4-stage fault-isolated batch execution (`transform_batch`). | Raises `InvalidCoordinateErr`, `WrongDimensionErr`, `OutOfGrid2DErr`, `OutOfGrid1DErr`, or `NoDataGridErr` (or returns `(coords, error)` tuples in `transform_batch`). |
| **API Service Adapter** | [`api/services/transformer.py`](../../api/services/transformer.py) | Coordinate tuple dimensionality ($2 \le \text{dim} \le 3$), Legacy `split(";")` tokenization, radian $\leftrightarrow$ degree conversions, and single/batch domain exception translation (`execute_point_transformation`, `execute_batch_transformation`). | Maps domain exceptions to `(coos, warning)` tuples for Legacy and Modern contracts. |
| **API HTTP Routers & Schemas** | [`api/schemas/point.py`](../../api/schemas/point.py)<br>[`api/schemas/batch.py`](../../api/schemas/batch.py)<br>[`api/routes/legacy.py`](../../api/routes/legacy.py)<br>[`api/routes/point.py`](../../api/routes/point.py)<br>[`api/routes/batch.py`](../../api/routes/batch.py)<br>[`api/routes/file.py`](../../api/routes/file.py) | HTTP query/body syntax (`is None` parameter checks), eager `cooOp` verification, `AngleUnitEnum` vs `FileAngleUnitEnum` separation, Option B per-point batch tolerance (`List[Any]`), `max_batch_points` and `max_file_upload_bytes` limits, balanced CSV quote stripping, and file delimiter validation. | Fast-fails structural errors with HTTP `400` / `413` / `422`; isolates point/row domain warnings with HTTP `200`. |

### 2. Endpoint Contract & Warning Matrix
| Scenario | Legacy `/cooOpService` (`GET` & `POST`) | Modern `/api/v1/transform/point` & `/batch` | Modern `/api/v1/transform/file` |
| :--- | :--- | :--- | :--- |
| **Valid 2D / 3D Point** | `{"coos": [...]}` *(warning omitted)* | `{"coos": [...], "warning": null}` | Transformed row (`[id,]c1,c2[,c3]`) |
| **1D ($<2$) or $\ge 4$D ($>3$) Point** | `{"coos": [...], "warning": "Wrong dimension"}` *(input echoed)* | `{"coos": [...], "warning": "Wrong dimension"}` *(input echoed)* | `<raw_line> # Wrong dimension` (for $<2$ or $>4$ tokens) |
| **Non-numeric / `NaN` / `Inf` / `null` / Bad DMS** | `{"coos": [], "warning": "Invalid coordinate data"}` | `{"coos": [], "warning": "Invalid coordinate data"}` | `<raw_line> # Invalid coordinate data` |
| **Empty Semicolon / CSV Token (`500000;;500000` or `P1,500000,,100`)** | `{"coos": [], "warning": "Invalid coordinate data"}` (`500000;500000;` allowed) | N/A | `<raw_line> # Invalid coordinate data` |
| **Outside 2D or 1D Grid Bounds (`OutOfGrid2DErr` / `OutOfGrid1DErr`)** | `{"coos": [...], "warning": "Out of grid"}` *(input echoed)* | `{"coos": [...], "warning": "Out of grid"}` *(input echoed)* | `<raw_line> # Out of grid` |
| **No-Data Grid Cell (`NaN` / `999`)** | `{"coos": [...], "warning": "No data on grid"}` *(input echoed)* | `{"coos": [...], "warning": "No data on grid"}` *(input echoed)* | `<raw_line> # No data on grid` |
| **Obsolete `Stereo30*` Operation** | `{"coos": [...], "warning": "Stereo30 transformation is obsolete and unsupported"}` | HTTP `422` (`VALIDATION_ERROR`) | HTTP `422` (`VALIDATION_ERROR`) |
| **Unknown `cooOp` / `op`** | HTTP `400` (`There is no defined transformation with name ...`) | HTTP `422` (`VALIDATION_ERROR`) | HTTP `422` (`VALIDATION_ERROR`) |

---

## Implementation Details

1. **6-Point Extended Dual-CRS Bounding Boxes & Grid Exceptions ([`pytransdatro/exceptions.py`](../../pytransdatro/exceptions.py), [`pytransdatro/trans_ro.py`](../../pytransdatro/trans_ro.py), [`pytransdatro/trans_grid.py`](../../pytransdatro/trans_grid.py))**:
   - Added `OutOfGrid2DErr(OutOfGridErr)` and `OutOfGrid1DErr(OutOfGridErr)` subclasses so callers can distinguish horizontal (`geodetic_shifts`) vs. vertical (`geoid_heights`) boundary violations while preserving `except OutOfGridErr:` compatibility.
   - Because Oblique Stereographic (`Stereo70`) parallels bow outward relative to a 4-corner geographic rectangle (reaching $\varphi_{\text{mid,N}} \approx 48.5458^\circ > \varphi_{\text{corner,N}} \approx 48.4620^\circ$ at $\lambda_0 = 25^\circ\text{ E}$ for `Grid2D`, and $N_{\text{mid,S}} \approx 203,319.9\text{ m} < N_{\text{corner,S}} \approx 212,213.9\text{ m}$ for `Grid1D`), `TransRO._init_dual_crs_extents()` samples 6 points per grid (4 corners + North/South central meridian intersections) to form a true outer bounding envelope in both CRSs.
   - Pre-checks (`_precheck_st70_point`, `_precheck_etrs89_point`) run immediately at the start of `st70_to_etrs89`, `etrs89_to_st70`, and `transform_batch` (also blocking $2\pi$ periodic radian aliasing and pole singularities before projection), while `Grid2D.covered_by_grid` and `Grid1D.covered_by_grid` catch points falling inside the curved border wedges.

2. **4-Stage Fault-Isolated Batch Pipeline ([`pytransdatro/trans_ro.py`](../../pytransdatro/trans_ro.py), [`api/services/transformer.py`](../../api/services/transformer.py))**:
   - Implemented `TransRO.transform_batch()` to process batches in 4 stages: (1) per-point syntax, dimensionality, and 6-point extended envelope pre-check; (2) batch Stereographic projection and 2D Helmert transformation on survivors; (3) per-point fault-isolated `Grid2D` and `Grid1D` interpolation; and (4) single-call batch telemetry recording (`_collector.record`) and index-order result merging.

3. **Offline SPG Metadata vs. Array Boundary Audit ([`tools/extract_spg_metadata.py`](../../tools/extract_spg_metadata.py), [`tools/extract_spg_grid_csv.py`](../../tools/extract_spg_grid_csv.py))**:
   - Added Section 5 (`Metadata vs. Actual Grid Boundary Integrity & Dual-CRS Extents`) to both offline tools and regenerated [`reports/spg_grid_metadata.txt`](../../reports/spg_grid_metadata.txt) and [`reports/spg_grid_extraction_summary.txt`](../../reports/spg_grid_extraction_summary.txt).
   - **Audit Discovery**: Uncovered a real metadata discrepancy in `rom_grid3d_25.09.spg` where the 1D Geoid Grid (`geoid_heights`) header declares `maxla = 29.83333333°` and `maxphi = 48.43333333°`, whereas the actual $154 \times 293$ array ends at `29.80000000°` and `48.40000000°` (off by 1 grid step of $2'$). Because `SpgReader` computes `maxla`/`maxphi` directly from `data.shape`, `PyTransDatRO` uses the true array boundaries and is immune to the header error.

---

## Verification & Testing
- **Test Modules**: [`tests/test_trans_ro.py`](../../tests/test_trans_ro.py), [`tests/test_api_legacy.py`](../../tests/test_api_legacy.py), [`tests/test_api_modern.py`](../../tests/test_api_modern.py)
- **Results**:
  - Python test suite: **`118/118 passed`** (`54` core tests, `17` legacy API tests, `28` modern API tests, `19` rate-limiter/config/feedback/telemetry tests).
  - ANCPI / Romgeo benchmark verification: **`PASS`** (`9,046` control points, zero numerical regression).
  - Frontend Vitest suite: **`67/67 passed`**.
  - Static analysis (`ruff check` & `ruff format --check`): **`All checks passed`**.
