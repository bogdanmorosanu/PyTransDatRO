# REST API Test Matrix Completion, Live Legacy Differential Parity & Non-Transformation Service Verification

## Objectives
- [x] Fix subtle validation edge cases in the modern REST API (`PointRequest.coos: List[Any]` inline warning consistency in [`api/schemas/point.py`](../../api/schemas/point.py) and `math.isfinite` guard on the 3-token Point ID heuristic in [`api/routes/file.py`](../../api/routes/file.py)).
- [x] Complete the local Legacy `TransDatOnline` compatibility test matrix in [`tests/test_api_legacy.py`](../../tests/test_api_legacy.py) across `GET` and `POST`, `Stereo70ToETRS89` and `ETRS89ToStereo70`, `2D` and `3D`, `"Out of grid"` (`OutOfGrid2DErr`, `OutOfGrid1DErr`, curved-border pre-check wedge), `"No data on grid"`, `"Invalid coordinate data"`, `"Wrong dimension"`, string-encoded numbers in `POST`, obsolete `Stereo30*` operations, and `multipart/form-data`.
- [x] Implement live **Differential Parity Tests** (`@pytest.mark.legacy_live`) in [`tests/test_api_legacy.py`](../../tests/test_api_legacy.py) comparing live responses from `https://www.geo-spatial.org/transdatonline/cooOpService` against the local FastAPI service within a $5\text{ mm}$ planar tolerance (`0.005 m` / $\approx 7.85 \times 10^{-10}\text{ rad}$) and $1.0\text{ m}$ vertical tolerance (EGG97 vs. modern SPG quasigeoid).
- [x] Implement an automatic reachability probe and decommissioning skip mechanism (`_resolve_live_legacy_url`, `PYTRANSDAT_LEGACY_LIVE_ENABLED`, `PYTRANSDAT_LEGACY_LIVE_URL`) and register the `legacy_live` marker in [`pyproject.toml`](../../pyproject.toml) so the test suite skips cleanly once the legacy service is decommissioned.
- [x] Complete the Modern API (`/api/v1/transform/point`, `/batch`, `/file`) test matrix in [`tests/test_api_modern.py`](../../tests/test_api_modern.py) covering all engine-level exceptions (`OutOfGrid2DErr`, `OutOfGrid1DErr`, curved-border wedge points, northern parallel bulge points, `NoDataGridErr`, `InvalidCoordinateErr`, `WrongDimensionErr`), file upload size (`HTTP 413`) and row truncation limits.
- [x] Expand non-transformation service test coverage across `/health` (degraded states), `/api/v1/telemetry/*` (seeded SQLite `usage.db` with 4-quadrant spatial data), `/api/v1/feedback` (mocked SMTP TLS/auth delivery and failure fallback in [`tests/test_feedback.py`](../../tests/test_feedback.py)), and `SlidingWindowRateLimiter` (per-IP isolation, shared batch bucket, unmetered routes, and window expiration in [`tests/test_api_rate_limiter.py`](../../tests/test_api_rate_limiter.py)).
- [x] Update living documentation in [`docs/rest_api.md`](../../docs/rest_api.md) and [`docs/testing.md`](../../docs/testing.md).

---

## Implementation Details

### 1. API Validation Refinements
1. **Modern Single-Point Inline Warning Consistency ([`api/schemas/point.py`](../../api/schemas/point.py))**:
   - Changed `PointRequest.coos` from `List[float]` to `List[Any]`.
   - Previously, Pydantic coerced booleans (`True -> 1.0`, `False -> 0.0`) and rejected `null` or non-numeric strings at the outer HTTP layer with `422 VALIDATION_ERROR` on `/api/v1/transform/point`, while `/api/v1/transform/batch` returned HTTP `200` with `{"coos": [], "warning": "Invalid coordinate data"}`.
   - With `coos: List[Any]`, all coordinate elements pass through `execute_point_transformation()` $\to$ `parse_numeric_coord_list()`, rejecting `bool`, `None`, non-numeric strings, `NaN`, and `Infinity` uniformly with HTTP `200` and `{"coos": [], "warning": "Invalid coordinate data"}`.
2. **Streaming File 3-Token Non-Finite Guard ([`api/routes/file.py`](../../api/routes/file.py))**:
   - Added `math.isfinite()` checks inside the 3-token Point ID disambiguation heuristic in `_transform_line()`.
   - Prevents lines such as `-Inf,500000,500000` (where `-Inf < 50000`) or `Inf,46.0,25.0` (where `abs(Inf) > 90`) from misclassifying `-Inf`/`Inf` as a numeric Point ID and transforming the remaining two tokens as a valid 2D point.

### 2. Legacy API Matrix & Live Differential Parity ([`tests/test_api_legacy.py`](../../tests/test_api_legacy.py))
- **Complete Local Emulation Matrix**:
  - Verified both `/transdatonline/cooOpService` and root alias `/cooOpService` across `GET` and `POST`, `Stereo70ToETRS89` and `ETRS89ToStereo70`, `2D` and `3D`.
  - Tested all warning states: `"Out of grid"` (`OutOfGrid2DErr`, `OutOfGrid1DErr` at `[778000.0, 500000.0, 100.0]` and `[48.45°, 25.00°, 100.0]`, and curved-border wedge points `[205000.0, 260000.0, 100.0]` and `[43.50°, 19.90°]`), `"No data on grid"` (`[549364.0, 147553.0]` and `[0.8087707898617313, 0.35742127559921274]`), `"Invalid coordinate data"`, `"Wrong dimension"` (`0D`, `1D`, `4D`), string-encoded numbers in `POST` (`{"coos": ["500000.0", "500000.0", "100.0"]}`), and obsolete `Stereo30ToETRS89` / `ETRS89ToStereo30`.
- **Live Differential Parity Suite (`@pytest.mark.legacy_live`)**:
  - Compares live HTTP responses from `https://www.geo-spatial.org/transdatonline/cooOpService` against local `TestClient(app)` responses across 16 `GET` scenarios and 2 multi-point `POST` batches.
  - Validates 2D coordinates within $5\text{ mm}$ (`PLANAR_TOL_M = 0.005 m`, `RAD_TOL = 0.005 / 6371000.0 rad`) and 3D elevations within $1.0\text{ m}$ (`ELEV_TOL_M = 1.0 m`).
  - Accounts for the legacy Java in-place `Coordinate` mutation bug on `ETRS89ToStereo70` warning responses (where the legacy Java code mutates input radians into degrees or Stereo70 meters before catching the grid error, whereas `PyTransDatRO` cleanly echoes the original input coordinates).
  - Guards live tests with `_resolve_live_legacy_url()`, automatically calling `pytest.skip(...)` when `PYTRANSDAT_LEGACY_LIVE_ENABLED=0` or when the remote endpoint becomes unreachable/decommissioned.

### 3. Modern API & Non-Transformation Service Tests
- **Modern Transformation Suite ([`tests/test_api_modern.py`](../../tests/test_api_modern.py))**:
  - Added end-to-end engine-to-API warning translation tests for `/api/v1/transform/point`, `/api/v1/transform/batch`, and `/api/v1/transform/file` across `NoDataGridErr`, `OutOfGrid1DErr`, SW curved-border wedge (`43.50°N, 19.90°E`), northern parallel bulge (`48.50°N, 25.00°E`), 3-token `-Inf`/`Inf`/`NaN` CSV lines, HTTP `413` upload size limits, and `max_file_points` stream truncation.
  - Added degraded `/health` probe verification (`app.state.trans = None`) and seeded SQLite `usage.db` verification for `/api/v1/telemetry/stats`, `/api/v1/telemetry/timeline` (daily/weekly/monthly rollups and filters), and `/api/v1/telemetry/spatial-matrix` (4-quadrant percentage math and `metric=points` vs. `metric=calls`).
- **Feedback & Mailer Suite ([`tests/test_feedback.py`](../../tests/test_feedback.py))**:
  - Added mocked `smtplib.SMTP` tests verifying TLS/auth email dispatch (`status="success"`, `delivered_email=True`), `X-Forwarded-For` client IP extraction, resilient fallback to `status="logged"` when `smtplib.SMTP` raises `SMTPException`, and Pydantic schema bounds (`report_type` and `max_length=5000`).
- **Sliding-Window Rate Limiter Suite ([`tests/test_api_rate_limiter.py`](../../tests/test_api_rate_limiter.py))**:
  - Added tests for `/api/v1/transform/point` (60 req/min), shared batch bucket across `/api/v1/transform/batch`, `/api/v1/transform/file`, and `POST /transdatonline/cooOpService` (10 req/min), per-IP isolation via `X-Forwarded-For`, unmetered route pass-through (`/health`, `/api/v1/grid/info`), and sliding-window expiration recovery.

---

## Verification & Testing
- **Python Test Suite (`pytest`)**: `157/157 passed` (`54` core in `test_trans_ro.py`, `42` legacy local + live differential in `test_api_legacy.py`, `36` modern API in `test_api_modern.py`, `7` feedback in `test_feedback.py`, `7` config in `test_config.py`, `6` telemetry in `test_telemetry.py`, `5` rate limiter in `test_api_rate_limiter.py`).
- **ANCPI Official Benchmark**: `PASS` across all 9,046 national control points (`tools/generate_ancpi_comparison_reports.py`).
- **Frontend Unit Tests (`vitest`) & Build**: `67/67 passed` (`npm --prefix webapp test`) and zero TypeScript errors (`npm --prefix webapp run build`).
- **Static Analysis (`ruff`)**: Zero lint or formatting issues across all 51 Python files.
