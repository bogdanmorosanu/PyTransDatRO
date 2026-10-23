# Testing in PyTransDatRO

This document summarizes the testing structure, datasets, tolerances, and test execution patterns implemented in `PyTransDatRO` (primarily in `tests/test_trans_ro.py` and the `tests/` directory).

---

## 1. Test Data & Fixtures

The test suite relies on official reference coordinates and verification thresholds based on the official ROMPOS documentation (*Help_TransDatRO_code_source_EN.pdf*) and outputs from the official `TransDatRO v4.08` software.

### A. Verification Tolerances

Tolerances define the maximum acceptable discrepancy against official reference outputs, matching the criteria specified in *Help_TransDatRO_code_source_EN.pdf* (page 3) and sub-millimeter round-trip precision:

| Tolerance Fixture | Value | Physical Meaning | Description |
|-------------------|-------|------------------|-------------|
| `coo_rad_tol` | `0.00000000014544410` rad | $\approx 0.00003''$ | Angular tolerance for geographic latitude and longitude against official software. |
| `coo_plan_tol` | `0.003` m | $3\text{ mm}$ | Planar tolerance for Stereo70 Northing and Easting coordinates against official software. |
| `coo_elev_tol` | `0.003` m | $3\text{ mm}$ | Vertical tolerance for elevation / height components ($Z$ and $h$). |
| `roundtrip_plan_tol` | `0.0005` m | $0.5\text{ mm}$ | Planar tolerance for forward-and-back round-trip conversions (accommodating the algebraic linear inverse in Helmert 2D). |
| `roundtrip_rad_tol` | `1e-10` rad | $\approx 2.06 \times 10^{-5}''$ | Angular tolerance for forward-and-back round-trip conversions (matching the 6th decimal of an arcsecond). |

### B. External Benchmark Files

For extensive surface verification across Romania, benchmark datasets are organized in a dedicated directory under `tests/data/` using a standardized naming convention (`<source>_to_<target>_<role>.csv`):

- **Stereo70 to ETRS89 Transformation**:
  1. **`tests/data/st70_to_etrs89_input.csv`**:
     - Contains 9,046 Stereo70 coordinates formatted as CSV (`id,n,e,z`).
     - Generated across a half-grid (nodes and mid-points) spanning the entire 2D shift grid extent.
  2. **`tests/data/st70_to_etrs89_expected.csv`**:
     - Contains the corresponding 9,046 ETRS89 coordinates processed through official `TransDatRO v4.08`.
     - Encoded in ANSI/Windows-1252, with coordinates in sexagesimal DMS strings (e.g., `47°42'56.40000"N, 22°28'31.99998"E`).
- **ETRS89 to Stereo70 Transformation**:
  3. **`tests/data/etrs89_to_st70_input.csv`**: Contains ETRS89 coordinates in sexagesimal DMS format for reverse transformation inputs.
  4. **`tests/data/etrs89_to_st70_expected.csv`**: Reference Stereo70 outputs for reverse transformation verification.

---

## 2. Test Execution Patterns

All tests are implemented using `pytest` and execute against the top-level coordinator class `pytransdatro.TransRO`.

### A. Round-Trip Forward-and-Back Testing

To ensure mathematical reversibility without relying on hardcoded discrete reference points, two comprehensive round-trip tests process the bulk coordinate datasets:

- **`test_st70_to_etrs89_to_st70_roundtrip`**:
  - Transforms $(N, E, Z) \rightarrow (\text{lat}, \text{lon}, h) \rightarrow (N', E', Z')$.
  - Verifies planar coordinates match within `roundtrip_plan_tol` ($0.5\text{ mm}$) and vertical elevation matches within $10^{-9}\text{ m}$.
- **`test_etrs89_to_st70_to_etrs89_roundtrip`**:
  - Transforms $(\text{lat}, \text{lon}, h) \rightarrow (N, E, Z) \rightarrow (\text{lat}', \text{lon}', h')$.
  - Verifies geographic coordinates match within `roundtrip_rad_tol` ($10^{-10}\text{ rad}$) and ellipsoidal height matches within $10^{-9}\text{ m}$.

### B. Bulk Grid Validation (`test_st70_to_etrs89_2D_fromfile`)

- Ingests the 7,866 points from `tests/data/st70_to_etrs89_input.csv` and `tests/data/st70_to_etrs89_expected.csv`.
- Converts DMS strings to radians using `pytransdatro.utils.sexa_dms_chars_to_rad()`.
- Executes `st70_to_etrs89(n, e)` on all points.
- Validates that every point across the grid matches the official `TransDatRO` output within `coo_rad_tol`.

### C. Boundary & Exception Handling (Parametrized)

Edge cases and out-of-boundary behavior are tested using `@pytest.mark.parametrize`:

- **Out of Grid Boundaries (`test_st70_to_etrs89_extent_inner_out_of_grid`)**:
  - Evaluates a 16-point Cartesian product ($4 \times 4$ of $N$ and $E$) placed at or outside the outer grid bounding box.
  - Asserts that `pytransdatro.exceptions.OutOfGridErr` is raised.
- **No-Data Zones Within Grid Extents (`test_st70_to_etrs89_no_grid_inner_no_data`)**:
  - Evaluates 16 points located inside the grid bounding box, but positioned where the binary grid cells contain no valid interpolation data.
  - Asserts that `pytransdatro.exceptions.NoDataGridErr` is raised.
- **National Border Transition Areas (`test_st70_to_etrs89_just_out_ro_border_no_data`)**:
  - Tests 4 coordinates located immediately outside the Romanian territorial boundary.
  - Asserts that `pytransdatro.exceptions.NoDataGridErr` is raised.

### D. Test Data Generation Helper

- **`create_transdatro_test_file(file, template='st70_half_grid')`**:
  - A utility function contained within `test_trans_ro.py` (not a test itself).
  - Iterates through the active 2D grid (`_t_gr2d`), generating points at grid nodes and half-step intervals.
  - Outputs a CSV file that can be loaded directly into the official TransDatRO desktop software to produce new benchmark files.

---

## 3. Web API, Configuration & Telemetry Test Suite

The API layer is thoroughly tested using `pytest` and `starlette.testclient.TestClient`:

| Module | Tests | Focus Areas |
| :--- | :--- | :--- |
| **`tests/test_trans_ro.py`** | 54 | Core 2D/3D transformations, round-trip accuracy, 6-point extended dual-CRS bounding boxes, `OutOfGrid2DErr` vs `OutOfGrid1DErr`, pole singularity/periodic aliasing guards, 4-stage `transform_batch` fault isolation, sequence length/finite number validation (`InvalidCoordinateErr`, `WrongDimensionErr`), and strict DMS parsing. |
| **`tests/test_api_modern.py`** | 36 | Point (`/api/v1/transform/point`), Option B fault-isolated batch (`/transform/batch`), file streaming (`/transform/file`), `No data on grid` / `OutOfGrid1DErr` / curved-border wedge / northern bulge warnings, 3-token `-Inf`/`Inf`/`NaN` file guards, HTTP `413` size & row truncation limits, grid metadata (`/grid/info`, `/grid/metadata`, `/grid/comparison`), seeded SQLite telemetry (`/telemetry/stats`, `/telemetry/timeline`, `/telemetry/spatial-matrix`), health (`/health` healthy/degraded), and structured validation errors (`ErrorResponse`). |
| **`tests/test_api_legacy.py`** | 42 | Full local and live differential parity suite for `/transdatonline/cooOpService` and `/cooOpService` across `GET` & `POST`, `Stereo70ToETRS89` & `ETRS89ToStereo70`, `2D` & `3D`, `Out of grid` (2D/1D/wedge), `No data on grid`, `Invalid coordinate data`, `Wrong dimension`, string-encoded numbers in `POST`, obsolete `Stereo30*`, `multipart/form-data`, and live comparison against `https://www.geo-spatial.org/transdatonline/cooOpService` ($5\text{ mm}$ planar / $1.0\text{ m}$ vertical tolerance) with automatic decommissioning skip (`@pytest.mark.legacy_live`). |
| **`tests/test_api_rate_limiter.py`** | 5 | Sliding-window in-memory rate limiting enforcement (60/min point, 10/min batch shared across `/batch`, `/file`, and legacy `POST`), per-IP `X-Forwarded-For` isolation, unmetered route pass-through, `Retry-After` headers, and window expiration recovery. |
| **`tests/test_config.py`** | 8 | Centralized `Settings` resolution, environment variable overrides (`.env`), `webapp_dist_dir` default and relative override anchoring, CORS parsing, and `dictConfig` logging dictionary generation. |
| **`tests/test_feedback.py`** | 7 | User feedback submission (`/api/v1/feedback`), honeypot bot misdirection, schema validation, file persistence (`feedback.jsonl`), mocked SMTP TLS/auth delivery, and SMTP failure fallback. |
| **`tests/test_telemetry.py`** | 6 | Asynchronous background SQLite WAL logging, transaction queues, rolling dynamic window retention, and quadrant coordinate aggregation. |

**Execution Commands**:
```powershell
# Run full test suite (including live legacy differential parity if reachable)
D:\anaconda3\envs\pytransdat\python.exe -m pytest -q

# Run offline only (explicitly skipping live legacy server tests)
D:\anaconda3\envs\pytransdat\python.exe -m pytest -m "not legacy_live" -q
```
*(Runs all 158 Python unit, integration, contract, and live differential parity tests. Live legacy tests automatically skip via `pytest.skip` if `PYTRANSDAT_LEGACY_LIVE_ENABLED=0` or once the legacy server is decommissioned).*

---

## 4. Frontend Utility & UI Test Suite (`webapp/`)

Client-side parsing, serialization, and UI factory helpers are tested with **Vitest**:

| Test File | Tests | Coverage Scope |
| :--- | :--- | :--- |
| **`src/utils/dmsParser.test.ts`** | 13 | Unspaced DMS parsing (`DD°MM'SS.SSSSS"`), edge cases, whitespace tolerance, sexagesimal to decimal conversions. |
| **`src/utils/pointParser.test.ts`** | 28 | Auto-delimiter detection (tab, space, comma, semicolon), Point ID extraction, 2D vs 3D tuples, coordinate ordering (`NE(H)` vs `EN(H)`). |
| **`src/utils/gridBounds.test.ts`** | 10 | Romania national bounding envelope checks, coordinate bounds validation, out-of-grid pre-flight filters. |
| **`src/utils/exporters.test.ts`** | 9 | Pure client-side export generation for CSV, GeoJSON, and AutoCAD DXF R12 ASCII formats (including out-of-grid and no-data warnings). |
| **`src/ui/ui.test.ts`** | 6 | Accessible UI factory helpers (`renderDropdownMenu`, `renderSegmentedControl`, `renderSelect`, `renderCopyBadge`, `renderCopyButton`). |

**Execution Commands**:
```powershell
# Run Vitest suite once
npm --prefix webapp test

# Run Vitest in interactive watch mode
npm --prefix webapp run test:watch
```
*(Runs all 80 TypeScript unit tests in $< 350\text{ ms}$).*

---

## 5. Official ANCPI / Romgeo Benchmark Verification

Coordinate transformations must maintain 100% byte-for-byte fidelity against ANCPI benchmark reference sets:

```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/generate_ancpi_comparison_reports.py
```

- **Reference Dataset**: 9,046 national geodetic control points across Romania.
- **Criteria for Success**:
  - Horizontal planar discrepancy: $< 10^{-6}\text{ m}$ ($< 0.001\text{ mm}$).
  - Vertical ellipsoidal height discrepancy: $< 10^{-6}\text{ m}$.
  - Overall status: `PASS` exported to `pytransdatro/grids/ancpi_comparison_summary.json`.

---

## 6. End-to-end UI tests (Playwright)

Browser-level flow tests live in `webapp/e2e/*.spec.ts` (Chromium only, no retries). Vitest remains
the default for parsers and `src/ui/` factories.

**One-time setup**:
```bash
npm --prefix webapp install
npx --prefix webapp playwright install chromium
```

**Run Options**:

- **Recommended Single-Command Runners (Auto-detects Python & Conda)**:
  ```powershell
  # Windows PowerShell
  .\test-e2e.ps1         # Run all 59 E2E browser tests across 11 suites (headless)
  .\test-e2e.ps1 -Ui     # Open interactive visual runner with time-travel debugging
  .\test-all.ps1         # Run Python suite + Vitest unit tests + Playwright E2E in one go
  ```
  ```bash
  # macOS / Linux (POSIX)
  ./test-e2e.sh          # Run all 59 E2E browser tests across 11 suites (headless)
  ./test-e2e.sh --ui     # Open interactive visual runner
  ./test-all.sh          # Run full verification suite
  ```

- **Direct NPM Command** (Python is auto-detected by `webapp/playwright.config.ts`, or overridden via `PYTHON`):
  ```bash
  npm --prefix webapp run test:e2e
  ```
The script typechecks `webapp/e2e` (`tsc -p e2e`, separate from the app `tsconfig.json`) and then runs Playwright.

**Server**: `playwright.config.ts` starts one process from the repo root: it builds the frontend
with `VITE_E2E=true` into `webapp/dist-e2e/` (`npm run build:e2e`) and then serves web + API with
`python -m uvicorn api.main:app` on `http://127.0.0.1:8100`. `reuseExistingServer` is false, so the
normal dev server on :8000 is never touched.

**Isolation** (env set by the config): `PYTRANSDAT_WEBAPP_DIST_DIR=webapp/dist-e2e`,
`PYTRANSDAT_LOG_DIR` and `PYTRANSDAT_FEEDBACK_LOG_PATH` under `webapp/e2e/.tmp/` (git-ignored),
`PYTRANSDAT_SMTP_HOST=""` (no mail), and rate limits raised to 100000.

**Fixtures** (`webapp/e2e/fixtures.ts`; import `test`/`expect` from there, not from Playwright):
- External hosts (map tiles, Google Fonts, ...) are intercepted and answered with a stub (1x1
  transparent PNG / empty CSS) carrying an `x-e2e-stub` header; `externalRequests` lists them.
- Language is pinned to `en` via `localStorage['pytransdat_lang']` (only if unset, so a test that
  switches language keeps it across reloads).
- `consoleErrors` collects console errors and uncaught page errors.
- `mapView()` / `readMapView(page)` return `window.__e2e.mapView()`.

**`window.__e2e`**: exposed only when the bundle is built with `VITE_E2E=true` (flag
`E2E_ENABLED` in `webapp/src/config.ts`); absent from the normal production build. `mapView()`
returns a plain snapshot: `center`, `zoom`, `bounds`, `maxBounds` (or `null`), `markerCount` and
`markersByLayer`. Use it to assert that transforms never auto-pan or zoom the map.

**Reference data**: `webapp/e2e/referencePoints.ts` holds ANCPI benchmark points (decimal degrees
and metres) and `matchesAtDisplayPrecision()` for comparisons at UI precision (8 decimals degrees,
3 decimals metres, +-1 in the last digit).

**Failures**: traces (`retain-on-failure`) and screenshots (`only-on-failure`) are written to
`webapp/test-results/`, the HTML report to `webapp/playwright-report/` (both git-ignored). Open a
trace with `npx --prefix webapp playwright show-trace <trace.zip>`.

---

## 7. Code Hygiene & Static Analysis Audits

Before committing or releasing changes, execute the static quality audit:

```powershell
# Ruff linter and formatter check
D:\anaconda3\envs\pytransdat\python.exe -m ruff check pytransdatro api tests tools
D:\anaconda3\envs\pytransdat\python.exe -m ruff format --check pytransdatro api tests tools

# TypeScript compilation
npm --prefix webapp run build
```

