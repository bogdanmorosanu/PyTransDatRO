# Dynamic SPG Grid Update Lifecycle & Unified Management Tooling

## Objectives
- [x] Analyze the entire grid lifecycle across core, API, web frontend, and benchmark testing suites when an official ANCPI `.spg` update is released.
- [x] Create a unified orchestration CLI tool ([`tools/update_grid.py`](../../tools/update_grid.py)) providing `status`, `inspect`, `apply`, `benchmark`, `verify`, and `rollback` actions.
- [x] Eliminate hardcoded grid names and stale fallbacks across API endpoints ([`api/routes/info.py`](../../api/routes/info.py), [`api/schemas/info.py`](../../api/schemas/info.py)) and frontend panels ([`webapp/src/components/FeedbackPanel.ts`](../../webapp/src/components/FeedbackPanel.ts), [`webapp/src/components/GeodeticHelpPanel.ts`](../../webapp/src/components/GeodeticHelpPanel.ts)).
- [x] Implement dynamic grid metadata synchronization, active grid vs. benchmark baseline comparison detection, and UI notice alerts.
- [x] Author comprehensive living documentation ([`docs/grid_lifecycle.md`](../../docs/grid_lifecycle.md)) detailing code mechanics, operator instructions, and the critical requirement for updating official ANCPI benchmark reference data.
- [x] Future-proof API and geodetic test assertions ([`tests/test_api_modern.py`](../../tests/test_api_modern.py)).
- [x] Create an agent skill ([`.agents/skills/grid-update/SKILL.md`](../skills/grid-update/SKILL.md)) that guides any future agent or operator through the grid update workflow, specifically enforcing awareness of the ANCPI reference test data refresh.
- [x] Execute full test suites: core geodetic precision, full pytest suite (207 tests), frontend Vitest suite (123 tests), and Playwright browser E2E suite (71 tests).

## Implementation Details

1. **Unified Orchestration Tool ([`tools/update_grid.py`](../../tools/update_grid.py))**:
   - Implemented `status` to audit active `.spg` file, metadata JSON consistency, comparison benchmark status, and test dataset timestamps.
   - Implemented `inspect` to validate candidate `.spg` arrays, dimensions, nodata masks, and display side-by-side release comparisons.
   - Implemented `apply` to safely move the active grid to [`pytransdatro/grids/archive/`](../../pytransdatro/grids/archive/) (preventing `AmbiguousGridError`), install the candidate, re-extract [`pytransdatro/grids/spg_metadata.json`](../../pytransdatro/grids/spg_metadata.json), audit client envelope bounds against `gridBounds.ts`, and emit operator reminders for expected test datasets.
   - Implemented `benchmark` to execute the full 9,046-point ANCPI validation against [`tests/data/*_expected.csv`](../../tests/data/), record active grid provenance (`grid_file`, `grid_version`, `grid_sha256`) in [`ancpi_comparison_summary.json`](../../pytransdatro/grids/ancpi_comparison_summary.json), and generate Leaflet HTML maps and text reports in `reports/`.
   - Implemented `verify` to run pytest and Vitest verification suites with a single command.
   - Implemented `rollback` to restore the latest archived grid from `pytransdatro/grids/archive/` and regenerate metadata.

2. **API Dynamic Resilience ([`api/routes/info.py`](../../api/routes/info.py), [`api/schemas/info.py`](../../api/schemas/info.py))**:
   - Refactored `_load_or_extract_metadata()` to check whether cached `spg_metadata.json` matches the active `.spg` file on disk. If missing or mismatched, dynamically extracts metadata on-the-fly and saves it.
   - Removed the hardcoded `"rom_grid3d_25.09.spg"` fallback dictionary.
   - Enhanced `GridComparisonResponse` schema with optional provenance fields (`grid_file`, `grid_version`, `grid_sha256`).

3. **Web Application Help & Feedback Panels ([`webapp/src/`](../../webapp/src/))**:
   - Updated [`webapp/src/components/GeodeticHelpPanel.ts`](../../webapp/src/components/GeodeticHelpPanel.ts) to display a baseline grid provenance chip and an amber notice banner if the active grid differs from the benchmark baseline dataset.
   - Added internationalized strings (`baselineGridLabel`, `benchmarkMismatchTitle`, `benchmarkMismatchDesc`) in [`webapp/src/i18n/en.ts`](../../webapp/src/i18n/en.ts) and [`webapp/src/i18n/ro.ts`](../../webapp/src/i18n/ro.ts).
   - Cleaned up [`webapp/src/components/FeedbackPanel.ts`](../../webapp/src/components/FeedbackPanel.ts) initial state to avoid hardcoding `rom_grid3d_25.09.spg`.
   - Updated [`webapp/src/utils/apiClient.ts`](../../webapp/src/utils/apiClient.ts) interface definitions.

4. **Living Documentation & Operator Guides ([`docs/grid_lifecycle.md`](../../docs/grid_lifecycle.md))**:
   - Authored a comprehensive technical reference detailing the `.spg` internal format, memory caching via `SpgReader`, FastAPI pre-warming, and dynamic API endpoints.
   - Provided an explicit step-by-step operator guide for ingesting new ANCPI grids.
   - Included a prominent notice detailing why `tests/data/*_expected.csv` must be regenerated via the official ANCPI desktop application on `tests/data/*_input.csv`.
   - Cross-linked in [`docs/grid.md`](../../docs/grid.md), [`docs/architecture.md`](../../docs/architecture.md), [`docs/testing.md`](../../docs/testing.md), and [`README.md`](../../README.md).

5. **Future-Proofing Test Suites ([`tests/test_api_modern.py`](../../tests/test_api_modern.py), [`tools/export_spg_metadata_json.py`](../../tools/export_spg_metadata_json.py))**:
   - Replaced version-specific assertions (e.g. `meta["release"]["major"] == 25`, node counts `2689`, `43071`) with dynamic threshold and structural assertions.
   - Updated `tools/export_spg_metadata_json.py` to auto-discover any single `.spg` file rather than falling back to a hardcoded filename.

6. **Official Grid Update Agent Skill ([`.agents/skills/grid-update/SKILL.md`](../skills/grid-update/SKILL.md))**:
   - Created the dedicated agent skill `grid-update` to orchestrate grid replacements for any AI agent or operator.
   - Specifically enforces awareness that updating the binary `.spg` grid requires refreshing official ANCPI ground truth coordinates in [`tests/data/*_expected.csv`](../../tests/data/).
   - Instructs the agent to prompt the user about test dataset readiness before running benchmarks, outlines pre-flight inspection via `tools/update_grid.py inspect`, safe archiving, dynamic metadata extraction, report updates, quality verification, and disaster rollback.
   - Formally registered in the project workflow in [`GEMINI.md`](../../GEMINI.md).

## Verification & Testing
- **Unified Tooling**:
  - `python tools/update_grid.py status`: PASS (Active grid, metadata, comparison summary, and test datasets all verified).
  - `python tools/update_grid.py inspect pytransdatro/grids/rom_grid3d_25.09.spg`: PASS (Validated 2D and 1D layers, extents, and Helmert parameters).
  - `python tools/update_grid.py apply pytransdatro/grids/rom_grid3d_25.09.spg`: PASS (Detected identical active grid, refreshed metadata).
  - `python tools/update_grid.py benchmark`: PASS (9,046 points forward and inverse with zero flag 4 residuals).
  - `python tools/update_grid.py rollback`: PASS (Proper error handling when archive is empty).
- **Python Regression & API Suite**:
  - `pytest`: 207 passed in 19.35s (100% pass rate).
- **Frontend Unit Suite**:
  - `npm --prefix webapp test`: 123 passed across 10 test files.
- **Frontend Production Build**:
  - `npm --prefix webapp run build`: 0 TypeScript errors, clean Vite compilation.
- **Browser E2E Suite**:
  - `npm --prefix webapp run test:e2e`: 71 passed across all 8 workers in 57.1s.
- **Code Audit & Ruff**:
  - `ruff check` & `ruff format --check`: All checks passed.
  - Machine-specific path audit: 0 hardcoded absolute paths.
