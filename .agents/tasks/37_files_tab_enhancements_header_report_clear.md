# Task 37: Files Tab Enhancements (Header Ignore, Multi-Line Quality Report & Remove All)

## Objectives
- [x] **End-to-End "Ignore header / first line" Toggle**:
  - Wired `#file-skip-header-toggle` via `apiClient.ts` to backend endpoint `POST /api/v1/transform/file`.
  - Added `skip_header: bool = Form(default=False)` to `api/routes/file.py`. When enabled, the first non-comment data line is emitted verbatim into the output stream without transforming it, without appending `# Invalid coordinate data` / warnings, and without logging server errors.
  - Excluded the header line from quality metrics (`valCount`, `outGridCount`, `errCount`) and coordinate dimensions (`count2D`, `count3D`) in `FileUpload.ts`.
  - Updated `parsePairedPoints` in `FileUpload.ts` to skip the header line when extracting map coordinates.
- [x] **Multi-Line Quality Report in Transformation Subsection**:
  - Upgraded Card 3 (Report) in `FileUpload.ts` to detect multi-file batches (`outputFiles.length > 1`).
  - Rendered a structured list displaying each processed file on its own line (`.quality-report-subrow`) with its truncated file name and semantic status pills (`Valid`, `Outside Grid`, `Invalid Input`, `Error`).
  - Added a summary line (`.quality-report-subrow.total-row`) separated by a dotted divider (`.quality-report-divider-dotted`) displaying grand total counts.
  - Preserved the compact single-line layout when exactly 1 file is transformed or ready to transform.
- [x] **"Remove All" Output List Action & Map Clearing**:
  - Renamed `clearOutput` button text from "Clear List" $\to$ "Remove All" (`Elimină Tot` in Romanian).
  - Aligned button styling to `class="btn btn-secondary btn-secondary-action"` matching `Download All` exactly.
  - Implemented `clearOutputFiles()` in `FileUpload.ts` to revoke blob object URLs (`URL.revokeObjectURL`) to prevent memory leaks, filter out all completed and error files from `fileQueue`, trigger `onClearFiles()` to purge markers from the Leaflet map layer, and reset the output view to the empty state.
- [x] **Automated Tests & Quality Assurance**:
  - Added `test_modern_file_skip_header_option` in `tests/test_api_modern.py` covering `skip_header=True` vs `skip_header=False`.
  - Added 2 new Playwright browser E2E scenarios in `webapp/e2e/file.spec.ts` covering header ignore with file download validation, multi-line report breakdown with total line, and "Remove All" clearing output cards and map markers.
  - Verified 100% pass on core geodetic suite (`tests/test_trans_ro.py`), API tests (`tests/test_api_modern.py`), Vitest suite (80 tests), and Ruff linter.
- [x] **Living Documentation**:
  - Updated `docs/webapp/actionable_inventory.md` with "Remove All", "Ignore Header" API integration, and multi-line report.
  - Updated `docs/rest_api.md` with `skip_header` parameter documentation.

---

## Verification & Testing
- **FastAPI Modern Suite (`pytest`)**:
  - `python -m pytest tests/test_api_modern.py` passed all 38 tests in 2.47s.
- **Core Geodetic Precision (`pytest`)**:
  - `python -m pytest tests/test_trans_ro.py` passed all 54 tests in 3.12s ($< 10^{-9}$ tolerance maintained).
- **Vitest Unit Suite**:
  - `npm --prefix webapp test` passed all 80 tests in 303ms.
- **Playwright End-to-End Suite**:
  - `npm --prefix webapp run test:e2e -- e2e/file.spec.ts` passed all 3 tests in 8.0s.
- **Ruff Linter**:
  - `python -m ruff check --no-cache api/ tests/ pytransdatro/` passed with 0 errors.
- **TypeScript & Vite Build**:
  - `npm --prefix webapp run build` passed with 0 errors in 224ms.
