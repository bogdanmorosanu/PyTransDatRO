# Task 50: Enforce 1-to-1 Line Parity & Centralized Coordinate/Angle Parsing

## Objectives
- [x] Resolve GitHub Issue #2: enforce strict 1-to-1 line parity in `/api/v1/transform/file` streaming transformations so empty lines and comments are yielded verbatim rather than dropped.
- [x] Resolve GitHub Issue #4: centralize coordinate tokenization and angle unit identification so place names with letters `d`, `m`, or `s` (e.g. *Bucuresti*, *Timisoara*, *Constanta*, *Iasi*) are correctly parsed as point IDs instead of DMS angles.
- [x] Create centralized coordinate and angle utilities in `webapp/src/utils/pointParser.ts` (`extractPointIdAndCoords`, `hasDmsSymbols`, `isPotentialAngle`) shared by Text and Files tabs.
- [x] Add automated test coverage in Python (`tests/test_api_modern.py`) and TypeScript (`webapp/src/utils/pointParser.test.ts`).
- [x] Audit GitHub issues #1, #3, and #5, confirming their resolution in earlier refactoring phases.
- [x] Synchronize living documentation in `docs/rest_api.md`, `docs/webapp/architecture.md`, and record task history.

## Implementation Details

1. **API Streaming File Parity (`api/routes/file.py`, `tests/test_api_modern.py`)**:
   - In `line_generator()`: Replaced silent omission of empty lines (`if not line: continue`) with `yield b"\n"`.
   - Streaming transformations of CSV and TXT files now maintain exact row alignment: input row $N$ corresponds strictly to output row $N$.
   - Updated `test_modern_file_streaming_csv` to verify empty line preservation.
   - Added `test_modern_file_streaming_preserves_blank_lines_and_alignment` asserting 1-to-1 line correspondence across a dataset with interleaved blank lines and `#` comments.

2. **Centralized Coordinate & Angle Tokenizer (`webapp/src/utils/pointParser.ts`)**:
   - Implemented `hasDmsSymbols(s)`: requires digits before letters `d`, `m`, and `s` (`/[°'"\u2018\u2019\u201C\u201D`]|(?:\d\s*[dms])/i`), preventing Romanian city names from triggering angle detection.
   - Implemented `isPotentialAngle(s)`: checks whether a token is a valid decimal number or DMS representation.
   - Implemented `extractPointIdAndCoords(tokens, isSourceGeo, order)`: reusable token classifier handling 2-column, 3-column, and 4-column data with or without leading Point IDs across both geodetic and Stereo 70 coordinate regimes.
   - Refactored `parseInputLines()` to use the shared extractor.

3. **Files Tab Input-Output Pairing (`webapp/src/components/FileUpload.ts`)**:
   - Replaced ad-hoc regex `/[°'"dms]/i` in `FileUpload.ts` with the centralized `extractPointIdAndCoords` and `detectDelimiter` from `pointParser.ts`.
   - Ensures uploaded files with city names containing `d`, `m`, `s` (such as `Bucuresti,325824.54,587904.09`) are plotted accurately on the map without being silently dropped.

4. **GitHub Issue Audit & Resolution Status**:
   - **Issue #1** (*No-data-zone points shown as successful conversions*): Already resolved in commit `860061a` (`isUntransformable = isOutOfGrid || isNoData`).
   - **Issue #3** (*Single-point result shows green "Valid coordinates" pill for out-of-grid*): Already resolved in commits `860061a` and `d131eb7` (`qualityStatus` status badge hierarchy).
   - **Issue #5** (*Transform clicks during running calculation dropped*): Already resolved in commit `d131eb7` (`this.pendingCalculation`).
   - **Issue #2** (*File endpoint drops blank lines*): Resolved and verified in this task.
   - **Issue #4** (*Uploaded points whose name contains d, m, s missing from map*): Resolved and verified in this task.

## Verification & Testing

- **Vitest Unit Suite (`webapp/src/utils/pointParser.test.ts`)**:
  - Added tests for `extractPointIdAndCoords` across Romanian city names, DMS symbol checks, and 2D/3D extractions.
  - All 88 tests in 6 test files passed in 371ms.
- **Frontend Compilation**:
  - `npm --prefix webapp run build` compiled cleanly with 0 TypeScript or Vite build errors.
- **Python Pytest Suite**:
  - `pytest tests/test_api_modern.py -k "streaming"` passed (6 tests).
  - Full suite passed: 206 tests passing in 21.5s with zero regressions.
- **Ruff Code Audit**:
  - `ruff check api tests` reported all checks passed.
