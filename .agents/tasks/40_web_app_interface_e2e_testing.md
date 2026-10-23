# Task 40: Web App Interface End-to-End Testing & Negative Interaction Test Suite

## Objectives
- [x] **Living Documentation Alignment**: Enforce the Empty Fields Rule across all living documentation (`actionable_inventory.md`, `architecture.md`, `design_system.md`, `archetypes.md`), establishing that invalid, cleared, or untransformable coordinate output fields must remain empty (`""`), not em-dash (`—`).
- [x] **Test Fixtures & Programmatic File Helpers**: Create static test fixtures in `webapp/e2e/fixtures/files/` and programmatic ingestion helpers in `webapp/e2e/helpers/files.ts` covering valid 2D/3D datasets, syntax corruptions, impossible angles (>90°), out-of-grid bounds, mixed errors, and headers.
- [x] **Point Mode Negative & Edge Verification (`point.spec.ts`)**: Expand Point test suite from 7 to 18 tests covering syntax resilience, out-of-grid boundary rejection, impossible angles, 2D vs 3D elevation toggling, quasigeoid undulation $\zeta$ with explicit sign, smart multi-token pasting, Enter key triggering, clipboard copying, and order/unit reformatting.
- [x] **Text / Batch Mode Lifecycle & Error Integrity (`batch.spec.ts`)**: Expand Text mode suite to 9 tests verifying auto-delimiter detection (tabs, spaces, commas, semicolons), clean error formatting (`Point ID, # [Reason]`), 1-to-1 line preservation, vertical height persistence across transformations, settings-shift immutability, and clean CAD/GIS export filtering (CSV, GeoJSON, DXF).
- [x] **Files Mode Staged Batch & Queue Lifecycle (`file.spec.ts`)**: Implement 7 E2E tests for staged multi-batch transformations with mid-session settings changes, 1-to-1 error row preservation, header skipping, queue item removal, queue clearing, and large dataset volume warning modals (>2,000 points).
- [x] **Map Interactions & Spatial Tooling (`map.spec.ts`)**: Expand Map suite to 12 tests validating viewport stability, out-of-grid click protection, untransformable coordinate erasing, interactive ruler 3-click measurement lifecycle with Escape cancellation, layer drawer styling and swatch selection persisted in `localStorage`, and marker geodetic attribute popups with clipboard actions.
- [x] **App Lifecycle, Health & Diagnostic Shell (`lifecycle.spec.ts`, `errors.spec.ts`)**: Author tests covering Brand Logo fresh reset, soft title navigation, offline alert banner reactivity, floating panels search filtering and live benchmark cards, bot challenge verification, and backend HTTP 500 / HTTP 429 rate limit error handling.
- [x] **Quality Gates & Suite Verification**: Achieve 100% pass rate across the full Playwright browser suite (59/59 tests), Vitest unit suite (80/80 tests), clean TypeScript build, and Ruff linting.

---

## Implementation Details

### 1. Living Documentation Harmonization (`docs/webapp/`)
- **Empty Fields Output Invariant**:
  - [`docs/webapp/actionable_inventory.md`](../../docs/webapp/actionable_inventory.md): Updated lines 77, 78, 124, 272 to specify empty fields (`""`) instead of em-dash (`—`); updated test references to `lifecycle.spec.ts`.
  - [`docs/webapp/architecture.md`](../../docs/webapp/architecture.md): Standardized output coordinate failure behavior to empty strings.
  - [`docs/webapp/design_system.md`](../../docs/webapp/design_system.md) & [`docs/webapp/archetypes.md`](../../docs/webapp/archetypes.md): Documented that output input fields clear to empty on transform failures.
- **Testing & Webapp Documentation**:
  - [`webapp/README.md`](../../webapp/README.md): Documented 59 E2E tests organized into 11 suites with updated execution commands.
  - [`docs/testing.md`](../../docs/testing.md): Updated E2E runner specifications to 59 tests and Vitest unit suite to 80 tests.

### 2. E2E Test Fixtures & Helpers (`webapp/e2e/fixtures/`, `webapp/e2e/helpers/`)
- **Static Fixture Files (`webapp/e2e/fixtures/files/`)**:
  - `valid_etrs_3d.csv`: Clean ETRS89 3D points for forward transformations.
  - `valid_stereo70_3d.csv`: Clean Stereo 70 3D points for reverse transformations.
  - `corrupted_syntax.txt`: Malformed rows with non-numeric tokens and insufficient coordinates.
  - `impossible_angles.csv`: Geographic coordinates exceeding physical bounds ($\ge 90^\circ$ latitude).
  - `out_of_grid.csv`: Coordinates outside the Romanian national bounding box.
  - `mixed_errors.csv`: Interleaved valid points, syntax errors, and out-of-grid coordinates.
  - `headed_points.csv`: Files with column headers to verify the skip-header toggle.
- **Programmatic Helpers (`webapp/e2e/helpers/files.ts`, `webapp/e2e/helpers/ui.ts`)**:
  - `uploadMemoryFiles`: Programmatically constructs and uploads in-memory file payloads.
  - `uploadDiskFixtures`: Uploads committed static fixture files from disk.
  - `readDownloadContent`: Awaits download events and reads text content for CAD/GIS export verification.
  - `setDirection`, `setOrder`, `setUnit`: Modular UI state setters for direction, coordinate order, and angular units.
- **Map Click Integration (`webapp/e2e/fixtures.ts`)**:
  - Added `fireMapClick(page, lat, lng)` backed by the `window.__e2e.fireMapClick` bridge for deterministic Leaflet map interaction without canvas drift.

### 3. Comprehensive E2E Test Suites (`webapp/e2e/`)
- **Point Suite (`point.spec.ts` - 18 tests)**:
  - Validates forward and reverse conversions against ANCPI reference values.
  - Tests negative boundary conditions: syntax errors (`abc`, `45.0,notanumber`), impossible angles ($95^\circ, 25^\circ$), out-of-grid bounds ($52^\circ\text{ N}, 10^\circ\text{ E}$), and ensures output fields remain empty.
  - Verifies 2D vs 3D elevation input disabling/enabling and quasigeoid undulation $\zeta$ sign preservation.
  - Tests Enter key submission, copy coordinates button, and order/unit reformatting.
- **Text / Batch Suite (`batch.spec.ts` - 9 tests)**:
  - Verifies multi-delimiter parsing across tabs, spaces, commas, and semicolons.
  - Tests clean error formatting: invalid rows omit coordinate numbers and output only `Point ID, # [Reason]`.
  - Asserts 1-to-1 line correspondence preservation between input and output.
  - Verifies settings-shift invariant: changing direction or unit does not mutate existing output until `Transform` is clicked.
  - Validates export formats: CSV, GeoJSON, and DXF R12 ASCII export clean geometric entities while omitting syntax errors.
- **Files Suite (`file.spec.ts` - 7 tests)**:
  - Implements staged multi-batch transformations: File 1 under Settings A $\to$ transform $\to$ switch to Settings B $\to$ File 2 uploaded and transformed $\to$ File 1 preserves Settings A output while File 2 uses Settings B.
  - Tests mixed error handling with 1-to-1 row preservation.
  - Tests queue lifecycle: individual file deletion and `Clear queue` button.
  - Tests volume warning modal for files exceeding 2,000 points.
- **Map Suite (`map.spec.ts` - 12 tests)**:
  - Validates viewport stability: transformations do not auto-pan or zoom the map.
  - Tests out-of-grid click protection: clicking outside Romania displays `Outside Grid` warning and clears output fields.
  - Tests interactive ruler measurement lifecycle: 3-click distance accumulation, Escape key cancellation, and reset button.
  - Tests layer drawer: swatch selection, visibility toggling, and persistence in `localStorage`.
  - Tests marker interaction: geodetic attribute popup displays accurate coordinates with clipboard copy action.
- **App Lifecycle & Errors (`lifecycle.spec.ts`, `errors.spec.ts` - 7 tests)**:
  - Tests Brand Logo click: triggers fresh first-load state reset, restoring default coordinates, direction, order, and unit.
  - Tests soft navigation between workbench tabs and subpages.
  - Tests offline alert banner and health dot reactivity.
  - Tests floating panels search filtering and live benchmark cards.
  - Tests backend error resilience: HTTP 500 error display and HTTP 429 rate limiting toast notification.

### 4. Application Bug Fixes & Refinements (`webapp/src/`)
- **`MapView.ts` Ruler Tooltip Initialization**:
  - Ensured measurement tooltip is created immediately on the second click during ruler operation even if no `mousemove` event was dispatched, preventing null reference errors during programmatic testing.
- **`main.ts` Fresh Reset State**:
  - Strengthened `resetToFreshState()` to reset projection direction (`ETRS89ToStereo70`), coordinate order (`NE`), and angular unit (`degrees`), resetting Point inputs to default projection center.
- **`.gitignore`**:
  - Configured gitignore to preserve `webapp/e2e/fixtures/files/` while properly ignoring Playwright test results and reports.

---

## Verification & Quality Gates
- **Playwright E2E Suite**: 59 passed / 59 tests across 11 suites (`npx playwright test --workers=1`).
- **Vitest Unit Suite**: 80 passed / 80 tests (`npm --prefix webapp test`).
- **TypeScript & Vite Build**: Clean build without errors (`npm --prefix webapp run build`).
- **Python Ruff Checks**: Passed with 0 errors (`python -m ruff check pytransdatro api tests tools`).
- **ANCPI Numerical Benchmark**: Zero functional regression ($< 10^{-9}$ tolerance maintained).
