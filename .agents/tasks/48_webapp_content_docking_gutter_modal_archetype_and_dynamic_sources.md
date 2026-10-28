# Web App Content Standardization, Yahoo-Style Docking Boundary, Modal Archetype 6.2, Dynamic Sources, and File Option Layout

## Objectives
- [x] Standardize coordinate widget readout in `MapView.ts` to use `Lon:` instead of `Long:`, establishing universal abbreviation consistency (`Lat:` and `Lon:`) across the web app.
- [x] Replace protruding 20×48px rectangular tab toggle on the desktop transformation panel with a sleek Yahoo Mail / Linear-style full-height vertical boundary gutter (`.workbench-dock-toggle`), rendering a top-to-bottom accent highlight line on hover, custom high-contrast SVG single-directional arrow cursors (`←` to dock left, `→` to expand right), and complete elimination of visual center boxes/arrows.
- [x] Move the "Ignore header / first line" toggle on the Files tab to the right side of the card, swapping the internal order so the text label is on the left and the sliding toggle switch is on the right (`FileUpload.ts`, `forms.css`).
- [x] Filter map points in bulk file parsing (`FileUpload.ts`) and map rendering (`MapView.ts`) to valid successfully transformed coordinates only.
- [x] Suppress the 1,000-point high-volume warning modal completely when zero valid points exist (e.g. when 5,000 points fail with out-of-grid or syntax errors due to wrong coordinate system input).
- [x] Redesign `VolumeWarningModal.ts` according to Archetype 6.2 (Secondary & Utility Action Buttons): compact modal geometry (`max-width: 440px`), structured horizontal `.modal-footer`, and compact 30px buttons (`.btn-sm`: `#vw-skip-btn`, `#vw-sample-btn`, `#vw-all-btn`) with `Escape` key and backdrop dismissal.
- [x] Dynamically filter client sources on the Usage / Telemetry dashboard (Computation and Spatial sub-views) by querying `usage.db` in `get_telemetry_timeline` (`api/routes/info.py`) to only return sources with $>0$ recorded transactions, omitting sources (like `desktop_cli`) when they have zero records.
- [x] Synchronize living documentation in `docs/webapp/archetypes.md`, `docs/webapp/actionable_inventory.md`, `docs/rest_api.md`, and `webapp/README.md`.
- [x] Verify full regression suites: Vitest unit tests (83/83), Vite production build, Playwright E2E browser suite (71/71), and Python pytest suite (187/187).

---

## Implementation Details

1. **Coordinate Label Standardization (`webapp/src/components/MapView.ts`, `webapp/e2e/map.spec.ts`)**:
   - Replaced all 4 instances of `Long:` with `Lon:` in `updateCursorBadge()`.
   - Updated Playwright map tests in `webapp/e2e/map.spec.ts` to assert against `/^Lat:.*Lon:/` and `/^Lon:.*Lat:/`.

2. **Yahoo Mail & GitHub Primer Neutral Workbench Docking Boundary (`webapp/src/style/layout.css`, `webapp/src/main.ts`)**:
   - Transformed `.workbench-dock-toggle` into a full-height interactive boundary rail (`width: 13px`, `height: 100%`, `top: 0; bottom: 0; left: 473px;` when expanded, `left: 0; width: 14px;` when collapsed).
   - Engineered a dual-layer Yahoo Mail-inspired shifted guide rail strictly using GitHub Primer Neutral tokens from `variables.css`:
     - `::before` renders a subtle 5px tinted docking channel (`var(--color-primary-glow)`) anchored against the structural panel border (`var(--border-medium)`), while showing a subtle 2px docked lip along `left: 0` when collapsed so the docked edge remains discoverable.
     - `::after` renders a crisp `1.5px` top-to-bottom vertical guide line in `var(--border-focus)` (`#7c93ab` dark / `#57606a` light, brightening to `var(--text-main)` on `:active`) that smoothly shifts `4px` inward (`translateX(-4px)`) on hover when expanded and `4px` outward (`translateX(4px)`) on hover when collapsed.
   - Designed precision OS-grade filled single-direction SVG cursors (`fill="#ffffff" stroke="#0e1013" stroke-width="2.4" paint-order="stroke fill"` with `feDropShadow`):
     - Expanded state: Crisp solid-white left-pointing arrow (`◀▬`) with centered hotspot (`11 12`), indicating panel docking to the left.
     - Collapsed state: Crisp solid-white right-pointing arrow (`▬▶`) with left-anchored hotspot (`7 12` so it never clips off the left viewport edge), indicating panel expansion to the right.
   - Removed the visual rectangular pill in the middle entirely (`.workbench-dock-toggle-icon` is styled as `sr-only` for screen readers and Playwright tests without overlapping the line).
   - Preserved accessible keyboard shortcuts (`Ctrl+\` / `Alt+W`) and localStorage persistence.

3. **Files Tab "Ignore Header" Layout Alignment (`webapp/src/components/FileUpload.ts`, `webapp/src/style/forms.css`)**:
   - In `forms.css`: Updated `.file-options-row` to `justify-content: flex-end;`, positioning the option row on the right side of the card.
   - In `FileUpload.ts`: Restructured `.toggle-switch-container` children so the label `<span class="toggle-switch-label">` appears first on the left, followed by the adjacent `<input class="toggle-switch-input">` and `<span class="toggle-switch-track">` on the right.

4. **High-Volume Point Warning Filtering & Archetype 6.2 Redesign (`webapp/src/components/FileUpload.ts`, `webapp/src/components/MapView.ts`, `webapp/src/components/VolumeWarningModal.ts`, `webapp/src/style/layout.css`)**:
   - In `FileUpload.ts` (`parsePairedPoints`): Skipped error and out-of-grid output lines (`outRaw.includes('#') || isOutOfGrid`) from being converted into map items.
   - In `MapView.ts` (`plotBatchPoints`, `plotFilePoints`): Gated volume threshold check strictly on valid coordinates (`validPoints = points.filter(...)`). If `validPoints.length === 0`, clears the layer and immediately exits without showing the volume modal. If `validPoints.length > 1000`, displays the modal passing the count of valid points.
   - In `VolumeWarningModal.ts`: Restructured from vertically stacked 42px primary buttons in the body to a compact dialog (`.modal-dialog-compact`, `max-width: 440px`), explanatory body text, and horizontal `.modal-footer` housing compact Archetype 6.2 buttons (`.btn-sm`: Cancel `#vw-skip-btn`, Sample `#vw-sample-btn`, Render All `#vw-all-btn`).
   - Added automated Playwright test in `webapp/e2e/file.spec.ts` verifying that 2,050 invalid/out-of-grid points suppress the volume modal.

5. **Dynamic Telemetry Client Sources Filter (`api/routes/info.py`, `webapp/src/components/TelemetryDashboard.ts`, `tests/test_api_modern.py`)**:
   - In `api/routes/info.py` (`get_telemetry_timeline`): Replaced static `SELECT id, name, description FROM sources` with a query filtering `WHERE EXISTS` in `computation_usage` or `spatial_usage` with count $>0$.
   - In `webapp/src/components/TelemetryDashboard.ts`: Updated dynamic source resolution to ingest `res.available_sources` directly, rendering filter buttons in Computation and Spatial views only for sources with actual database records.
   - Added backend test in `tests/test_api_modern.py` verifying that sources with 0 records (e.g. `desktop_cli`) are omitted from `available_sources`.
   - Updated E2E assertions in `webapp/e2e/computation.spec.ts` and `webapp/e2e/spatial.spec.ts` to assert dynamic source count $\ge 1$.

---

## Verification & Testing

1. **Vitest Unit Test Suite**:
   - `npm --prefix webapp test`
   - Result: 83/83 unit tests passing (`336ms`).
2. **Vite Production Compilation**:
   - `npm --prefix webapp run build`
   - Result: Clean TypeScript compilation (`tsc -b`) and minified production build in `418ms`.
3. **Playwright E2E Browser Suite (71/71 Passed)**:
   - `powershell -ExecutionPolicy Bypass -File .\test-e2e.ps1`
   - Result: All 71 tests passed across 12 spec files (`1.1m`).
4. **Full Python Test & Regression Suite (187/187 Passed)**:
   - `pytest tests/test_trans_ro.py tests/test_api_modern.py tests/test_api_legacy.py tests/test_api_rate_limiter.py tests/test_telemetry.py tests/test_feedback.py tests/test_config.py`
   - Result: 187/187 tests passed (`15.43s`) with zero geodetic regressions ($< 10^{-9}$ tolerance).
5. **Ruff Code Hygiene**:
   - `ruff check --no-cache api/ pytransdatro/ tests/ tools/`: Clean (0 errors).
   - `ruff format --check --no-cache api/ pytransdatro/ tests/ tools/`: All 56 files formatted.
