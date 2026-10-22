# Task 33: UI Actionable Inventory Audit Round 2 & Alignment Implementation

## Objectives
- [x] Conduct a complete second-round audit evaluating all 100 actionable elements in `docs/ui_actionable_inventory.md` against the web application implementation.
- [x] Reconcile the 4 specific user remarks regarding batch reporting layout, Romanian number/percentage formatting, scrollbar cursor, and map coordinate order.
- [x] Refactor batch quality reporting (`BatchConverter.ts`) into a two-tier layout:
  - Top row: `Total points: [val]`, `2D Coordinates: [val]`, `3D Coordinates: [val]`, `✓ In Grid: [val]`.
  - Bottom row: `⚠ Outside Grid: [val]`, `⚠ No Grid Data: [val]`, `✕ Syntax Errors: [val]` (on their own row underneath, only when diagnostic issues occur).
- [x] Implement locale-aware decimal and percentage formatters (`formatTelemetryDecimal`, `formatTelemetryPercent` in `numberFormat.ts`) to ensure quadrant share percentages (e.g. `(24,5%)`) and batch intensity ratios use comma in Romanian while preserving CAD/GIS dot-decimal coordinates.
- [x] Fix mouse cursor over `#batch-input-area` scrollbar by enforcing `cursor: default` on `::-webkit-scrollbar*` pseudo-elements in `webapp/src/style/base.css`.
- [x] Extend `MapView.ts` with `setCoordinateOrder` and coordinate order awareness so `.map-coords-badge` displays `Long: ...  Lat: ...` when `EN` is active and `Lat: ...  Long: ...` when `NE` is active.
- [x] Attach client-side SPA navigation handler to `.header-brand-link` in `Header.ts` to switch to Studio tab without full browser document reloads.
- [x] Ensure 100% DOM selector parity with living test specifications (`#tab-point`, `#tab-text`, `#tab-files`, `#angle-format-example-badge`, `#legacy-toggle-total`).
- [x] Pass all automated test suites: Vitest (79/79), Playwright E2E browser tests (25/25), Python core geodetic tests (62/62), Ruff linter, and official ANCPI benchmarks.

## Implementation Details

1. **Two-Tier Batch Quality Report (`webapp/src/components/BatchConverter.ts`)**:
   - Structured `.quality-report-card` into `.quality-report-row.stats-row` and `.quality-report-row.diagnostics-row` (separated by a subtle border).
   - Replaced combined `X 2D / Y 3D` pill with distinct labels: `2D Coordinates: [N]` and `3D Coordinates: [N]`.
   - Added translation keys `quality.coords2D` and `quality.coords3D` in `webapp/src/i18n/en.ts` and `webapp/src/i18n/ro.ts`.

2. **Locale-Aware Decimal & Percentage Telemetry Formatting (`webapp/src/utils/numberFormat.ts`)**:
   - Implemented `formatTelemetryDecimal(num, fractionDigits = 1)` and `formatTelemetryPercent(num, fractionDigits = 1)`.
   - Formatted quadrant cards in `TelemetryDashboard.ts` using `formatTelemetryPercent(q.percentage, 1)` (e.g., `(24,5%)` in RO vs. `(24.5%)` in EN).
   - Formatted `avgRatio` (Batch Intensity) in `TelemetryDashboard.ts` and `StatsView.ts` with `formatTelemetryDecimal(ratio, 1)`.
   - Applied locale formatting to file sizes in `FileUpload.ts` and comparison point counts in `GeodeticHelpPanel.ts`.
   - Extended `numberFormat.test.ts` to 11 Vitest unit tests verifying decimal and percentage conversions across locales.

3. **Custom Scrollbar Hover Cursor Fix (`webapp/src/style/base.css`)**:
   - Added `cursor: default;` to `::-webkit-scrollbar`, `::-webkit-scrollbar-track`, and `::-webkit-scrollbar-thumb`.
   - Resolved Chromium inheritance where editable textarea scrollbars showed the text I-beam cursor instead of the standard arrow.

4. **Map Coordinates Badge Coordinate Order Synchronization (`webapp/src/components/MapView.ts`, `main.ts`)**:
   - Added `currentOrder: CoordinateOrder = 'NE'` and `setCoordinateOrder(order)` to `MapView.ts`.
   - Formatted `.map-coords-badge` text as `Long: ${lonStr}  Lat: ${latStr}` in `EN` mode and `Lat: ${latStr}  Long: ${lonStr}` in `NE` mode.
   - Wired `this.mapView.setCoordinateOrder(...)` during initial boot and coordinate order toggle events in `main.ts`.

5. **Brand Logo Client SPA Reset (`webapp/src/components/Header.ts`)**:
   - Intercepted click on `.header-brand-link` to prevent native browser page reload and trigger `callbacks.onNavigate('studio')`.

6. **DOM Selector Parity & Legacy Total Toggle (`main.ts`, `TelemetryDashboard.ts`)**:
   - Added `id="tab-point"`, `id="tab-text"`, `id="tab-files"` to the transformation mode tab buttons.
   - Added `id="angle-format-example-badge"` to the angle format example chip.
   - Added `#legacy-toggle-total` interactive toggle button in `TelemetryDashboard.ts` to allow users to toggle the cumulative total series line on/off.

7. **Automated Verification & Living Documentation**:
   - Updated `webapp/e2e/i18n.spec.ts` asserting Romanian comma formatting for quadrant percentages and SPA brand logo navigation.
   - Updated `webapp/e2e/batch.spec.ts` asserting the two-tier quality report structure and separate 2D/3D labels.
   - Updated `webapp/e2e/map.spec.ts` asserting `.map-coords-badge` reflects `NE` vs. `EN` coordinate order.
   - Updated `webapp/e2e/legacy.spec.ts` asserting `#legacy-toggle-total` visibility and toggle states.
   - Synchronized `docs/ui_actionable_inventory.md` and `webapp/CHANGELOG.md`.

## Verification Results
- **Vitest Unit Suite**: 79 passed across 6 test files (100%).
- **Playwright E2E Suite**: 25 passed across 7 specs in 29.1s (100%).
- **Python Core Geodetic Suite**: 62 passed in 3.54s (100%).
- **Ruff Code Linter**: Zero errors.
- **Production Bundle**: Built cleanly in 310ms (`dist/assets/index-BLuVOKfS.js` 444 kB).
