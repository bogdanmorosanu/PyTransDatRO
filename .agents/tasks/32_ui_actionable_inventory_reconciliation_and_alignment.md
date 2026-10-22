# Task 32: UI Actionable Elements Inventory Reconciliation & Implementation

## Objectives
- [x] Compare living specification in `docs/ui_actionable_inventory.md` against codebase implementation across `webapp/src/` and `webapp/e2e/`.
- [x] Classify divergences and gaps into UI/UX Recommended vs. Not Recommended changes based on surveying workflows and ergonomic best practices.
- [x] Implement URL deep-linking parameter hydration (`?op=`, `?unit=`, `?order=`, `?c1=`, `?c2=`, `?c3=`) across top-level application state (`main.ts`) and child converters.
- [x] Implement batch textarea input locking (`disabled readonly`) during active calculation in `BatchConverter.ts`.
- [x] Extend batch quality report badges with Total Points, In Grid, Out of Grid, No Data on Grid, Syntax Errors, and 2D vs. 3D breakdown counts.
- [x] Add `.dat` file format extension support to native file picker in `FileUpload.ts`.
- [x] Centralize `batchSampleSize` alongside `batchWarningThreshold` in `webapp/src/config.ts`.
- [x] Introduce locale-aware telemetry number formatting (`numberFormat.ts`) respecting Romanian (`,`) and English (`.`) decimal rules while preserving strict dot-decimal geodetic coordinate purity.
- [x] Eliminate asynchronous race conditions in `PointConverter.ts` using `pendingCalculation` queue.
- [x] Add Playwright E2E test accompaniment for deep-linking, batch badge breakdown, measurement cancel on Escape, and locale-aware telemetry formatting.
- [x] Pass 100% of Vitest unit tests (75/75), Playwright E2E browser tests (24/24), Python core geodetic tests (54/54), REST API tests (127/127), and ANCPI official benchmarks.

## Implementation Details

1. **Centralized Configuration (`webapp/src/config.ts`, `MapView.ts`, `VolumeWarningModal.ts`)**:
   - Added `batchSampleSize: number` to `AppConfig` and `APP_CONFIG` (defaults to `1,000`, overridable via `VITE_BATCH_SAMPLE_SIZE`).
   - Replaced hardcoded `1000` slices in `MapView.plotBatchPoints` and `MapView.plotFilePoints` with `APP_CONFIG.batchSampleSize`.
   - Updated `VolumeWarningModal.ts` button label to dynamically interpolate `{sample}` from `APP_CONFIG.batchSampleSize`.

2. **Global URL Deep-Linking Hydration (`webapp/src/main.ts`, `PointConverter.ts`)**:
   - Implemented `readInitialQueryParams()` in `App` (`main.ts`) to parse `op`, `unit`, `order`, `c1`, `c2`, `c3` from `window.location.search`.
   - Synchronized `currentOp`, `currentAngularUnit`, and `currentOrder` on application boot, updating segmented sliders, dropdowns, and child component states.
   - Enhanced `PointConverter.ts` to parse `unit` and `order` and maintain synchronized browser history replacement via `updateQueryParams()`.
   - Added `pendingCalculation` queue mechanism in `PointConverter.ts` to prevent dropped calculations when direction switches or rapid inputs race existing requests.

3. **Batch Converter Input Locking & Extended Quality Report (`webapp/src/components/BatchConverter.ts`)**:
   - Bound `disabled readonly` attributes to `#batch-input-area` and `disabled` to `#batch-transform-btn` during active API transit (`isProcessing = true`).
   - Extended `lastSummary` to track `totalCount`, `inGridCount`, `outOfGridCount`, `noDataCount`, `syntaxErrorCount`, `count2D`, and `count3D`.
   - Rendered standard quality badges matching Single Point mode (`Total points: X`, `✓ In Grid: Y`, `⚠ Out of Grid: Z`, `⚠ No Grid Data: W`, `✕ Syntax Errors: V`, `X 2D / Y 3D`).

4. **File Upload Mode Extension Support (`webapp/src/components/FileUpload.ts`)**:
   - Added `.dat` extension to hidden file picker: `accept=".csv,.txt,.xyz,.dat"`.

5. **Locale-Aware Telemetry Formatting (`webapp/src/utils/numberFormat.ts`)**:
   - Authored `formatTelemetryNumber(num, options)` and `formatTelemetryCompact(num, fractionDigits)`.
   - Enforces `ro-RO` (comma decimal) when active language is Romanian and `en-US` (dot decimal) when English.
   - Integrated into `TimelineChart.ts`, `LegacyTimelineChart.ts`, `SpatialMatrixChart.ts`, and `TelemetryDashboard.ts`.
   - Added Vitest suite `webapp/src/utils/numberFormat.test.ts` (7 unit tests).
   - Guarded `typeof localStorage !== 'undefined'` in `webapp/src/i18n/index.ts` for safe SSR/testing environments.

6. **Authoritative Living Inventory Reconciliation (`docs/ui_actionable_inventory.md`)**:
   - Updated Section 3.3 to classify smart-paste auto-calculation as an approved convenience feature.
   - Updated Section 3.4 to document textarea locking and clarify that the volume warning dialog guards Leaflet DOM rendering rather than blocking text math.
   - Documented `webapp/src/config.ts` configuration parameters in Section 3.5.

7. **Automated E2E Test Accompaniment (`webapp/e2e/`)**:
   - `point.spec.ts`: Added test asserting full query parameter deep-linking hydration (`?op=Stereo70ToETRS89&unit=dms&order=EN&c1=500000&c2=500000&c3=100`).
   - `batch.spec.ts`: Added test verifying extended quality badges and 2D/3D count breakdown pills.
   - `i18n.spec.ts`: Added test verifying telemetry formatting while asserting coordinate inputs strictly preserve dot-decimal numbers.
   - `map.spec.ts`: Added test verifying distance measurement tool deactivates on `Escape` key.

## Verification & Testing
- **Vitest Unit Suite**: 75/75 passed across 6 test files (`npm --prefix webapp test`).
- **Production Build**: Clean bundle compilation with zero TypeScript errors (`npm --prefix webapp run build`).
- **Playwright E2E Suite**: 24/24 passed headless via Chromium in 28.4s (`.\test-e2e.ps1`).
- **Core Geodetic Suite**: 54/54 passed in 3.50s (`pytest tests/test_trans_ro.py`).
- **REST API Suite**: 127/127 passed in 14.06s (`pytest tests/test_api_*.py`).
- **ANCPI Benchmarks**: Byte-for-byte fidelity confirmed via `tools/generate_ancpi_comparison_reports.py`.
