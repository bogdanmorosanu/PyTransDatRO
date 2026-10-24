# Telemetry Period Snapping, Sentinel Windows, Spatial Calls Parity, and Active Period UI

## Objectives
- [x] Investigate and resolve the mismatch between user-selected query date ranges and periodic database aggregations in telemetry tabs (`computation_usage` and `spatial_usage`).
- [x] Implement period snapping: database queries match period intervals using overlap (`c.period_end >= from_date AND c.period_start <= to_date`), and compute effective bounds `effective_from_date = MIN(period_start)` (snapped backwards) and `effective_to_date = MAX(period_end)` (snapped forwards).
- [x] Expose `effective_from_date` and `effective_to_date` in `/api/v1/telemetry/timeline` and `/api/v1/telemetry/spatial-matrix` endpoints.
- [x] Display "Active Period" badge (`effective_from_date` $\to$ `effective_to_date`) within the date limits cluster above both Computation and Spatial dashboards in the web application.
- [x] Support explicit calendar sentinel window values in `telemetry_config`:
  - `window_days > 0`: Mathematical $N$-day step.
  - `window_days = -1`: Calendar month (from day 1 to last day of month).
  - `window_days = -2`: Calendar year (from Jan 1 to Dec 31).
- [x] Add `count_calls INTEGER NOT NULL DEFAULT 1` to `spatial_usage` so spatial calls metric is `SUM(count_calls)` instead of the row-counting hack `COUNT(*)`.
- [x] Implement seamless transition helper `SqliteUsageLogger.update_window_config(...)` that clamps active period `period_end` to today and starts new periods tomorrow (zero overlap, zero data loss).
- [x] Rewrite `tools/populate_dummy_telemetry.py` to target modern tables across 3 distinct period regimes (Yearly, Monthly, 10-day) with 100% mathematical parity between computation points and spatial points.
- [x] Update and expand unit tests (`tests/test_telemetry.py`, `tests/test_api_modern.py`) and browser E2E specs (`webapp/e2e/computation.spec.ts`, `webapp/e2e/spatial.spec.ts`).

---

## Implementation Details

1. **Telemetry Core & Dynamic Window Sentinels (`pytransdatro/telemetry.py`)**:
   - Enhanced `_calculate_period_boundaries(cur_date, window_days, anchor_date)`:
     - When `window_days == -1` (month sentinel): determines month start `cur_date.replace(day=1)` and month end using `calendar.monthrange(cur_date.year, cur_date.month)[1]`.
     - When `window_days == -2` (year sentinel): determines year start `cur_date.replace(month=1, day=1)` and year end `cur_date.replace(month=12, day=31)`.
     - When `window_days == 1`: start and end are identical to `cur_date`.
     - When `window_days > 1`: steps $N$ days forward from `anchor_date` modulo $N$.
   - Schema auto-migration: Added check on `PRAGMA table_info(spatial_usage)` to dynamically add `count_calls INTEGER NOT NULL DEFAULT 1` if missing on older databases.
   - Batch Aggregation: Updated `_process_batch(...)` spatial dictionary to key by `(lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)` and store `{'count': points, 'count_calls': calls}`.
   - UPSERT logic: Updated to `DO UPDATE SET count = count + excluded.count, count_calls = count_calls + excluded.count_calls`.
   - Window Transition Helper: Implemented `update_window_config(metric_type, window_days)` which clamps active `telemetry_config.valid_to` and any active `period_end` to today, setting `valid_from = tomorrow` and `anchor_date = tomorrow` for the new configuration.

2. **Synthetic Data Generator with Mathematical Parity (`tools/populate_dummy_telemetry.py`)**:
   - Modernized script to populate `computation_usage`, `spatial_usage`, and `telemetry_config` directly.
   - Created three distinct period regimes spanning 3 years (~1,073 days):
     - Regime 1: Calendar Year (`window_days = -2`).
     - Regime 2: Calendar Month (`window_days = -1`).
     - Regime 3: 10-day Mathematical Step (`window_days = 10`).
   - Parity guarantee: Every single point logged to `computation_usage` is mapped directly to a 1 km spatial cell in `spatial_usage`, ensuring `SUM(c.count_points) == SUM(s.count)` across any time window.
   - Generated 6,565 computation records and 18,398 spatial records in 1.8s with residual discrepancy $= 0$.

3. **REST API & Schema Contracts (`api/schemas/info.py`, `api/routes/info.py`)**:
   - Extended `TimelineSummary` and `SpatialMatrixSummary` models with `effective_from_date: Optional[str] = None` and `effective_to_date: Optional[str] = None`.
   - In `/api/v1/telemetry/timeline`:
     - Changed SQL filter from `period_start BETWEEN :from_date AND :to_date` to interval overlap `c.period_end >= :from_date AND c.period_start <= :to_date`.
     - Added boundary query `SELECT MIN(period_start), MAX(period_end) ...` to extract `effective_from_date` and `effective_to_date`.
   - In `/api/v1/telemetry/spatial-matrix`:
     - Applied interval overlap filtering and extracted `effective_from_date` and `effective_to_date`.
     - Calculated calls as `SUM(s.count_calls)` instead of `COUNT(*)`, ensuring spatial call counts are invariant to bucket size.

4. **Web Application & UI Indicators (`webapp/`)**:
   - `webapp/src/utils/apiClient.ts`: Added `effective_from_date` and `effective_to_date` to `TimelineSummary` and `SpatialMatrixSummary` interfaces.
   - `webapp/src/i18n/en.ts` & `webapp/src/i18n/ro.ts`: Added dictionary key `telemetry.effectivePeriod`: `'Active Period'` / `'Perioad\u0103 activ\u0103'`.
   - `webapp/src/components/TelemetryDashboard.ts`: Appended third indicator subgroup to `.telemetry-date-limits-cluster`:
     `<span>${t('telemetry.effectivePeriod')}:</span> <strong>${summary.effective_from_date} &rarr; ${summary.effective_to_date}</strong>`.
   - Backwards compatibility: Maintained `nth(0)` ("Earliest Available") and `nth(1)` ("Current Date") positions so existing E2E tests and layout expectations remained stable.

---

## Verification & Testing

- **Core Geodetic Precision Tests**:
  - `pytest tests/test_trans_ro.py`: 54/54 passed (0 numerical regressions, $< 10^{-9}$ tolerance maintained).
- **Python Telemetry & API Tests**:
  - `pytest tests/test_telemetry.py`: 31/31 passed (including new sentinel tests, transition tests, and spatial call/point count tests).
  - `pytest tests/test_api_modern.py`: 10/10 passed (including `test_telemetry_effective_date_snapping_and_parity`).
- **Frontend Unit Tests (Vitest)**:
  - `npm --prefix webapp test`: 80/80 passed across 6 test suites.
- **Frontend TypeScript & Production Build**:
  - `npx --prefix webapp tsc --noEmit -p webapp/tsconfig.json`: 0 errors.
  - `npm --prefix webapp run build`: Minified production bundle generated cleanly.
- **Browser End-to-End Tests (Playwright)**:
  - `webapp/e2e/computation.spec.ts`: Passed (verified Active Period badge rendered as `nth(2)` with `YYYY-MM-DD → YYYY-MM-DD`).
  - `webapp/e2e/spatial.spec.ts`: Passed (verified Active Period badge rendered as `nth(2)` with `YYYY-MM-DD → YYYY-MM-DD`).
