# Legacy TransdatOnline Statistics Integration & Usage Subpage

## Objectives
- [x] Audit data integrity of `legacy_stats.db` to check for missing periods and potential replacement by zero counts.
- [x] Diagnose and resolve MySQL race conditions that created duplicate calendar day records in the legacy dump.
- [x] Update ETL script `tools/extract_legacy_stats.py` to aggregate counts per day, enforce unique constraints, and read directly from `logs/localhost.sql`.
- [x] Introduce FastAPI REST API endpoints `/api/v1/telemetry/legacy/stats` and `/api/v1/telemetry/legacy/timeline` with SQLite daily/weekly/monthly grouping, date bounds slicing, and faceted source/operation filtering.
- [x] Integrate "Legacy" subpage tab into the Web App Usage view replicating the "Computation" tab styling and behaviors:
  - Custom SVG `LegacyTimelineChart` with smooth cubic bezier curve interpolation, gridlines, axis text, interactive crosshair, and tooltips.
  - Interactive legend swatches with color picker popovers and `localStorage` persistence.
  - Real-time dark/light theme switching support via `themechange` event dispatching.
  - Date range presets relative to database archive end date (`2026-10-14`).
  - Clear label copy: "Archive Earliest" / "Cea Mai Veche Dată".
- [x] Full automated test accompaniment across unit tests (`pytest`, `vitest`), TypeScript build, and Playwright browser E2E test.

## Implementation Details

1. **Integrity Audit Tool (`tools/audit_legacy_stats_integrity.py`)**:
   - Analyzed timeline across 5,758 calendar days (2011-01-09 to 2026-10-14).
   - Produced 3 structured reports in `reports/legacy_stats_integrity_report.txt`: missing date intervals (781 lost days across 322 gaps), zero-count days (20 days across 9 periods), and combined void days (801 days across 330 periods).
   - Confirmed no multi-year outage exists (longest gap was 56 days from Dec 2019 to Feb 2020).

2. **Deduplication & Extraction (`tools/extract_legacy_stats.py`)**:
   - Identified root cause in legacy Java server (`DBtransDatGrid.java`): lack of synchronization/locking caused TOCTOU race conditions where concurrent requests executed `SELECT` simultaneously and followed with separate `INSERT` statements.
   - Migrated dump input path to `logs/localhost.sql`.
   - Aggregated MySQL rows by calendar date (`day`), summing all 4 operation metrics (`cooOpAppSt70`, `cooOpAppSt30`, `cooOpServiceSt70`, `cooOpServiceSt30`).
   - Hardened SQLite schema with `day TEXT NOT NULL UNIQUE` and post-extraction validation (`GROUP BY day HAVING COUNT(*) >= 2`).

3. **Backend API Endpoints (`api/config.py`, `api/dependencies.py`, `api/schemas/info.py`, `api/routes/info.py`)**:
   - `Settings.legacy_db_path`: configurable path defaulting to `./logs/legacy_stats.db`.
   - `GET /api/v1/telemetry/legacy/stats`: aggregate totals (27,745,527 transformed points).
   - `GET /api/v1/telemetry/legacy/timeline`: time-series endpoint with daily, weekly, and monthly bucket expressions (`date(day, 'weekday 0', '-6 days')`, `substr(day, 1, 7)`), source filtering, and operation filtering.

4. **Frontend Architecture (`webapp/`)**:
   - `webapp/src/utils/apiClient.ts`: TypeScript interfaces and fetch helpers `getLegacyStats()` and `getLegacyTimeline()`.
   - `webapp/src/i18n/en.ts` & `webapp/src/i18n/ro.ts`: Full bilingual translation support.
   - `webapp/src/components/LegacyTimelineChart.ts`: Standalone SVG chart component matching `TimelineChart.ts` visual tokens, supporting points-only mode, custom swatches, and theme responsiveness.
   - `webapp/src/components/TelemetryDashboard.ts`: Integrated third subpage navigation tab "Legacy" under "Usage", complete with presets, granularity toggles, calendar inputs clamped to database bounds, and source/operation faceted buttons.
   - `webapp/src/components/Header.ts`: Emits `themechange` event upon theme toggling to notify telemetry charts.

## Verification & Testing
- **Backend Tests**:
  - `tests/test_api_legacy_stats.py`: 4 tests passed verifying summary metrics, daily timeline, monthly grouping, and source/op filters.
  - `tests/test_trans_ro.py`: 54 passed verifying zero regression in core geodetic precision.
- **Frontend Unit Tests**:
  - `npm --prefix webapp test`: 68 Vitest tests passed across 5 test suites.
  - `npm --prefix webapp run build`: TypeScript compiler and Vite production build succeeded with zero errors.
- **Browser End-to-End Suite**:
  - `webapp/e2e/legacy.spec.ts`: Playwright test passed, validating tab navigation, KPI cards, archive date boundaries, faceted toggles, and chart rendering.
