# Telemetry Upgrade, Dynamic Windows & REST API Timeline

## Objectives
- [x] Restructure telemetry database schema to support dynamic time-window granularities and self-describing date boundaries (`period_start`, `period_end`).
- [x] Introduce normalized `operations` dimension table (`1: st70_to_etrs89`, `2: etrs89_to_st70`), replacing the previous `is_st70_to_etrs89` boolean convention.
- [x] Track vertical dimensionality (`is_3d`) across both volume and spatial distribution records.
- [x] Separate metrics into `computation_usage` (function calls & processed points volume) and `spatial_usage` (~1 km geographic cell hits).
- [x] Implement dynamic window registry table (`telemetry_config`) supporting Slowly Changing Dimension (SCD Type 2) validity intervals.
- [x] Swap client source IDs in `sources` table to set `1: rest_api` and `2: web_map`.
- [x] Implement automated database migration routine preserving legacy data from `daily_stats` and `usage` in `logs/usage.db`.
- [x] Connect `X-Client-Source` header resolution across all FastAPI transformation routes.
- [x] Implement the timeline visualization REST API endpoint (`GET /api/v1/telemetry/timeline`) per `specs/telemetry_visualization.md`.
- [x] Update documentation in `docs/logging.md` and `docs/rest_api.md`, and verify test suite pass rate.

---

## Implementation Details

### 1. Pure Python Core & Telemetry Subsystem (`pytransdatro/`)

#### Dynamic Window Registry & Schema Overhaul (`pytransdatro/telemetry.py`)
- **`sources` table**: Updated client IDs to standardize REST API as primary backend and Web Map as primary frontend:
  - `0`: `unknown` (Unspecified client)
  - `1`: `rest_api` (REST API Service / backend microservices)
  - `2`: `web_map` (Web Map Application)
  - `3`: `desktop_cli` (Desktop / CLI GIS applications)
- **`operations` dimension table**:
  - `1`: `st70_to_etrs89` (Stereo 70 -> ETRS89)
  - `2`: `etrs89_to_st70` (ETRS89 -> Stereo 70)
- **`telemetry_config` registry**:
  - Tracks `metric_type` (`computation` or `spatial`), `window_days`, `anchor_date`, `valid_from`, and `valid_to`.
  - Configured by default to 1-day windows for `computation` and 30-day windows for `spatial`.
- **`computation_usage` table** (formerly `daily_stats`):
  - Primary Key: `(period_start, operation, is_3d, source)`.
  - Tracks `count_calls` and `count_points` per period window.
- **`spatial_usage` table** (formerly `usage`):
  - Primary Key: `(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)`.
  - Stores `period_start` and `period_end` date strings, preserving exact provenance regardless of future window length modifications.
- **Self-Describing Snapping Window Algorithm**:
  - Automatically maps timestamps to `period_start` and `period_end` using anchor-aligned delta calculations (pure standard library `datetime` and `calendar`).
- **Seamless Database Migration**:
  - Automatically migrates existing legacy `daily_stats` and `usage` tables into the new schema without downtime or data loss, and migrates existing source IDs `1` and `2` if detected in old order.

#### Core Engine Parameter Alignment (`pytransdatro/trans_ro.py`)
- Updated telemetry collector calls to pass `operation = 1` for `st70_to_etrs89` and `operation = 2` for `etrs89_to_st70` along with `is_3d`.

---

### 2. REST API Web Layer (`api/`)

#### Client Origin Detection (`api/dependencies.py`)
- Introduced `get_client_source(request: Request) -> int` dependency:
  - Reads `X-Client-Source` HTTP header.
  - Returns `2` (`web_map`) when header is `web_map` or `web_app`.
  - Returns `3` (`desktop_cli`) when header is `desktop_cli`.
  - Defaults to `1` (`rest_api`) for standard API requests.
- Injected dependency into:
  - `api/routes/point.py`
  - `api/routes/batch.py`
  - `api/routes/file.py`
  - `api/routes/legacy.py`

#### Timeline Visualization Endpoint (`api/routes/info.py` & `api/schemas/info.py`)
- Fixed existing `GET /api/v1/telemetry/stats` to query `computation_usage`.
- Added `GET /api/v1/telemetry/timeline` implementing the specification in `specs/telemetry_visualization.md`:
  - **Query Parameters**: `from_date`, `to_date`, `granularity` (`daily` or `monthly`), `source`, `operation`, `is_3d`.
  - **Summary Metrics**: `total_calls`, `total_points`, `avg_points_per_call`, `from_date`, `to_date`, `days_count`.
  - **Time Series**: Aggregates calls and points per bucket with breakdowns for `web_map`, `rest_api`, `desktop_cli`, `unknown`, and `total`.

---

### 3. Documentation Updates (`docs/`)
- **`docs/logging.md`**: Fully refreshed to describe the 5-table schema, self-describing window intervals, operations dimension, and sample analytical SQL queries.
- **`docs/rest_api.md`**: Added endpoint documentation for `GET /api/v1/telemetry/timeline` and updated metadata descriptions.

---

## Verification & Testing

### Automated Test Suite
- Executed via `D:\anaconda3\envs\pytransdat\python.exe -m pytest`:
  - `tests/test_trans_ro.py`: 46 core geodetic mathematical tests.
  - `tests/test_api_legacy.py`: 11 legacy servlet compatibility tests.
  - `tests/test_api_rate_limiter.py`: 2 rate-limiter tests.
  - `tests/test_telemetry.py`: 6 telemetry tests covering queue flushing, operation tracking, 3D flags, window date ranges, out-of-grid failure isolation, and legacy schema migration.
  - `tests/test_api_modern.py`: 14 modern REST API tests covering transformation endpoints, file streaming, health, grid info, stats, and timeline visualization endpoints.
- **Total Pass Rate**: **79 passed**, 0 failed (100% success).

### Database Integrity Verification
- Verified migration on actual `logs/usage.db`:
  - `sources`: `[(0, 'unknown'), (1, 'rest_api'), (2, 'web_map'), (3, 'desktop_cli')]`.
  - `computation_usage`: 3,229 rows safely migrated.
  - `spatial_usage`: 99,385 rows safely migrated with calculated `period_start` and `period_end` date bounds.
