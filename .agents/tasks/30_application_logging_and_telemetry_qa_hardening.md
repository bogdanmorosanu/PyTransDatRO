# Application Logging Infrastructure & Comprehensive Telemetry QA Test Suite

## Objectives
- [x] Rename `docs/logging.md` to `docs/telemetry.md` (specifically focused on usage statistics and spatial density SQLite telemetry) and create a dedicated `docs/log.md` detailing general application logging, log rotation, and failure dumps.
- [x] Configure standard library rotating file logging (`logs/app.log`) in `api/config.py` with environment-driven settings (`PYTRANSDAT_APP_LOG_*`), 10 MB file caps, 5 backups, and handlers for `api`, `pytransdatro`, and root loggers.
- [x] Harden `SqliteUsageLogger` in `pytransdatro/telemetry.py`:
  - Preserve pure Python standard library purity (zero third-party dependencies).
  - Implement retry logic with exponential backoff (up to 3 attempts) for transient SQLite errors (e.g. database locks).
  - Add fallback error logging of unwritten computation and spatial payloads to `logs/app.log` upon persistent write failure to avoid silent data loss.
  - Implement atomic dropped event counter and warning log on `queue.Full`.
  - Validate coordinate inputs in `_process_batch()` against non-finite floats (`NaN`, `Inf`) to ensure worker threads never terminate unexpectedly.
  - Add `get_diagnostics()` exposing worker thread status, queue depth, dropped counts, total points, and last error.
- [x] Optimize streaming file upload telemetry in `api/routes/file.py` with `BufferedTelemetryCollector` chunking points into 500-point batches, preventing queue saturation during 50,000+ line transformations.
- [x] Enhance `/health` endpoint in `api/routes/health.py` to monitor telemetry worker diagnostics, reporting degraded states when telemetry fails.
- [x] Implement comprehensive automated test coverage covering all 31 scenarios across Categories A through H in `tests/test_logging.py`, `tests/test_telemetry.py`, and `tests/test_api_modern.py`.

## Implementation Details

1. **Documentation Separation & Disambiguation (`docs/telemetry.md`, `docs/log.md`)**:
   - Renamed `docs/logging.md` to `docs/telemetry.md` to clarify that it exclusively covers usage metrics and spatial density binned telemetry in `logs/usage.db`.
   - Authored new `docs/log.md` documenting the application-wide rotating file logging infrastructure (`logs/app.log`), log formatting, hierarchy, environment variables, write failure fallbacks, and test accompaniment.
   - Updated cross-references across `.agents/skills/`, `CLAUDE.md`, and living architecture documentation.

2. **Rotating File Logging Subsystem (`api/config.py`, `.env.example`)**:
   - Added settings: `app_log_enabled` (default `True`), `app_log_path` (default `./logs/app.log`), `app_log_max_bytes` (10 MB), `app_log_backup_count` (5).
   - Configured `logging.handlers.RotatingFileHandler` alongside `StreamHandler` in `get_logging_config()`.

3. **Telemetry Collector Hardening & Fallback (`pytransdatro/telemetry.py`)**:
   - Added atomic dropped counters and last-error tracking.
   - Added 3-attempt exponential retry loop in `_process_batch()` for `sqlite3.OperationalError`.
   - Added unwritten payload dump to `logger.error(...)` if writes fail after all attempts.
   - Added `math.isfinite()` sanity checks in worker loop to isolate corrupt or non-finite inputs.
   - Implemented `BufferedTelemetryCollector` with slicing loops to guarantee chunks stay strictly $\le chunk\_size$ (500).

4. **Streaming File Upload Ingestion (`api/routes/file.py`)**:
   - Integrated `BufferedTelemetryCollector` into `line_generator`, wrapping transformations and ensuring all buffered points are flushed upon generator exit.

5. **Health Check Observability (`api/routes/health.py`)**:
   - Query `telemetry.get_diagnostics()`. Transition status to `"degraded"` with diagnostic error details if worker thread terminates or database write failures occur.

## Verification & Testing

- **New Logging Test Suite (`tests/test_logging.py`)**:
  - `test_rotating_file_handler_creation_and_formatting`: PASSED
  - `test_log_file_rotation`: PASSED
  - `test_telemetry_queue_saturation_logged_to_file`: PASSED
  - `test_telemetry_db_write_failure_fallback_logging`: PASSED
- **Telemetry Robustness Suite (`tests/test_telemetry.py`)**:
  - 28/28 tests passed covering Categories A–E, G.
- **API Integration Suite (`tests/test_api_modern.py`)**:
  - 37/37 tests passed including streaming file 500-point chunking and health probe degraded states.
- **Core Geodetic Suite (`tests/test_trans_ro.py`)**:
  - 54/54 tests passed with zero numerical regression against ANCPI benchmarks.
- **Full Python Test Suite**:
  - 171 passed, 18 skipped, 0 failed in 13.64s.
- **Frontend Vitest Suite**:
  - 68/68 passed in 314ms.
- **Ruff Code Audit**:
  - 0 errors, 42 files cleanly formatted.
