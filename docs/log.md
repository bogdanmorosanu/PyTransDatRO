# PyTransDatRO: Application Logging (`logs/app.log`)

## 1. Overview & Purpose
`PyTransDatRO` provides a robust, standard-library-based application logging mechanism designed to track application lifecycle events, request processing issues, diagnostics, warnings, and fallback dumps when subsystems experience operational difficulties.

> [!NOTE]
> This document specifically describes **application file and console logging** (e.g. `logs/app.log`, logger hierarchies, configuration parameters, and error fallback captures). For **usage metrics & spatial density telemetry** persisted in SQLite (`logs/usage.db`), refer to [`docs/telemetry.md`](telemetry.md).

---

## 2. Architectural Principles & Layer Discipline
- **Pure Python Standard Library (`pytransdatro/`)**: Core mathematical processing uses Python's standard `logging` library (`logging.getLogger("pytransdatro")`). It does not import or depend on web frameworks or third-party logging packages.
- **Centralized Configuration (`api/config.py`)**: Environment variables configure application-wide logging through Python's standard `dictConfig` format.
- **Rotating File Logging**: Log entries are automatically written to rotating disk files (`logs/app.log`) with size caps and backup retention to prevent unbounded disk utilization.
- **Unified Handlers**: Handlers simultaneously direct formatted messages to standard output (`console`) and to disk (`file`), ensuring compatibility with Docker/container logs (`stdout`) as well as persistent volume mounts.

---

## 3. Configuration & Environment Variables

Logging is configured through `Settings` in `api/config.py` and supports `.env` overrides with the `PYTRANSDAT_` prefix:

| Setting Key | Environment Variable | Default Value | Description |
|---|---|---|---|
| `log_level` | `PYTRANSDAT_LOG_LEVEL` | `INFO` | Logging threshold level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `app_log_enabled` | `PYTRANSDAT_APP_LOG_ENABLED` | `true` | Enables or disables rotating file logging. |
| `app_log_path` | `PYTRANSDAT_APP_LOG_PATH` | `./logs/app.log` | Destination file path for application logs. |
| `app_log_max_bytes` | `PYTRANSDAT_APP_LOG_MAX_BYTES` | `10485760` (10 MB) | File size in bytes before rotating to a new archive file. |
| `app_log_backup_count` | `PYTRANSDAT_APP_LOG_BACKUP_COUNT` | `5` | Maximum number of rotated backup files to retain (`app.log.1`, `app.log.2`, etc.). |
| `feedback_log_path` | `PYTRANSDAT_FEEDBACK_LOG_PATH` | `./logs/feedback.log` | Dedicated append log for user feedback submissions. |

### Sample `.env` Configuration
```ini
# Logging Settings
PYTRANSDAT_LOG_LEVEL=INFO
PYTRANSDAT_APP_LOG_ENABLED=true
PYTRANSDAT_APP_LOG_PATH=./logs/app.log
PYTRANSDAT_APP_LOG_MAX_BYTES=10485760
PYTRANSDAT_APP_LOG_BACKUP_COUNT=5
```

---

## 4. Log Formatting & Logger Hierarchy

### 4.1 Format Specification
The standard formatter produces structured, human-readable entries with exact timestamps, log level, and logger name:
```
%(asctime)s [%(levelname)s] %(name)s: %(message)s
```
*Date format*: `YYYY-MM-DD HH:MM:SS`

*Example output*:
```text
2026-10-21 01:40:35 [INFO] api.main: PyTransDatRO API initialized successfully.
2026-10-21 01:40:35 [WARNING] pytransdatro.telemetry: Telemetry queue saturated (size=10000/10000). Dropped 1 points for operation 1 (source=1).
2026-10-21 01:40:35 [ERROR] pytransdatro.telemetry: SQLite write failed after 3 attempts: database is locked. Logging unwritten batch payload (1 comp, 1 spatial items) to application log.
```

### 4.2 Logger Namespace Hierarchy
- `api`: Root logger for FastAPI endpoints, middleware, dependency injectors, rate limiters, and routers.
- `pytransdatro`: Logger for the core geodetic engine, coordinate parsers, and grid loaders.
- `pytransdatro.telemetry`: Logger for background SQLite usage telemetry processing, queue monitoring, and error diagnostics.
- `root`: Default fallback logger capturing unhandled application and third-party library events.

---

## 5. Integration with Telemetry & Error Fallbacks

One of the primary roles of the application log is serving as a reliable failsafe when structured database writes or telemetry components encounter errors:

1. **Queue Saturation Warnings**:
   When client traffic floods the in-memory queue beyond capacity, `SqliteUsageLogger` emits a `WARNING` entry detailing the count of dropped coordinates, operation code, and client source.
2. **Database Write Retries & Fallback Payloads**:
   If `usage.db` is locked or disk errors occur, the worker thread retries up to 3 times with exponential backoff. If writing still fails, it emits an `ERROR` containing the full unwritten computation and spatial batch payloads directly into `logs/app.log`. This ensures that coordinate metrics can be recovered and reviewed even during database downtime.
3. **Worker Thread Health & Diagnostics**:
   Any unexpected thread exceptions update `last_error` and are logged with complete tracebacks without crashing the main application.

---

## 6. Verification & Automated Testing

The logging system is validated via automated tests in [`tests/test_logging.py`](../tests/test_logging.py):
- `test_rotating_file_handler_creation_and_formatting`: Verifies `configure_logging()` creates `logs/app.log` and formats messages across multiple logger hierarchies.
- `test_log_file_rotation`: Emits messages exceeding `max_bytes` and asserts that backup files (`app.log.1`, `app.log.2`) are generated without exceeding `backup_count`.
- `test_telemetry_queue_saturation_logged_to_file`: Verifies that queue drops write structured warning entries into the file log.
- `test_telemetry_db_write_failure_fallback_logging`: Simulates SQLite operational failures and asserts that unwritten batches are logged to `app.log`.
