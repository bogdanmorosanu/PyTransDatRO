# PyTransDatRO: Application Logging (`logs/app.log`)

## 1. Overview & Purpose
`PyTransDatRO` provides a robust, standard-library-based application logging mechanism designed to track application lifecycle events, request processing issues, diagnostics, warnings, and fallback recovery when subsystems experience operational difficulties.

> [!NOTE]
> This document specifically describes **application file and console logging** (e.g. `logs/app.log`, logger hierarchies, configuration parameters, rate limiting alerts, and error captures). For **usage metrics & spatial density telemetry** persisted in SQLite (`logs/usage.db`), refer to [`docs/telemetry.md`](telemetry.md).

---

## 2. Architectural Principles & Layer Discipline

### 2.1 The "Errors & Warnings Only" Philosophy
To maintain high performance and avoid log clutter:
- **Silence on Healthy Operations**: Normal HTTP `200 OK` requests generate **zero** log entries.
- **Zero Coordinate Math Overhead**: Mathematical algorithms in `pytransdatro/` execute without per-point logging overhead.
- **High-Signal Diagnostics**: Logging is reserved strictly for operational anomalies (rate limiting, file validation errors, unhandled 500 exceptions, and subsystem delivery failures).

### 2.2 Symmetrical Status vs. Data Separation
Operational diagnostics and machine-readable data payloads are kept strictly segregated across three dedicated files:

```
                            +-------------------------------------------+
                            |               logs/app.log                |
                            |  (Operational Health & Error Log Only)    |
                            |  • 500 Unhandled Server Exceptions        |
                            |  • 429 Rate Limiter Throttling Warnings   |
                            |  • Streaming File Upload Parse Failures   |
                            |  • 1-line Telemetry DB Failure Warnings   |
                            |  • 1-line SMTP Email Delivery Warnings    |
                            +-------------------------------------------+
                                           |             |
                      If email fails:      |             | If SQLite fails:
                                           v             v
             +-------------------------------+         +-------------------------------------+
             |      logs/feedback.jsonl      |         |    logs/telemetry_fallback.jsonl    |
             |  (Persistent User Submissions)|         |   (Unwritten Aggregated DB Batches) |
             +-------------------------------+         +-------------------------------------+
```

- **`logs/app.log`**: Human-readable operational narrative. Never contains raw coordinate arrays or database parameter dumps.
- **`logs/feedback.jsonl`**: Machine-readable JSONL audit archive of user feedback submissions.
- **`logs/telemetry_fallback.jsonl`**: Emergency JSONL spillover containing unwritten aggregated telemetry batches if SQLite is locked or unreachable.

### 2.3 Pure Python Invariant
Per architectural rules in `GEMINI.md`, the logging system relies exclusively on Python's built-in standard library (`logging`, `logging.handlers`, `json`, `pathlib`). No third-party logging frameworks are used.

---

## 3. Configuration & Environment Variables

Logging is configured through `Settings` in `api/config.py` and supports `.env` overrides with the `PYTRANSDAT_` prefix:

| Setting Key | Environment Variable | Default Value | Description |
|---|---|---|---|
| `log_level` | `PYTRANSDAT_LOG_LEVEL` | `INFO` | Logging threshold level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `app_log_enabled` | `PYTRANSDAT_APP_LOG_ENABLED` | `true` | Enables or disables rotating file logging. |
| `app_log_path` | `PYTRANSDAT_APP_LOG_PATH` | `./logs/app.log` | Destination file path for application logs. |
| `app_log_max_bytes` | `PYTRANSDAT_APP_LOG_MAX_BYTES` | `10485760` (10 MB) | File size in bytes before rotating to a new archive file. |
| `app_log_backup_count` | `PYTRANSDAT_APP_LOG_BACKUP_COUNT` | `5` | Maximum number of rotated backup files to retain (`app.log.1` - `app.log.5`). |
| `feedback_log_path` | `PYTRANSDAT_FEEDBACK_LOG_PATH` | `./logs/feedback.jsonl` | Dedicated append log for user feedback submissions. |
| `telemetry_fallback_path` | `PYTRANSDAT_TELEMETRY_FALLBACK_PATH` | `./logs/telemetry_fallback.jsonl` | Dedicated emergency JSONL spillover for unwritten telemetry batches. |

### Sample `.env` Configuration
```ini
# Logging & Storage Settings
PYTRANSDAT_LOG_LEVEL=INFO
PYTRANSDAT_APP_LOG_ENABLED=true
PYTRANSDAT_APP_LOG_PATH=./logs/app.log
PYTRANSDAT_APP_LOG_MAX_BYTES=10485760
PYTRANSDAT_APP_LOG_BACKUP_COUNT=5
PYTRANSDAT_FEEDBACK_LOG_PATH=./logs/feedback.jsonl
PYTRANSDAT_TELEMETRY_FALLBACK_PATH=./logs/telemetry_fallback.jsonl
```

---

## 4. Log Formatting & Logger Hierarchy

### 4.1 Format Specification
The standard formatter produces structured, human-readable entries with exact timestamps, log level, and logger name:
```
%(asctime)s [%(levelname)s] %(name)s: %(message)s
```
*Date format*: `YYYY-MM-DD HH:MM:SS`

### 4.2 Logger Namespace Hierarchy
- `api.main`: Worker initialization, SPG grid memory pre-warming, worker shutdown, and unhandled 500 server crashes.
- `api.middleware.rate_limit`: Throttling warnings when clients exceed points/batch rate limits (HTTP 429).
- `api.routes.file`: Diagnostic warnings for rejected file uploads (invalid delimiter, file size/line limit exceeded, corrupt lines).
- `api.services.mailer`: Warnings when SMTP email delivery fails (preserving user payload in `logs/feedback.jsonl`).
- `pytransdatro`: Core geodetic engine milestones.
- `pytransdatro.telemetry`: Queue saturation warnings and SQLite write failure alerts (routing unwritten batches to `logs/telemetry_fallback.jsonl`).
- `root`: Default fallback logger capturing unhandled application and third-party library events.

### 4.3 Concrete Output Examples

#### Server Boot & Pre-Warming
```text
2026-10-22 20:00:00 [INFO] api.main: Initializing PyTransDatRO Web Service (Log Level: INFO)
2026-10-22 20:00:00 [INFO] api.main: Loading SPG grid into RAM...
2026-10-22 20:00:01 [INFO] api.main: TransRO transformation engine pre-warmed and ready.
```

#### Unhandled Server Crash (HTTP 500)
```text
2026-10-22 20:05:12 [ERROR] api.main: Unhandled server exception on POST /api/v1/transform/point: Simulated unexpected database crash
Traceback (most recent call last):
  File "api/routes/point.py", line 23, in transform_point
    out_coos, warn = execute_point_transformation(...)
RuntimeError: Simulated unexpected database crash
```

#### Rate Limiting Throttling (HTTP 429)
```text
2026-10-22 20:10:45 [WARNING] api.middleware.rate_limit: Rate limit exceeded for client 192.168.1.100 on POST /api/v1/transform/point (limit: 60 req/60s, retry_after: 42s)
```

#### File Upload Error Diagnostics
```text
2026-10-22 20:15:30 [WARNING] api.routes.file: File upload 'parcels.csv' issue at line 14: Invalid coordinate data
2026-10-22 20:16:01 [WARNING] api.routes.file: File upload 'survey.txt' rejected: invalid delimiter 'pipe'
```

#### SMTP Feedback Delivery Failure
```text
2026-10-22 20:20:15 [WARNING] api.services.mailer: SMTP feedback delivery failed (SMTPConnectError: Connection refused). Feedback safely preserved in logs/feedback.jsonl.
```

#### Telemetry Database Lock & Spillover
```text
2026-10-22 20:25:00 [ERROR] pytransdatro.telemetry: SQLite write failed after 3 attempts: database is locked. Saved unwritten batch (1 comp, 4 spatial items) to logs/telemetry_fallback.jsonl.
```

---

## 5. File Rotation & Disk Safety

To ensure bounded disk consumption:
- Uses `logging.handlers.RotatingFileHandler`.
- Rotates when file size reaches `app_log_max_bytes` (10 MB).
- Retains at most `app_log_backup_count` (5 backups: `app.log.1` through `app.log.5`).
- **Maximum Lifetime Footprint**: $10\text{ MB} \times (1 + 5) = \mathbf{60\text{ MB}}$ maximum total disk space. Oldest archives are purged automatically.

---

## 6. Verification & Automated Testing

The logging system is validated via 10 automated test scenarios in [`tests/test_logging.py`](../tests/test_logging.py):
1. `test_rotating_file_handler_creation_and_formatting`: Verifies `configure_logging()` creates `logs/app.log` and formats messages across multiple logger hierarchies.
2. `test_log_file_rotation`: Emits messages exceeding `max_bytes` and asserts that backup files (`app.log.1`, `app.log.2`) are generated without exceeding `backup_count`.
3. `test_telemetry_queue_saturation_logged_to_file`: Verifies that queue drops write structured warning entries into the file log.
4. `test_telemetry_db_write_failure_fallback_logging`: Simulates SQLite operational failures and asserts that unwritten batches are saved to `logs/telemetry_fallback.jsonl` without raw parameter dumps polluting `app.log`.
5. `test_feedback_smtp_failure_logging_and_jsonl_persistence`: Simulates SMTP failure and asserts warning in `app.log` while verifying submission persistence in `feedback.jsonl` without password leaks.
6. `test_unhandled_500_exception_logged_with_traceback`: Simulates server crash and asserts complete Python stack trace is captured in `app.log`.
7. `test_rate_limit_warning_logged_to_app_log`: Asserts 429 rate limit triggers emit structured warnings with client IP and route path.
8. `test_normal_200_ok_request_does_not_log_to_app_log`: Verifies that healthy `200 OK` requests remain completely silent in `app.log`.
9. `test_file_upload_invalid_delimiter_logged_to_app_log`: Verifies invalid delimiter uploads log structured warnings with filename.
10. `test_file_upload_corrupt_line_warning_logged_to_app_log`: Verifies invalid coordinate lines log warnings detailing filename, line number, and error reason.
