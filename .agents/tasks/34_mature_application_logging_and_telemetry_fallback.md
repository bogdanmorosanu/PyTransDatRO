# Mature Application Logging Architecture & Telemetry Emergency Fallback

## Objectives
- [x] Implement symmetrical separation between operational health logging (`logs/app.log`) and data storage payloads (`logs/feedback.jsonl`, `logs/telemetry_fallback.jsonl`).
- [x] Enforce the **"Errors & Warnings Only"** logging discipline across the API: normal `200 OK` operations remain completely silent with zero performance penalty.
- [x] Route unwritten telemetry batches to `logs/telemetry_fallback.jsonl` upon persistent SQLite write failure without dumping raw parameter arrays into `logs/app.log`.
- [x] Add actionable WARNING logging in `api/services/mailer.py` when SMTP feedback email delivery fails, while safely preserving user feedback in `logs/feedback.jsonl` with zero PII/password leakage.
- [x] Implement global unhandled 500 error exception logging with full Python stack tracebacks in `api/main.py`.
- [x] Add structured WARNING logs in `api/middleware/rate_limit.py` when clients exceed request rate limits (HTTP 429).
- [x] Add diagnostic WARNING logs in `api/routes/file.py` when streaming CSV/TXT uploads fail due to invalid delimiters, oversized payloads, or corrupt coordinate lines (with line number and error reason).
- [x] Add configuration setting `telemetry_fallback_path` in `api/config.py` and `.env.example`.
- [x] Expand automated test coverage in `tests/test_logging.py` to 10 comprehensive scenarios validating all logging behaviors.
- [x] Update living documentation in `docs/log.md` and `docs/telemetry.md`.

## Implementation Details

1. **Telemetry Emergency Fallback (`pytransdatro/telemetry.py`)**:
   - Added `_save_fallback_payload()` method appending unwritten aggregated computation and spatial stats as structured JSON lines to `self.fallback_path` (`logs/telemetry_fallback.jsonl`).
   - Replaced raw parameter dumps in `logger.error` with clean 1-line messages indicating failure count and fallback destination path.
   - Updated `get_diagnostics()` to expose `fallback_path`.

2. **Feedback & Mailer Error Logging (`api.services.mailer`, `api/config.py`)**:
   - In `dispatch_feedback()`, emit structured `WARNING` to `logs/app.log` when SMTP fails (`SMTPConnectError`, timeouts, etc.), noting that feedback was safely preserved in `logs/feedback.jsonl`.
   - Updated default `feedback_log_path` to `./logs/feedback.jsonl` and ensured fallback resolution supports both direct custom settings and relative `log_dir` configurations.

3. **Core API Unhandled Exception & Rate Limiting (`api/main.py`, `api/middleware/rate_limit.py`)**:
   - Added global exception handler in `create_app()` logging `logger.exception("Unhandled server exception on %s %s: %s", request.method, request.url.path, exc)` for 500 responses.
   - Added `Settings` parameter injection to `create_app(settings=...)` enabling isolated test environments.
   - Added structured `logger.warning("Rate limit exceeded for client %s on %s %s ...")` in `SlidingWindowRateLimiter`.

4. **Streaming File Upload Diagnostics (`api/routes/file.py`)**:
   - Emits structured warnings when file uploads exceed size or line limits, when delimiter auto-detection fails, or when a coordinate line contains corrupt data (logging filename, row number, and reason with a cap of 5 warnings per upload to prevent log flood).

5. **Documentation & Configuration Synchronization (`docs/log.md`, `docs/telemetry.md`, `.env.example`)**:
   - Documented the three-file architecture (`app.log`, `feedback.jsonl`, `telemetry_fallback.jsonl`).
   - Detailed the "Errors & Warnings Only" discipline, log rotation parameters (10 MB cap, 5 backups, 60 MB maximum disk footprint), and sample log entries.

## Verification & Automated Testing

- **Logging Test Suite (`tests/test_logging.py`)**:
  - 10/10 passed:
    - `test_rotating_file_handler_creation_and_formatting`: PASSED
    - `test_log_file_rotation`: PASSED
    - `test_telemetry_queue_saturation_logged_to_file`: PASSED
    - `test_telemetry_db_write_failure_fallback_logging`: PASSED
    - `test_feedback_smtp_failure_logging_and_jsonl_persistence`: PASSED
    - `test_unhandled_500_exception_logged_with_traceback`: PASSED
    - `test_rate_limit_warning_logged_to_app_log`: PASSED
    - `test_normal_200_ok_request_does_not_log_to_app_log`: PASSED
    - `test_file_upload_invalid_delimiter_logged_to_app_log`: PASSED
    - `test_file_upload_corrupt_line_warning_logged_to_app_log`: PASSED
- **Telemetry Test Suite (`tests/test_telemetry.py`)**:
  - 28/28 passed.
- **Feedback Test Suite (`tests/test_feedback.py`)**:
  - 7/7 passed.
- **Configuration Test Suite (`tests/test_config.py`)**:
  - 8/8 passed.
- **Core Geodetic Suite (`tests/test_trans_ro.py`)**:
  - 54/54 passed with zero numerical regression against ANCPI benchmarks.
- **Full Python Test Suite (`pytest`)**:
  - 195 passed, 0 failed in 14.89s.
- **Frontend Vitest Suite (`npm test`)**:
  - 80/80 passed in 373ms.
- **Code Quality & Typing Audit (`ruff check`, `ruff format`)**:
  - 0 errors across 42 files.
