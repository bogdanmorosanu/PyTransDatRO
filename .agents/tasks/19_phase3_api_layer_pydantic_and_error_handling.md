# Phase 3: API Layer, Pydantic & Error Handling

## Objectives
- [x] Configure standard library logging via `logging.dictConfig` driven by centralized settings in `api/config.py`.
- [x] Replace raw `print()` statements with structured loggers across routers and services.
- [x] Audit and eliminate all silent failures (`except Exception: pass`) across `api/routes/info.py`, `api/routes/file.py`, `api/routes/legacy.py`, and `api/services/mailer.py`.
- [x] Introduce uniform `ErrorResponse` Pydantic model (`api/schemas/error.py`) and global exception handlers in `api/main.py` while preserving 100% backward compatibility for existing `detail` consumers.
- [x] Decouple hardcoded `./logs` filesystem paths via FastAPI dependency injection (`get_db_path`).
- [x] Verify zero regression across all API endpoints, rate limiting, telemetry, and official ANCPI geodetic benchmarks.

## Implementation Details

1. **Standardized Logging Configuration (`api/config.py`, `api/main.py`)**:
   - Implemented `get_logging_config(settings: Optional[Settings] = None) -> dict` returning a standard `dictConfig` dictionary.
   - Implemented `configure_logging(settings: Optional[Settings] = None) -> None` called during app factory creation and worker lifespan boot.
   - Added centralized SMTP parameters to `Settings` (`smtp_host`, `smtp_port`, `smtp_user`, `smtp_password`, `smtp_use_tls`, `feedback_recipient_email`, `feedback_sender_email`).
   - Standardized feedback JSONL log path configuration.

2. **Mailer Service Refactor (`api/services/mailer.py`)**:
   - Replaced unmanaged `os.environ` lookups with `Settings` dependency.
   - Replaced raw `print()` outputs with module-level `logger = logging.getLogger("api.services.mailer")`.
   - Trapped specific exceptions (`(smtplib.SMTPException, OSError, TimeoutError)` and `(OSError, TypeError, ValueError)`) with structured warning logs instead of generic catches.

3. **Silent Error Elimination (`api/routes/info.py`, `api/routes/legacy.py`, `api/routes/file.py`)**:
   - `api/routes/info.py`: Replaced bare `except Exception: pass` in `_load_comparison_summary()` and `_load_or_extract_metadata()` with `(json.JSONDecodeError, OSError)` trapping and warning logs.
   - `api/routes/info.py`: Replaced broad catch-alls in `get_telemetry_stats()`, `get_telemetry_timeline()`, and `get_telemetry_spatial_matrix()` with explicit `sqlite3.Error` trapping and structured error logging.
   - `api/routes/legacy.py`: Replaced generic `except Exception:` on JSON body and array parsing with explicit `json.JSONDecodeError`.
   - `api/routes/file.py`: Replaced generic `except Exception:` on coordinate token extraction with explicit `(ValueError, TypeError, IndexError)`.
   - `api/config.py`: Narrowed CORS origins JSON parse exception to `(json.JSONDecodeError, TypeError)`.

4. **Pydantic Schemas & Uniform Error Response (`api/schemas/error.py`, `api/main.py`)**:
   - Created `ErrorResponse` schema with `detail: str` and `code: Optional[str] = None`.
   - Exported `ErrorResponse` in `api/schemas/__init__.py`.
   - Registered custom exception handlers for `StarletteHTTPException` and `RequestValidationError` in `api/main.py`, guaranteeing consistent responses and `"VALIDATION_ERROR"` error codes.

5. **FastAPI Dependency Injection (`api/dependencies.py`, `api/routes/`)**:
   - Added `get_db_path(settings: Settings = Depends(get_settings)) -> Path` in `api/dependencies.py`.
   - Injected `get_db_path` into `get_telemetry_stats`, `get_telemetry_timeline`, and `get_telemetry_spatial_matrix`.
   - Injected `get_trans_engine` and `get_client_source` into legacy router endpoints.
   - Injected `get_settings` into `api/routes/feedback.py`.

## Verification & Testing

- **Static Analysis & Linting**:
  - `ruff check api/`: 0 errors
  - `ruff format --check api/`: 20 files verified clean
- **Automated Test Results**:
  - `tests/test_api_legacy.py`: 11/11 passed
  - `tests/test_api_modern.py`: 22/22 passed (+2 new tests for `ErrorResponse` and validation errors)
  - `tests/test_api_rate_limiter.py`: 2/2 passed
  - `tests/test_config.py`: 7/7 passed (+2 new tests for `get_logging_config` and SMTP defaults)
  - `tests/test_feedback.py`: 4/4 passed (+1 new test for custom `Settings` injection)
  - `tests/test_telemetry.py`: 6/6 passed
  - `tests/test_trans_ro.py`: 46/46 passed
  - **Overall Pass Rate**: **98/98 tests passed (100%)**
- **Geodetic Benchmarks**:
  - `tools/generate_ancpi_comparison_reports.py`: 100% byte-for-byte numerical precision across 9,046 national control points ($< 10^{-9}$ tolerance).
- **Web App Compilation**:
  - `npm --prefix webapp run build`: Clean build in 482ms.
