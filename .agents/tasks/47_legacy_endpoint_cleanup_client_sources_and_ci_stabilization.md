# Legacy Endpoint Cleanup, Client Source Registry, and CI Stabilization

## Objectives
- [x] Investigate and eliminate redundant `/cooOpService` alias in favor of the canonical `/transdatonline/cooOpService` legacy endpoint.
- [x] Implement centralized client origin registry (`api/client_sources.py`) and refine `get_client_source()` to map missing headers to `1` (`rest_api`), recognized origins to explicit IDs, and unrecognized values to `0` (`unknown`).
- [x] Diagnose and resolve GitHub Actions CI failures on push to `main` across frontend and backend workflows.
- [x] Fix `Frontend Typecheck & Vitest` CI failure by synchronizing `webapp/package-lock.json` version with `webapp/package.json` for `npm ci`.
- [x] Fix Python CI test failures (`test_webapp_static_mount`, `test_telemetry_effective_date_snapping_and_parity`, and `test_api_legacy_stats.py` suite) by ensuring tests are hermetic and do not depend on uncommitted/gitignored local workstation artifacts.
- [x] Synchronize living documentation in `docs/rest_api.md`, `docs/deployment.md`, and reverse proxy `nginx.conf`.

---

## Implementation Details

1. **Legacy Endpoint Normalization (`api/routes/legacy.py`, `api/main.py`)**:
   - Dropped redundant `@router.get("/cooOpService")` and `@router.post("/cooOpService")` handlers.
   - Retained `/transdatonline/cooOpService` with mandatory `cooOp` and `coos` query parameter validation, preserving byte-for-byte fidelity with the historical ANCPI/Romgeo Java servlet.
   - Updated the custom Starlette `HTTPException` handler in `api/main.py` to target `/transdatonline/cooOpService` exclusively.
   - Updated `nginx.conf` and `docs/deployment.md` to route `/transdatonline/*` without the dead alias.

2. **Centralized Client Source Registry & Header Resolution (`api/client_sources.py`, `api/dependencies.py`)**:
   - Created `api/client_sources.py` housing `CLIENT_SOURCE_LOOKUP = {"rest_api": 1, "web_app": 2, "desktop_cli": 3, "unknown": 0}`.
   - Refactored `get_client_source()` in `api/dependencies.py`:
     - Empty or missing `X-Client-Source`: Resolves to `1` (`rest_api`) for seamless compatibility with scripts, QGIS, and third-party tools.
     - Known origin header: Resolves to designated integer source ID.
     - Unknown or unexpected header: Resolves to `0` (`unknown`), preventing unclassified clients from polluting REST API metrics.
   - Added integration test cases in `tests/test_api_modern.py` verifying unknown and unrecognized client header resolutions.

3. **Frontend Lockfile Synchronization (`webapp/package-lock.json`)**:
   - Synchronized `webapp/package-lock.json` root version to `"version": "1.0.0-beta"` matching `webapp/package.json`.
   - Resolves `npm ci` strict lockfile validation failures on GitHub Actions runner (`ubuntu-latest`).

4. **Hermetic Test Hardening (`tests/test_api_modern.py`, `tests/test_api_legacy_stats.py`)**:
   - `test_webapp_static_mount`: Added fallback test fixture with a stub `index.html` when `webapp/dist/` is not pre-built in Python-only CI environments.
   - `test_telemetry_effective_date_snapping_and_parity`: Isolated test execution by dynamically creating and seeding an in-memory/temporary SQLite database with March 2025 monthly data via `app.dependency_overrides[get_db_path]`.
   - `tests/test_api_legacy_stats.py`:
     - Added `test_legacy_stats_endpoints_when_db_missing` to explicitly verify empty responses when `logs/legacy_stats.db` is absent.
     - Added `synthetic_legacy_db` fixture allowing daily, monthly, and source filtering queries to run and verify SQL aggregation logic in CI.
     - Guarded real-archive assertions with `@pytest.mark.skipif(not HAS_REAL_LEGACY_DB, ...)`.

5. **Living Documentation & Reverse Proxy Sync (`docs/rest_api.md`, `docs/deployment.md`, `nginx.conf`)**:
   - Updated `docs/rest_api.md` to document canonical `/transdatonline/cooOpService`, remove `/cooOpService` alias mentions, and detail `X-Client-Source` header resolution.
   - Updated `nginx.conf` and `docs/deployment.md` reverse proxy location matchers.

---

## Verification & Testing

- **Ruff Code Hygiene**:
  - `ruff check --no-cache .`: Passed (0 errors).
  - `ruff format --check pytransdatro api tests tools`: Passed (56 files clean).
- **Core Geodetic Parity**:
  - `pytest tests/test_trans_ro.py tests/test_config.py`: **62 / 62 passed** (0 numerical regressions, $< 10^{-9}$ tolerance maintained).
- **Web API & Telemetry Integration**:
  - `pytest tests/test_api_*.py tests/test_telemetry.py tests/test_feedback.py tests/test_logging.py`: **142 / 142 passed** (100% pass rate).
- **Frontend TypeScript Build & Vitest**:
  - `npm --prefix webapp run build`: Passed (production bundle minified cleanly).
  - `npm --prefix webapp test`: **83 / 83 passed** across 6 test suites.
