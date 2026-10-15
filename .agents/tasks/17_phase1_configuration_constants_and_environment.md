# Phase 1: Configuration, Constants & Environment Grounding

## Objectives
- [x] Configure Ruff linter and formatter in `pyproject.toml` with PEP 8 import sorting and 100-character line lengths.
- [x] Harden `.gitignore` to prevent leaking `.env` environment files and local SQLite database files.
- [x] Update Python dependencies in `requirements-api.txt` and create `requirements-dev.txt`.
- [x] Provide template `.env.example` and `webapp/.env.example` documenting all configuration knobs.
- [x] Extract core geodetic and projection invariants into pure Python module `pytransdatro/constants.py` with zero regressions.
- [x] Implement centralized configuration `api/config.py` backed by `pydantic-settings.BaseSettings` with FastAPI dependency injection.
- [x] Create typed frontend configuration `webapp/src/config.ts` replacing hardcoded centroid, zoom, and volume warning thresholds.
- [x] Create comprehensive unit tests in `tests/test_config.py` verifying settings defaults, overrides, and constants integrity.

---

## Implementation Details

1. **Ruff Tooling & Project Grounding (`pyproject.toml`, `requirements-api.txt`, `requirements-dev.txt`, `.gitignore`)**:
   - Established `pyproject.toml` with project metadata, build requirements, and Ruff configuration targeting Python 3.10+, 100-character line length, and automated isort grouping for first-party imports (`pytransdatro`, `api`).
   - Hardened `.gitignore` with `.env`, `.env.*`, `*.db`, `*.sqlite`, `*.sqlite3`, `.ruff_cache/`, and `.pytest_cache/`.
   - Updated `requirements-api.txt` to include `pydantic-settings>=2.2.0` and added `requirements-dev.txt` for development dependencies (`ruff`, `pytest`, `httpx`).
   - Documented operational settings in `.env.example` (prefixed with `PYTRANSDAT_`) and `webapp/.env.example` (`VITE_*`).

2. **Core Geodetic Constants Extraction (`pytransdatro/constants.py`, `proj_stereo.py`, `utils.py`)**:
   - Extracted physical constants without altering mathematical precision:
     - WGS84 Ellipsoid: `WGS84_ELL_SMAJ_AXIS = 6378137.0`, `WGS84_ELL_INV_FLAT = 298.257223563`.
     - Stereo 70 Projection: Precomputed natural origin in radians (`46.0 * DEG_TO_RAD_FACTOR` ~0.8028514559173916, `25.0 * DEG_TO_RAD_FACTOR` ~0.4363323129985824, equivalent to 46° 00' 00.0" N and 25° 00' 00.0" E), eliminating runtime `sexa_to_rad` conversion upon instantiation; false coordinates (`500000.0`, `500000.0`), scale factor (`0.99975`), convergence tolerance (`0.0000000000484814`), and max iterations (`50`).
     - Angular conversion factors: `DEG_TO_RAD_FACTOR = math.pi / 180.0`, `RAD_TO_DEG_FACTOR = 180.0 / math.pi`.
     - Romania geographic bounding box: `ROMANIA_BBOX_LAT_MIN` (43.5), `ROMANIA_BBOX_LAT_MAX` (48.3), `ROMANIA_BBOX_LON_MIN` (20.2), `ROMANIA_BBOX_LON_MAX` (29.8).
   - Preserved backward-compatible class-level attributes on `StereoProj` and re-exports on `utils`.

3. **Centralized Backend Configuration (`api/config.py`, `dependencies.py`, `main.py`, `file.py`, `batch.py`)**:
   - Created `Settings` class using `pydantic-settings.BaseSettings` reading `PYTRANSDAT_*` environment variables with sensible defaults for local development.
   - Provided dependency injection via `get_settings` in `api/dependencies.py`.
   - Integrated `Settings` in `api/main.py` for CORS origins, telemetry logger directory, and rate limiting parameters.
   - Updated streaming file transformation (`api/routes/file.py`) and batch transformation (`api/routes/batch.py`) to enforce configurable limits (`max_file_points`, `max_batch_points`).

4. **Frontend Environment Configuration (`webapp/src/config.ts`, `MapView.ts`, `apiClient.ts`)**:
   - Created `APP_CONFIG` in `webapp/src/config.ts` reading from `import.meta.env.VITE_*` with defaults for `defaultMapCenter: [45.9432, 24.9668]`, `defaultMapZoom: 7`, `batchWarningThreshold: 2000`, `maxBatchPoints: 20000`, and `maxFileSizeMb: 50`.
   - Wired `APP_CONFIG` into `MapView.ts` for map initialization and point volume warning triggers.
   - Wired `APP_CONFIG.apiBaseUrl` into `apiClient.ts` through an `apiUrl()` URL-resolving helper function.

---

## Verification & Testing

- **Configuration Unit Tests (`tests/test_config.py`)**:
  - `5/5 passed` (0.36s): Default values, environment variable overrides, CORS string parsing, physical and projection constant values, backward compatibility re-exports.
- **Core Geodetic Unit Suite (`tests/test_trans_ro.py`)**:
  - `46/46 passed` (2.58s): Single/batch 2D and 3D transformations, roundtrips, out-of-grid handling, and dynamic grid interpolation binding.
- **API Functional & Rate Limiting Tests (`tests/test_api_legacy.py`, `tests/test_api_rate_limiter.py`, `tests/test_api_modern.py`)**:
  - `33/33 passed` (2.96s): 100% backward-compatible emulation of legacy service, rate-limiter enforcement, and modern JSON streaming endpoints.
- **Telemetry & Feedback Tests (`tests/test_telemetry.py`, `tests/test_feedback.py`)**:
  - `9/9 passed` (0.81s): Asynchronous SQLite WAL logging, aggregations, and feedback submission.
- **Code Audit & Ruff Lint**:
  - `ruff check`: All checks passed with 0 errors across Phase 1 files.
  - `ruff format`: Formatted 7 files to 100-character line length and PEP 8 import sorting.
  - Layer boundary check: Confirmed `pytransdatro/` has zero dependencies or imports from `api/` or web frameworks.
- **Frontend Build**:
  - `npm --prefix webapp run build`: 36 modules transformed, zero TypeScript errors, production build in 459ms.
