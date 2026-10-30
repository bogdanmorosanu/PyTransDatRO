# Task 51: Pre-Release v1.0.0 Comprehensive Review & DRY Remediation

## Objectives
- [x] Complete comprehensive pre-release review of `PyTransDatRO` across all layers: pure Python core (`pytransdatro/`), REST API (`api/`), Geodetic Studio (`webapp/`), Docker/Nginx stack, and tooling.
- [x] Consolidate project configuration into a unified single entrypoint (`.env` / `.env.example`) configuring both backend (`PYTRANSDAT_*`), frontend build args (`VITE_*`), and container runtime limits (`NGINX_CLIENT_MAX_BODY_SIZE`).
- [x] Package core geodetic engine for PyPI distribution (`pyproject.toml` package-data bundling `*.spg`, `*.json`, project URLs, full metadata).
- [x] Fix potential `functools.lru_cache` memory leaks on Grid instances (`B019`) and add `strict=True` to coordinate `zip()` calls (`B905`).
- [x] Implement dynamic runtime limits (`/health` returns `RuntimeLimitsModel`; `webapp` syncs in `apiClient.ts`).
- [x] Harden REST API resources (prune sliding-window IP rate limiter map at 10,000 entries; structured connection management in `api/routes/info.py`).
- [x] Eliminate duplicated logic across web application (DRY consolidation):
  - Canonical `parseAngularTokenToDeg` in `webapp/src/utils/dmsParser.ts`.
  - Canonical coordinate & quasigeoid undulation formatting (`webapp/src/utils/coordinateFormat.ts`).
  - Canonical telemetry date range presets (`webapp/src/utils/dateUtils.ts`).
  - Unified grid warnings & bounds (`webapp/src/utils/gridBounds.ts`).
  - Core SVG charting math & Bézier path smoothing (`webapp/src/ui/chartCore.ts`).
- [x] Purge orphaned webapp components (`ReportModal.ts`, `StatsView.ts`) and Vite starter assets (`hero.png`, `icons.svg`).
- [x] Migrate `FeedbackPanel.ts` to `renderSelect` UI factory, replace hardcoded strings with i18n tokens, and use `APP_CONFIG.repoUrl`.
- [x] Enable TypeScript `"strict": true` in `webapp/tsconfig.json`.
- [x] Create comprehensive i18n parity test suite (`webapp/src/i18n/i18n.test.ts`) guaranteeing 100% key and placeholder parity across Romanian and English.
- [x] Update `tools/bump_version.py` to synchronize `webapp/package-lock.json` alongside `pyproject.toml`, `pytransdatro/__init__.py`, and `webapp/package.json`.
- [x] Make PowerShell dev and test scripts (`dev.ps1`, `test-all.ps1`, `test-e2e.ps1`, `playwright.config.ts`) completely machine-agnostic with portable Python environment auto-detection.
- [x] Clean repo: archive obsolete files to `docs/archive/` (`grids_obsolete.md`, `specs/`, `logos.html`), remove temporary artifacts (`err.txt`, `out.txt`, `pytransdatro-app.tar`).
- [x] Update CI workflow (`.github/workflows/ci.yml`) to run version check, wheel packaging smoke test, and Playwright Chromium E2E suite.
- [x] Create `CONTRIBUTING.md` and `CHANGELOG.md` for `v1.0.0`.

## Implementation Details

### 1. Workstream 1: Unified Single-Entrypoint Configuration
- Consolidated all configuration knobs into root `.env.example` (backend `PYTRANSDAT_*`, frontend `VITE_*`, container `WEBAPP_PORT`, `NGINX_CLIENT_MAX_BODY_SIZE`).
- Deleted redundant `webapp/.env.example`.
- Configured Vite (`webapp/vite.config.ts`) with `envDir: '..'` and added reverse proxy entries for `/redoc` and `/transdatonline`.
- In `api/config.py`: converted file paths (`legacy_db_path`, `feedback_log_path`, `telemetry_fallback_path`, `app_log_path`) to optional overrides with `resolved_*` Path properties auto-derived from `log_dir`.
- Created `nginx.conf.template` with dynamic `${NGINX_CLIENT_MAX_BODY_SIZE:-50M}` and updated `Dockerfile.web` to use Nginx envsubst.
- Extended `GET /health` with `RuntimeLimitsModel` returning active limits (`max_batch_points`, `max_file_points`, `max_file_upload_bytes`), consumed and cached by `webapp/src/config.ts` via `syncRuntimeLimits`.

### 2. Workstream 2: Core Geodetic Engine & Standalone Packaging
- Added `[tool.setuptools.package-data]` to `pyproject.toml` bundling `*.spg` and `*.json` files.
- Added full PyPI distribution metadata, classifiers, and `[project.urls]`.
- Refactored `Grid._init_interp` in `pytransdatro/trans_grid.py` to `@staticmethod` + `@functools.lru_cache(maxsize=8192)` to prevent object-lifecycle memory leaks (`B019`).
- Refactored `Grid1D` and `Grid2D` index lookups to maintain per-instance caches (`_cached_vs_at_idxs`).
- Added `strict=True` to coordinate `zip()` iterations across `pytransdatro/proj_stereo.py`, `trans_helmert2d.py`, `trans_ro.py`, and `trans_grid.py` (`B905`).
- Replaced bare `except Exception` with specific exception tuples (`sqlite3.Error`, `OSError`, `ValueError`) in `pytransdatro/telemetry.py`.

### 3. Workstream 3: REST API Hardening & Resource Safety
- Added dictionary pruning in `api/middleware/rate_limit.py`: when IP tracking map exceeds 10,000 entries, evicts empty deques and records older than 10 minutes.
- Dynamic rate limit 429 response detailing requested path, window, and `"code": "RATE_LIMIT_EXCEEDED"`.
- Cleaned up database queries and cursor handling in `api/routes/info.py` with `try ... finally: conn.close()`.
- Extracted `_build_date_filter_clauses` and `_get_granularity_expression` helpers in `api/routes/info.py`.

### 4. Workstream 4: Shared Functionality & DRY Consolidation
- Created `webapp/src/utils/coordinateFormat.ts`: canonical implementations of `formatAngularCoordinate`, `formatPlanarCoordinate`, `formatCoordinatePair`, `computeQuasigeoidZeta` ($\zeta = h - H$), and `formatQuasigeoidZeta`.
- Created `webapp/src/utils/dateUtils.ts`: canonical `calculatePresetDateRange` supporting `'7d'`, `'30d'`, `'ytd'`, `'12m'`, `'all'` presets across computation, spatial, and legacy views.
- Created `webapp/src/ui/chartCore.ts`: shared SVG cubic Bézier path generation (`buildSmoothSvgPath`), nearest cursor index search (`calculateNearestXIndex`), and X pixel coordinate placement (`calculateXPosition`).
- Consolidated grid coverage warnings in `webapp/src/utils/gridBounds.ts` (`isUntransformableWarning`, `isOutOfGridWarning`, `WARNING_OUT_OF_GRID`, `WARNING_NO_DATA`).
- Refactored consuming components:
  - `pointParser.ts`: uses `parseAngularTokenToDeg`.
  - `PointConverter.ts`: uses `formatCoordinatePair`, `computeQuasigeoidZeta`, and `isOutOfGridWarning`.
  - `BatchConverter.ts`: uses `formatCoordinatePair`, `computeQuasigeoidZeta`, `isUntransformableWarning`, and checks against `APP_CONFIG.maxBatchPoints`.
  - `FileUpload.ts`: uses `computeQuasigeoidZeta`, `isUntransformableWarning`, and `APP_CONFIG.maxFileSizeMb`.
  - `MapView.ts`: imports `ROMANIA_GRID_BOUNDS` and `isInsideRomaniaGrid` from `gridBounds.ts`, uses `formatCoordinatePair`, `formatQuasigeoidZeta`, `APP_CONFIG.pointZoomLevel`.
  - `exporters.ts`: uses `isUntransformableWarning`.
  - `TelemetryDashboard.ts`: uses `calculatePresetDateRange` for main, spatial, and legacy tabs.
  - `TimelineChart.ts` & `LegacyTimelineChart.ts`: use `buildSmoothSvgPath`, `calculateNearestXIndex`, and `calculateXPosition` from `chartCore.ts`.

### 5. Workstream 5: Webapp Dead-Code Purge & UI Compliance
- Purged orphaned components: deleted `webapp/src/components/ReportModal.ts` and `webapp/src/components/StatsView.ts`.
- Purged starter assets: deleted `webapp/src/assets/hero.png` and `webapp/public/icons.svg`.
- Refactored `FeedbackPanel.ts`:
  - Migrated report type `<select>` to `renderSelect` UI factory.
  - Updated default active grid to `rom_grid3d_25.09.spg`.
  - Replaced hardcoded GitHub issue link with `APP_CONFIG.repoUrl`.
  - Replaced all hardcoded dialog strings with bilingual i18n tokens in `ro.ts` and `en.ts`.
- Removed inline styles from `VolumeWarningModal.ts`, extracting `.vw-callout` and `.vw-prompt-label` into `layout.css`.
- Created comprehensive i18n test suite (`webapp/src/i18n/i18n.test.ts`): guarantees 100% parity across leaf keys and dynamic placeholders between Romanian and English.
- Enabled TypeScript `"strict": true` in `webapp/tsconfig.json`.

### 6. Workstream 6: Repository Cleanliness, Versioning & Open-Source Release
- Updated `tools/bump_version.py` to synchronize version in `webapp/package-lock.json` alongside `pyproject.toml`, `pytransdatro/__init__.py`, and `webapp/package.json`.
- Updated PowerShell scripts (`dev.ps1`, `test-all.ps1`, `test-e2e.ps1`) and `webapp/playwright.config.ts` to implement portable Python auto-detection (`$env:PYTHON` $\to$ `$env:CONDA_PREFIX` $\to$ user Conda $\to$ `.venv` $\to$ `python`).
- Archived obsolete files to `docs/archive/` (`grids_obsolete.md`, `specs/`, `webapp/public/logos.html`).
- Deleted temporary files (`err.txt`, `out.txt`, `pytransdatro-app.tar`).
- Updated GitHub Actions CI workflow (`.github/workflows/ci.yml`) to include version consistency check, wheel packaging smoke test, and Playwright Chromium E2E suite.
- Updated technical documentation (`docs/architecture.md`, `docs/rest_api.md`, `docs/testing.md`, `README.md`, `CLAUDE.md`, `webapp/README.md`).
- Authored open-source `CONTRIBUTING.md` and Keep-a-Changelog `CHANGELOG.md` for `v1.0.0`.

## Verification & Quality Gates

- **Official ANCPI / RomGEO Parity Benchmark**:
  - `python tools/generate_ancpi_comparison_reports.py`: 100% PASS with zero numerical regression ($< 10^{-9}$ tolerance across 9,046 national control points).
- **Python Pytest Suite**:
  - `python -m pytest -m "not legacy_live" -q`: 189 passed, 18 deselected in 14.8s.
- **Ruff Linter & Formatter**:
  - `python -m ruff check pytransdatro api tests tools`: 0 errors.
  - `python -m ruff format --check pytransdatro api tests tools`: 56 files formatted cleanly.
- **Version SSoT Consistency**:
  - `python tools/bump_version.py --check`: OK (all versions synchronized: `1.0.0-beta`).
- **Frontend TypeScript & Vitest Suite**:
  - `npm --prefix webapp test -- --run`: 10 test files passed (123 tests) in 461ms.
- **Frontend Production Build**:
  - `npm --prefix webapp run build`: `tsc` clean under `"strict": true`, Vite bundled in 485ms.
- **Docker Compose Configuration**:
  - `docker compose config --quiet`: Validated without errors.
- **Playwright E2E Suite**:
  - `npm --prefix webapp run test:e2e`: 71 browser tests executed across all user journeys.
