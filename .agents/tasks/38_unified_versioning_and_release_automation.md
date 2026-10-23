# Task 38: Unified Application Versioning & Release Workflow Automation

## Objectives
- [x] Establish a strict Single-Source-of-Truth (SSoT) versioning architecture across the Python core library, FastAPI backend, and Vite/TypeScript web application.
- [x] Create an atomic version synchronization CLI tool (`tools/bump_version.py`) validating SemVer 2.0.0 compliance.
- [x] Dynamically inject the application version into the FastAPI backend (`api/main.py`, `/health` probe) from `pytransdatro.__version__`.
- [x] Dynamically inject the application version into the Vite frontend (`webapp/vite.config.ts`, `webapp/src/config.ts`, `Header.ts`, `FeedbackPanel.ts`) at build time from `package.json`.
- [x] Enhance the `task-commit` agent skill (`.agents/skills/task-commit/SKILL.md`) with release mode detection, pre-flight version consistency checks, conventional release commits (`chore(release): v<version>`), and interactive Git tagging.
- [x] Synchronize living documentation in `docs/architecture.md`.

---

## Implementation Details

1. **Version Synchronization CLI Tool (`tools/bump_version.py`)**:
   - Implemented a pure Python (standard library only) CLI utility to manage project versions.
   - Enforces strict SemVer 2.0.0 validation (`^(\d+)\.(\d+)\.(\d+)(?:-([0-9A-Za-z.-]+))?$`).
   - Supports:
     - `--check`: Verifies that `pytransdatro/__init__.py`, `pyproject.toml`, and `webapp/package.json` are synchronized and valid.
     - `--current`: Displays the active project version.
     - `<new_version>`: Atomically updates all three files, validating syntax before disk writes.
   - Initialized baseline version to `1.0.0-beta`.

2. **Backend API Dynamic Versioning (`api/main.py`, `api/routes/health.py`)**:
   - Replaced hardcoded version strings (`"1.0.0"`) with `pytransdatro.__version__`.
   - `FastAPI(..., version=pytransdatro_version)` now dynamically exposes the package version in OpenAPI `/openapi.json` and Swagger `/docs`.
   - The `/health` probe dynamically reports `"engine_version": pytransdatro_version`.
   - Updated `tests/test_api_modern.py` to assert that `/health` returns the exact engine and app version.

3. **Frontend Dynamic Versioning (`webapp/`)**:
   - In `webapp/vite.config.ts`, injected `__APP_VERSION__` via Vite's `define` by reading `package.json`.
   - In `webapp/src/config.ts`, added `appVersion` to `APP_CONFIG`, populated from `__APP_VERSION__` with fallback to `import.meta.env.VITE_APP_VERSION`.
   - In `webapp/src/components/Header.ts`, updated the `.brand-stage-badge` to render `v${APP_CONFIG.appVersion.replace(/^v/, '')}`.
   - In `webapp/src/components/FeedbackPanel.ts`, dynamically sets `app_version` from `APP_CONFIG.appVersion`.
   - In `webapp/src/i18n/ro.ts` and `en.ts`, updated `version` key to template `v{version}`.

4. **Agent Release Automation (`.agents/skills/task-commit/SKILL.md`)**:
   - Extended the `task-commit` skill to recognize release workflows (*"commit for a new release"*, *"release commit"*, *"commit v1.0.0"*, or `--release`).
   - Step 1 executes `python tools/bump_version.py --check` and verifies that precision and unit tests pass before committing.
   - Step 2 drafts conventional release commit messages: `chore(release): v<version>`.
   - Step 3 adds interactive approval for creating annotated Git tags (`v<version>`).
   - Step 4 executes targeted staging, committing, and tagging (`git tag -a v<version> -m "PyTransDatRO v<version>"`).
   - Step 5 provides exact remote push guidance (`git push origin <branch> --tags`).

5. **Living Documentation Synchronization (`docs/architecture.md`)**:
   - Documented the Single-Source-of-Truth pattern, runtime injection mechanism, and automated release lifecycle.

---

## Verification & Testing

- **Version Synchronization**:
  - `python tools/bump_version.py --check` -> `OK: All versions synchronized and SemVer compliant: 1.0.0-beta`
- **Python Quality & Tests**:
  - `ruff check --no-cache api tools tests` -> PASS (All checks passed)
  - `ruff format --check --no-cache api tools tests` -> PASS (44 files already formatted)
  - Full pytest suite (`pytest tests/`) -> PASS (196 passed, 2 warnings in 19.71s)
  - API Health Check test (`tests/test_api_modern.py::test_health_check`) -> PASS
- **Frontend Quality & Build**:
  - Vitest suite (`npm --prefix webapp test`) -> PASS (6 test files, 80 passed in 439ms)
  - Production bundle build (`npm --prefix webapp run build`) -> PASS (`dist/` built in 304ms without TypeScript errors)
