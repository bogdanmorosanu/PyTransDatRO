# CI Pipeline Stabilization & Packaging Fixes

## Objectives
- [x] Investigate failing GitHub Actions workflow jobs (Run `#41158920882`) across Python CI, Frontend CI, and Playwright E2E suites.
- [x] Resolve `ModuleNotFoundError: No module named build` in Python CI packaging smoke test by adding `build>=1.0.0` to `requirements-dev.txt` and `pyproject.toml`.
- [x] Resolve `ReferenceError: navigator is not defined` in Frontend Vitest suite by defensively guarding `navigator.clipboard` access and initializing mock navigator in `webapp/src/ui/ui.test.ts`.
- [x] Modernize `.github/workflows/ci.yml` by updating runner Node.js from `20` to `22` (Active LTS), filtering live external service calls (`-m "not legacy_live"`), adding `--no-isolation` to wheel packaging, and executing a module smoke import check.
- [x] Verify complete quality gates across Python ruff/pytest suites and frontend bundle/Vitest suites.

## Implementation Details

1. **Python Packaging Dependencies**:
   - Added `build>=1.0.0` to [`requirements-dev.txt`](../../requirements-dev.txt) and `pyproject.toml` (`[project.optional-dependencies] dev`).
   - Ensures `python -m build --wheel` can execute in headless continuous integration runners without missing dependencies.

2. **Frontend DOM Globals & Test Mock Hardening**:
   - In [`webapp/src/ui/copyButton.ts`](../../webapp/src/ui/copyButton.ts) and [`webapp/src/ui/copyBadge.ts`](../../webapp/src/ui/copyBadge.ts), replaced direct `navigator.clipboard` property accesses with `typeof navigator !== 'undefined' && navigator.clipboard?.writeText`.
   - In [`webapp/src/ui/ui.test.ts`](../../webapp/src/ui/ui.test.ts), safeguarded mock setup so that `(globalThis as any).navigator` is explicitly initialized if undefined before defining the `clipboard.writeText` spy, preventing silent failures on Node < 21.

3. **CI Workflow Modernization**:
   - In [`.github/workflows/ci.yml`](../../.github/workflows/ci.yml):
     - Added `-m "not legacy_live"` to `Web API Integration & Telemetry Tests` to prevent flaky failures caused by external network latency or downtime of legacy remote endpoints.
     - Enhanced `Core Library Wheel Packaging Smoke Test` with `--no-isolation` and added an explicit Python import smoke test (`python -c "import pytransdatro; print('PyTransDatRO wheel smoke test OK:', pytransdatro.__version__)"`).
     - Bumped `node-version` from `20` to `22` in `frontend-ci` and `e2e-ci` to resolve GitHub Actions Node 20 deprecation warnings.

## Verification & Quality Gates
- **Version Consistency**: Checked via `python tools/bump_version.py --check` (synchronized at 1.0.0-beta).
- **Ruff Lint & Format**: Passed cleanly across `pytransdatro`, `api`, `tests`, and `tools` (0 errors, 57 files formatted).
- **Core Geodetic Parity**: `pytest tests/test_trans_ro.py tests/test_config.py` passed (63/63 tests).
- **Web API & Telemetry**: `pytest -m "not legacy_live"` passed (126/126 selected tests, 18 deselected live tests).
- **Wheel Packaging & Smoke Import**: `python -m build --wheel --no-isolation` built wheel successfully; verified clean import of `pytransdatro`.
- **Frontend Production Build**: `npm --prefix webapp run build` succeeded (`dist/` generated).
- **Frontend Vitest Suite**: `npm --prefix webapp test` passed (10/10 test files, 123/123 tests).
- **Playwright E2E TypeScript Check**: `tsc -p e2e` passed (0 type errors).
