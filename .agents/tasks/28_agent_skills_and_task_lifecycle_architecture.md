# Task 28: Agent Skills and Task Lifecycle Architecture Upgrade

## Objectives
- [x] Compress `webapp/CHANGELOG.md` by consolidating post-refactor changes under release `## [1.1.0] - 2026-10-18` across `### Added`, `### Changed`, and `### Fixed`, pruning intermediate debugging narratives and trial hex codes, and leaving an empty `## [Unreleased]` section at top.
- [x] Create `.agents/skills/task-init/SKILL.md` implementing universal task initialization, dynamic living doc indexing across `docs/`, `webapp/design_guide.md`, and `.agents/tasks/`, runtime service discovery (ports 8000 and 5173), and a 4-bullet readiness briefing.
- [x] Upgrade `.agents/skills/task-record/SKILL.md` into a unified documentation orchestrator subsuming `doc-record` and `webapp-record` with a 5-step workflow (Git diff inspection, session task recording in `.agents/tasks/`, living doc sync via mapping matrix, changelog/README maintenance, and quality gates).
- [x] Create `.agents/skills/e2e-test/SKILL.md` dedicated to Playwright E2E browser test execution, visual UI runner (`.\test-e2e.ps1 -Ui`), targeted spec commands (`smoke.spec.ts`, `point.spec.ts`, etc.), and diagnostic failure trace inspection.
- [x] Update `.agents/skills/regression-verify/SKILL.md` Step 2 to include `tests/test_feedback.py` and `tests/test_config.py` alongside modern, legacy, rate_limiter, and telemetry test modules.
- [x] Retire obsolete skills by removing `.agents/skills/webapp-init/` and `.agents/skills/webapp-record/`, and replacing `.agents/skills/doc-record/SKILL.md` with an active forwarding redirect to `task-record`.
- [x] Synchronize documentation and rules:
  - `webapp/README.md`: Document `playwright.config.ts`, `e2e/`, and Playwright test commands.
  - `.agents/rules/webapp-design.md`: Update primary button styling to Option #8 Subtle Surface + Border (`#1f2630` / `#eaeef2`).
  - `GEMINI.md`: Update Section 5 workflow list to reference `task-init`, `code-audit`, `regression-verify`, `e2e-test`, `task-record`, `task-commit`.
- [x] Verify complete test suite and code quality:
  - TypeScript compilation and Vite build (`npm --prefix webapp run build` -> 0 errors).
  - Vitest unit test suite (`npm --prefix webapp test` -> 68/68 passed across 5 test files in 312ms).
  - Python Ruff linter & formatter (`ruff check`, `ruff format --check` -> all passed).
  - Core geodetic tests (`pytest tests/test_trans_ro.py tests/test_config.py` -> 62/62 passed).

## Implementation Details

1. **Changelog Compression (`webapp/CHANGELOG.md`)**:
   - Structured all 1.1.0 capabilities into clean architectural summaries:
     - `Added`: Playwright Chromium E2E suite, GitHub Primer dual-theme neutral palette, C2 Graticule brand logo, interactive timeline legend color pickers, status-displaying language dropdown, floating panels stylesheet (`panels.css`), modular CSS architecture, typed UI factory suite (`src/ui/`), Vitest test infrastructure, floating geodetic help panel, and decoupled map layers.
     - `Changed`: 2-line lowercase header attribution, Web GIS & spatial matrix theme reactivity, design system manual overhaul, decoupled exporters, and coordinate order / angular unit sync.
     - `Fixed`: "No data on grid" zone handling, unspaced DMS seconds vs South collision, point parser label regex fix, strict numeric validation, single point map recursion fix, file upload queue reset sync, and popup HTML quote escaping.
   - Cleaned out intermediate trial notes and placed an empty `## [Unreleased]` section at top for future tracking.

2. **Universal Task Initialization (`.agents/skills/task-init/SKILL.md`)**:
   - Implements dynamic context priming: maps user intent to target subsystems and queries the authoritative technical documentation in `docs/` (`architecture.md`, `core_engine.md`, `grid.md`, `rest_api.md`, `web_app.md`, `logging.md`, `testing.md`, `deployment.md`) and `.agents/tasks/`.
   - Proposes non-disruptive daemon orchestration for FastAPI (`127.0.0.1:8000`) and Vite (`localhost:5173`) using `D:\anaconda3\envs\pytransdat\python.exe`.
   - Concludes with a 4-bullet readiness briefing.

3. **Unified Documentation Orchestrator (`.agents/skills/task-record/SKILL.md`)**:
   - Consolidates documentation workflows into a single lifecycle:
     - Step 1: Git diff inspection across layers (`pytransdatro/`, `api/`, `webapp/`, `tests/`, `docs/`, `.agents/`).
     - Step 2: Session Task Recording (`.agents/tasks/NN_<slug>.md`).
     - Step 3: Living documentation synchronization using the subsystem mapping table.
     - Step 4: High-signal changelog and README updates.
     - Step 5: Test and quality gate execution.

4. **Playwright E2E Skill (`.agents/skills/e2e-test/SKILL.md`)**:
   - Documents headless execution, interactive time-travel runner (`.\test-e2e.ps1 -Ui`), HTML report viewing, trace inspection, and targeted execution matching the actual specs (`smoke.spec.ts`, `point.spec.ts`, `batch.spec.ts`, `file.spec.ts`, `map.spec.ts`, `i18n.spec.ts`, `errors.spec.ts`).

5. **Regression Verification Update (`.agents/skills/regression-verify/SKILL.md`)**:
   - Expanded Step 2 command and description to include `tests/test_feedback.py` and `tests/test_config.py` alongside `test_api_modern.py`, `test_api_legacy.py`, `test_api_rate_limiter.py`, and `test_telemetry.py`.

6. **Skill Retirement and Forwarding**:
   - Removed `.agents/skills/webapp-init/` and `.agents/skills/webapp-record/`.
   - Configured `.agents/skills/doc-record/SKILL.md` to redirect directly to `task-record`.

7. **Documentation and Design Rule Synchronization**:
   - `webapp/README.md`: Documented `e2e/`, `playwright.config.ts`, test commands, and 68 Vitest / 19 Playwright tests.
   - `.agents/rules/webapp-design.md`: Updated Section 4 primary button specification to Option #8 (`#1f2630` / `#eaeef2`).
   - `GEMINI.md`: Updated Section 5 workflow to the 6-stage lifecycle (`task-init` -> `code-audit` -> `regression-verify` -> `e2e-test` -> `task-record` -> `task-commit`).

## Verification & Testing
- **Vitest Unit Tests**: `npm --prefix webapp test` -> 68/68 passed (5 files, 312ms).
- **Vite Production Build**: `npm --prefix webapp run build` -> compiled in 260ms, 0 errors.
- **Python Ruff Linter**: `ruff check pytransdatro api tests tools` -> passed (0 errors).
- **Python Ruff Format**: `ruff format --check pytransdatro api tests tools` -> 51 files already formatted.
- **Python Core Geodetic Suite**: `pytest tests/test_trans_ro.py tests/test_config.py` -> 62/62 passed in 3.10s.
