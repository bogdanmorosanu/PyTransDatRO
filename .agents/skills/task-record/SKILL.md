---
name: task-record
description: Unified documentation orchestrator that records session outcomes in .agents/tasks/, synchronizes living documentation in docs/ (including docs/webapp/), maintains webapp/README.md, and verifies quality gates. Subsumes doc-record and webapp-record. Use whenever recording, documenting, or updating task results.
---

# Unified Task Record Skill (`task-record`)

## Purpose
This skill serves as the **unified documentation orchestrator** for the PyTransDatRO repository, subsuming the legacy `doc-record` and `webapp-record` skills into a single cohesive post-task workflow. It inspects all touched codebase layers, updates chronological session task logs, synchronizes authoritative living technical documents in `docs/` (including `docs/webapp/`), updates the web application single-source-of-truth README, and verifies quality gates.

---

## Operating Principles

1. **Single Unified Entrypoint**:
   - Replaces disparate documentation steps with a single structured process executed at the conclusion of every development session or significant milestone.
2. **Outcome-Oriented Technical Truth**:
   - Document the *final, verified state* of the work achieved.
   - Omit conversational back-and-forth, intermediate debugging dead-ends, or iterative exploratory trials.
3. **Living Documentation Fidelity (No Documentation Rot)**:
   - Authoritative documents in `docs/` define the current architectural reality, algorithms, schemas, and public contracts. They must never fall out of sync with code changes.
4. **Living Documentation & README Maintenance**:
   - `webapp/changelog.md` is retired; all features and specifications are maintained directly in the living documentation suite `docs/webapp/` (`README.md`, `design_system.md`, `archetypes.md`, `architecture.md`, `actionable_inventory.md`), while session histories are recorded in `.agents/tasks/`.
   - `webapp/README.md` must accurately reflect current component structures, directory trees, and features.
5. **Quality Gate Verification**:
   - Always verify that tests pass and bundles compile cleanly before recording or proposing a git commit.

---

## Documentation Mapping Matrix (`docs/`)

When inspecting modified files in Step 1, identify which authoritative documents in `docs/` must be synchronized in Step 3:

| Subsystem / Touched Layer | Authoritative Document | Contents to Maintain |
| :--- | :--- | :--- |
| **System Architecture** (`pytransdatro/`, `api/`, `webapp/`) | [`docs/architecture.md`](../../../docs/architecture.md) | High-level layer topology, class hierarchy, pipeline flow, grid discovery. |
| **Core Geodetic Engine** (`pytransdatro/`) | [`docs/core_engine.md`](../../../docs/core_engine.md) | Stereographic formulas, 2D Helmert, bicubic spline interpolation, quasigeoid undulation, dual-CRS bounding boxes, 4-stage pipeline, domain exceptions. |
| **Grid Management** (`pytransdatro/grids/`, `grids/`) | [`docs/grid.md`](../../../docs/grid.md) | `.spg` binary format specification, header offsets, planar shift subgrids. |
| **REST API Layer** (`api/`) | [`docs/rest_api.md`](../../../docs/rest_api.md) | Modern endpoints (`/api/v1/*`), Pydantic models, query parameters, error schemas, legacy `/transdatonline` emulation. |
| **Web Application & UI** (`webapp/src/`, `webapp/`) | [`docs/webapp/README.md`](../../../docs/webapp/README.md) | Web GIS Leaflet architecture, workbench tabs, UI factory helpers, modular CSS tokens, component archetypes, living actionable inventory. |
| **Usage Telemetry** (`pytransdatro/telemetry/`, `api/routes/`) | [`docs/telemetry.md`](../../../docs/telemetry.md) | SQLite asynchronous WAL usage logger schema, dynamic window aggregation, coordinate binning. |
| **Application Logging** (`api/config.py`, `logs/app.log`) | [`docs/log.md`](../../../docs/log.md) | Python `logging` dictConfig setup, rotating file handler, failure fallback dumps. |
| **Testing & Quality** (`tests/`, `tools/`, `webapp/e2e/`) | [`docs/testing.md`](../../../docs/testing.md) | Test matrix, numerical precision benchmarks, Vitest unit tests, Playwright E2E browser suite. |
| **Deployment & Ops** (`Dockerfile`, `docker-compose.yml`, `.github/`) | [`docs/deployment.md`](../../../docs/deployment.md) | Docker containerization, Uvicorn settings, reverse proxying, environment configurations. |

---

## Workflow Steps

### Step 1: Git Diff & Layer Inspection
Execute git inspection commands to detect all modified, added, or deleted files:
```powershell
git status -s
git diff --stat
```
Identify all touched layers:
- `webapp/`: Web frontend components, styles, utilities, tests, or E2E specs.
- `api/`: FastAPI routes, schemas, services, or configuration.
- `pytransdatro/`: Core geodetic library, algorithms, or grid loaders.
- `tests/` / `tools/`: Python tests, benchmark scripts, or test fixtures.
- `docs/`: Living documentation files.
- `.agents/`: Agent rules, tasks, or skills.

---

### Step 2: Session Task Recording (`.agents/tasks/NN_<slug>.md`)

1. **Determine Target Path**:
   - Inspect `.agents/tasks/` to determine the highest existing two-digit numeric prefix (e.g., `27` -> next is `28`).
   - If a task file was already created for the current session, update it in-place.
   - If not, format the file path as `.agents/tasks/<NN>_<descriptive_slug>.md`.

2. **Compose Structured Task Content**:
   ```markdown
   # <Clear, Descriptive Task Title>

   ## Objectives
   - [x] <Objective 1 achieved in this session>
   - [x] <Objective 2 achieved in this session>
   - [ ] <Optional: explicitly deferred or planned follow-up if discussed>

   ## Implementation Details

   1. **<Component / Focus Area 1> (`<path/to/file>`)**:
      - Detailed summary of changes, classes, or functions introduced or modified.
      - Design rationale and architectural invariants respected.

   2. **<Component / Focus Area 2> (`<path/to/file>`)**:
      - Details of data structures, algorithms, or UI patterns.

   ## Verification & Testing
   - **Test Modules**: `<tests/test_foo.py>`, `<webapp/src/ui/ui.test.ts>`, `<webapp/e2e/foo.spec.ts>`
   - **Scenarios Covered**: Edge cases, boundary tests, or browser workflows verified.
   - **Results**: Final test pass rate and benchmark residuals.
   ```

3. Write or update the file using `write_to_file` or `replace_file_content`.

---

### Step 3: Living Documentation Synchronization (`docs/`)

Using the **Documentation Mapping Matrix**, determine if any architectural contracts, algorithms, endpoints, or UI specifications changed:
- If backend routes, models, or query params changed: update [`docs/rest_api.md`](../../../docs/rest_api.md).
- If core geodetic formulas or boundary envelopes changed: update [`docs/core_engine.md`](../../../docs/core_engine.md).
- If webapp component hierarchy, design tokens, archetypes, or client algorithms changed: update the relevant living file in [`docs/webapp/`](../../../docs/webapp/README.md) (`architecture.md`, `design_system.md`, `archetypes.md`, or `actionable_inventory.md`).
- If test suites or E2E browser runner evolved: update [`docs/testing.md`](../../../docs/testing.md).
- If deployment, containerization, or environment variables changed: update [`docs/deployment.md`](../../../docs/deployment.md).

*Note*: Keep changes in `docs/` focused on the latest authoritative state without meta chatter about agents or session numbers.

---

### Step 4: Web Application Living Documentation & README Maintenance (If `webapp/` Touched)

If any files in `webapp/` were modified:

1. **Synchronize Living Web App Documentation ([`docs/webapp/`](../../../docs/webapp/README.md))**:
   - Update `docs/webapp/actionable_inventory.md` if actionable UI elements, triggers, or reactive controls changed.
   - Update `docs/webapp/archetypes.md` if component archetypes, CSS classes, or DOM locators changed.
   - Update `docs/webapp/design_system.md` if colors, layout dimensions, or theme tokens changed.
   - Update `docs/webapp/architecture.md` if client SPA topology, GIS layer mechanics, or i18n changed.


2. **Update [`webapp/README.md`](../../../webapp/README.md)**:
   - If files/directories were added, moved, or deleted: update Section 2 directory tree.
   - If architecture, guardrails, or testing commands changed: update Section 3 or Section 4.

---

### Step 5: Test & Quality Gate Check

Verify code health and test suite execution:
1. **Frontend Integrity** (if `webapp/` was touched):
   ```powershell
   npm --prefix webapp test
   npm --prefix webapp run build
   ```
2. **Python Code Hygiene & Tests** (if Python code was touched):
   ```powershell
   D:\anaconda3\envs\pytransdat\python.exe -m ruff check pytransdatro api tests tools
   D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_trans_ro.py
   ```
   *(Run full pytest with `BypassSandbox: true` if API TestClient tests are involved)*.

3. **Report to User**:
   - Provide a concise summary of:
     - Task file recorded in `.agents/tasks/`.
     - Living documentation updated in `docs/`.
     - Changelog & README updates in `webapp/`.
     - Quality gate test results.
   - Proactively suggest running `task-commit` to review and commit changes to Git.
