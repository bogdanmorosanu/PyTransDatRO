# PyTransDatRO Project Invariants & Agent Workflow

## Core Architecture Invariants
- **Pure Python Core Geodetic Library (`pytransdatro/`)**:
  - The core geodetic library is implemented strictly in pure Python (standard library only), with no external scientific C/compiled dependencies (e.g., NumPy, SciPy) for processing.
  - **Exception: `numpy` is allowed STRICTLY for reading and unpickling the binary grid files.** All data must be processed using pure Python.
  - Core data processing and mathematical pipelines in `pytransdatro/` must NEVER import from or depend on the web API layer (`api/`) or web frameworks.
- **Web API Layer (`api/`)**:
  - The web service layer is built on top of `pytransdatro` using standard modern web libraries (`fastapi`, `uvicorn`, `pydantic`, `python-multipart`, `httpx`).
  - These web dependencies are optional extras (`pytransdatro[api]`) and must never be imported inside `pytransdatro/`.
- **Zero Functional Regression Rule**:
  - Coordinate transformation numerical precision must remain byte-for-byte identical to official ANCPI/Romgeo benchmarks ($< 10^{-9}$ tolerance). Always verify using the `regression-verify` skill.
- **Binary Grid Files**: Do not modify the loading structure or content of the binary grid files in `grids/`.
- **Python Environment**: Always use the `pytransdat` conda environment located at `D:\anaconda3\envs\pytransdat` (Python executable: `D:\anaconda3\envs\pytransdat\python.exe`) for running all commands, tests, and scripts in this project.

## Development & Continuous Quality Workflow

### 1. Continuous Development & Maintenance Mode
- **Completed Refactor Baseline**: The multi-phase refactoring roadmap ([`refactor_roadmap.md`](refactor_roadmap.md)) is **100% COMPLETED** and serves as the immutable architectural foundation.
- The AI agent is fully prepared for general development (new features, optimizations, bug fixes, UI enhancements, and ongoing code hygiene) without being constrained by the completed roadmap phases.
- **Permanent Quality Standards**: Every subsequent code change must uphold the patterns and invariants established during the refactor:
  - Standard library purity in `pytransdatro/` (zero web dependencies, no scientific C-extensions).
  - Explicit error handling (zero bare `except: pass`; catch specific exceptions and log diagnostics).
  - Universal Python typing and reST/Sphinx docstrings with mathematical formulas.
  - Centralized configuration hierarchy (`constants.py`, `api/config.py`, `webapp/src/config.ts`, `.env.example`).
  - Modular CSS tokens and reusable UI factories in `webapp/src/ui/`.

### 2. Layer Discipline & Test Accompaniment
- **Layer Allocation**:
  - Geodetic math / spatial calculations $\to$ `pytransdatro/` (pure Python standard library).
  - Web endpoints / serialization / validation $\to$ `api/` (FastAPI / Pydantic).
  - User interface / interactions $\to$ `webapp/` (TypeScript / Leaflet).
- **Test Accompaniment**: Every new feature or bugfix must be accompanied by tests in `tests/` (`pytest`) and/or `webapp/src/` (`vitest`).
- **Config & Constants**: Never hardcode operational parameters. Place invariants in `pytransdatro/constants.py` and deployment knobs in `.env` / `api/config.py` / `webapp/src/config.ts`.

### 3. Windows Sandbox & Test Execution Rules
- **FastAPI / Starlette API Tests (`BypassSandbox: true`)**:
  - Starlette's `TestClient` uses `anyio` to bridge synchronous test calls to asynchronous FastAPI endpoints via background worker threads and local loopback sockets.
  - Under the Windows sandbox isolation (`BypassSandbox: false`), local loopback socket communication is blocked, causing tests importing `TestClient` or running full `pytest` to hang indefinitely.
  - **Rule**: Whenever running the full Python test suite (`pytest`) or any test module importing `TestClient` (`tests/test_api_*.py`, `tests/test_telemetry.py`), always run with `BypassSandbox: true`.
- **Standard Sandboxed Commands (`BypassSandbox: false`)**:
  - Core geodetic tests (`tests/test_trans_ro.py`), config tests (`tests/test_config.py`), ANCPI benchmark scripts, Ruff linting/formatting, and frontend Vitest unit tests (`npm --prefix webapp test`) run safely in standard sandbox mode.

### 4. Agent Rules (`.agents/rules/`)
Consult and follow the specialized rules across all development tasks:
- [`code-quality-standards.md`](.agents/rules/code-quality-standards.md) (layer decoupling, typing, error handling, formatting, Vitest)
- [`configuration-standards.md`](.agents/rules/configuration-standards.md) (invariants vs. `.env` settings)
- [`python-style.md`](.agents/rules/python-style.md) (typing, docstrings, no silent `pass`)
- [`webapp-design.md`](.agents/rules/webapp-design.md) (modular CSS, UI factories, no ad-hoc styling)

### 5. Standard Agent Tooling & Verification Workflow
At the end of any significant coding session (feature, bugfix, or maintenance):
1. **Audit**: Run the `code-audit` skill to ensure zero lint errors, syntax issues, or boundary breaches.
2. **Verify**: Run the `regression-verify` skill to ensure zero regressions against benchmarks, Python tests, and Vitest suite.
3. **Document**:
   - Use `doc-record` to update living documentation in `docs/` if architectural contracts or APIs changed.
   - Use `task-record` to document the chat session's outcome in `.agents/tasks/`.
4. **Commit**: Use `task-commit` to review diffs and create conventional git commits.
