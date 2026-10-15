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

## Development & Refactoring Workflow Modes

### 1. Refactoring & Code Quality Mode
- When working on planned code improvements or cleanup, consult and follow [`refactor_roadmap.md`](refactor_roadmap.md).
- Stick strictly to the active phase's scope. Do not alter unrelated subsystems or jump ahead to future phases.

### 2. New Features & Bugfix Mode
- **Layer Discipline**: If introducing a new capability, implement it in the correct layer:
  - Geodetic math / spatial calculations $\to$ `pytransdatro/` (pure Python).
  - Web endpoints / serialization / validation $\to$ `api/` (FastAPI / Pydantic).
  - User interface / interactions $\to$ `webapp/` (TypeScript / Leaflet).
- **Test Accompaniment**: Every new feature or bugfix must be accompanied by tests in `tests/` (`pytest`) or `webapp/src/utils/` (`vitest`).
- **Config & Constants**: Never hardcode operational parameters. Place constants in `constants.py` and deployment knobs in `.env` / `api/config.py`.

### 3. Agent Rules (`.agents/rules/`)
Consult and follow the specialized rules regardless of whether you are refactoring or building new features:
- [`code-quality-standards.md`](.agents/rules/code-quality-standards.md) (layer decoupling, typing, error handling, formatting)
- [`configuration-standards.md`](.agents/rules/configuration-standards.md) (invariants vs. `.env` settings)
- [`python-style.md`](.agents/rules/python-style.md) (typing, docstrings, no silent `pass`)
- [`webapp-design.md`](.agents/rules/webapp-design.md) (modular CSS, UI factories, no ad-hoc styling)

### 4. Standard Agent Tooling & Verification Workflow
At the end of any significant coding session (refactor, bugfix, or new feature):
1. **Audit**: Run the `code-audit` skill to ensure zero lint errors, syntax issues, or boundary breaches.
2. **Verify**: Run the `regression-verify` skill to ensure zero regressions against benchmarks and existing tests.
3. **Document**:
   - Use `doc-record` to update living documentation in `docs/` if architectural contracts or APIs changed.
   - Use `task-record` to document the chat session's outcome in `.agents/tasks/`.
4. **Commit**: Use `task-commit` to review diffs and create conventional git commits.
