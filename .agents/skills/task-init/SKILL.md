---
name: task-init
description: Initializes any PyTransDatRO development task or chat session. Discovers user intent, dynamically indexes living architecture docs (docs/), UI design guide, and task history (../../tasks/), verifies/launches FastAPI and Vite background services, and outputs a 4-bullet readiness briefing. Use at the start of any new session or task.
---

# Task Initialization Skill (`task-init`)

## Purpose
This skill serves as the **universal entrypoint** for all development, refactoring, bug fixing, deployment, and testing sessions in the PyTransDatRO repository. It dynamically primes the agent context with authoritative architecture guides and session history, verifies and orchestrates required background runtimes (Python FastAPI backend, Vite frontend dev server), and outputs a high-signal readiness briefing.

---

## Operating Principles

1. **Universal Applicability**:
   - Adapts to any task domain: core geodetic engine, REST API, web application, telemetry/logging, testing/benchmarking, Docker deployment, CI/CD pipelines, or isolated scripts.
2. **Dynamic Context Priming (Zero Guesswork)**:
   - Queries living architectural documentation (`docs/*.md`) and [`docs/webapp/README.md`](../../../docs/webapp/README.md) before proposing code edits.
   - Inspects existing task files in [`.agents/tasks/`](../../tasks/) to understand precedent design decisions and avoid duplicate work.
3. **Non-Disruptive Service Orchestration**:
   - Tests ports `8000` (FastAPI) and `5173` (Vite) before spawning daemons. Never spawns duplicate processes if servers are already active.
   - Always uses the official Conda Python executable: `D:\anaconda3\envs\pytransdat\python.exe`.

---

## Workflow Steps

### Step 1: Discover Task Intent & Subsystem Classification

1. **Analyze the Request**:
   - Read the user prompt or `<original_task>` to classify the affected domain(s):
     - **Core Geodetic Math & Grids**: `pytransdatro/` (pure Python standard library).
     - **Web API Layer**: `api/` (FastAPI, Pydantic, Uvicorn).
     - **Web Application & Web GIS**: `webapp/` (TypeScript, Leaflet, modular CSS).
     - **Tests & Quality Assurance**: `tests/`, `tools/`, `webapp/e2e/`.
     - **Deployment & Infrastructure**: `Dockerfile`, `docker-compose.yml`, GitHub Actions, environment configuration.
     - **Cross-Layer or Unexpected Tasks**: Docker containerization, log rotation, performance profiling, standalone utilities.

2. **Index Living Documentation (`docs/`)**:
   - Consult the authoritative documentation mapped to the classified domain:
     - **High-Level System & Boundaries**: [`docs/architecture.md`](../../../docs/architecture.md)
     - **Core Mathematical Algorithms**: [`docs/core_engine.md`](../../../docs/core_engine.md)
     - **Binary Grid Models & Shifts**: [`docs/grid.md`](../../../docs/grid.md)
     - **REST Endpoints & Schemas**: [`docs/rest_api.md`](../../../docs/rest_api.md)
     - **Web Application Architecture**: [`docs/webapp/architecture.md`](../../../docs/webapp/architecture.md)
     - **Web App Master Index & Guide**: [`docs/webapp/README.md`](../../../docs/webapp/README.md)
     - **UI Design System & Tokens**: [`docs/webapp/design_system.md`](../../../docs/webapp/design_system.md)
     - **Component Archetypes**: [`docs/webapp/archetypes.md`](../../../docs/webapp/archetypes.md)
     - **Actionable UI Elements Inventory**: [`docs/webapp/actionable_inventory.md`](../../../docs/webapp/actionable_inventory.md)
     - **Usage Telemetry Subsystem**: [`docs/telemetry.md`](../../../docs/telemetry.md)
     - **Application Logging Infrastructure**: [`docs/log.md`](../../../docs/log.md)
     - **Testing Strategy & Benchmarks**: [`docs/testing.md`](../../../docs/testing.md)
     - **Deployment & Production Runtime**: [`docs/deployment.md`](../../../docs/deployment.md)

3. **Inspect Task History ([`.agents/tasks/`](../../tasks/))**:
   - List the recent files in `.agents/tasks/` to understand the latest state, roadmaps, and patterns.

---

### Step 2: Service & Environment Orchestration

When a task involves the web API, web frontend, or browser E2E workflows, check whether backend and frontend servers are listening.

1. **Verify Backend Service (Port 8000)**:
   - Check if FastAPI is running:
     ```powershell
     try { (Invoke-WebRequest -Uri "http://127.0.0.1:8000/health" -UseBasicParsing -TimeoutSec 2).StatusCode } catch { 0 }
     ```
   - If response is `0` or fails:
     - Propose launching FastAPI backend as a daemon:
       - **CommandLine**: `D:\anaconda3\envs\pytransdat\python.exe -m uvicorn api.main:app --host 127.0.0.1 --port 8000`
       - **Cwd**: `.`
       - **IsDaemon**: `true`
       - **WaitMsBeforeAsync**: `2000`

2. **Verify Frontend Dev Server (Port 5173)**:
   - Check if Vite is running:
     ```powershell
     try { (Invoke-WebRequest -Uri "http://localhost:5173" -UseBasicParsing -TimeoutSec 2).StatusCode } catch { 0 }
     ```
   - If response is `0` or fails (and task touches frontend or UI):
     - Ensure dependencies exist in `webapp/node_modules/` (run `npm --prefix webapp install` if missing).
     - Propose launching Vite as a daemon:
       - **CommandLine**: `npm run dev`
       - **Cwd**: `webapp`
       - **IsDaemon**: `true`
       - **WaitMsBeforeAsync**: `2000`

3. **Browser Launch (Optional for UI tasks)**:
   - If interactive visual inspection is required by the user, launch the default browser:
     ```powershell
     Start-Process "http://localhost:5173"
     ```

---

### Step 3: Present 4-Bullet Readiness Briefing

Conclude initialization by providing a clean, authoritative 4-bullet briefing:

- **1. Intent & Domain Classification**: Identified task objective, target subsystems (`pytransdatro/`, `api/`, `webapp/`, `tests/`, `docs/`, or infra), and boundary invariants.
- **2. Documentation & History Primed**: Authoritative documents indexed (e.g. [`docs/webapp/README.md`](../../../docs/webapp/README.md), [`docs/webapp/architecture.md`](../../../docs/webapp/architecture.md)) and relevant prior tasks reviewed.
- **3. Environment & Runtime Status**: Status of Conda Python (`D:\anaconda3\envs\pytransdat\python.exe`), backend API (`http://127.0.0.1:8000`), and Vite frontend (`http://localhost:5173`).
- **4. Immediate Action Plan**: Exact next technical step to execute the planned work.
