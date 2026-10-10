---
name: webapp-init
description: Initializes the webapp development environment, verifies/starts FastAPI backend and Vite dev servers, launches the browser at http://localhost:5173, and primes the agent context with the latest web map status, changelog history, and API contracts. Use at the start of any conversation focusing on webapp or web map development.
---

# Webapp Init Skill

## Purpose
This skill primes the agent and development environment for working on the **PyTransDatRO Web Application** (`webapp/`), with a particular focus on the interactive Leaflet Web GIS map interface (`MapView.ts`), its component interactions, and its contract with the FastAPI backend (`api/`).

---

## Operating Principles

1. **Automate the Dev Environment**:
   - Ensure the developer doesn't have to manually launch servers or navigate to the browser.
   - Respect running processes (do not duplicate servers if already active on the designated ports).
2. **Fast Context Priming**:
   - Give the agent an immediate mental model of:
     - The latest web map capabilities and components.
     - Recent evolution via [webapp/CHANGELOG.md](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/CHANGELOG.md).
     - Single source of truth status via [webapp/README.md](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/README.md).
     - The underlying FastAPI web API contracts and endpoints (`api/`).
3. **Preserve Architectural Invariants**:
   - Respect pure Python core constraints (`pytransdatro/` has no web dependencies).
   - Keep API logic in `api/` and client logic in `webapp/`.

---

## Workflow Steps

### Step 1: Check and Start Services

Execute the following checks via `run_command` in PowerShell:

1. **Verify Backend (Port 8000)**:
   - Test if the FastAPI backend is running:
     ```powershell
     try { (Invoke-WebRequest -Uri "http://127.0.0.1:8000/health" -UseBasicParsing -TimeoutSec 2).StatusCode } catch { 0 }
     ```
   - If not running (returns 0 or errors):
     - Propose running the backend as a background daemon process using the project conda environment:
       - **CommandLine**: `D:\anaconda3\envs\pytransdat\python.exe -m uvicorn api.main:app --host 127.0.0.1 --port 8000`
       - **Cwd**: `D:\ProiecteRealizate\pyTransDatRO\repo\PyTransDatRO`
       - **IsDaemon**: `true`
       - **WaitMsBeforeAsync**: `2000`

2. **Verify Frontend (Port 5173)**:
   - Test if Vite dev server is running:
     ```powershell
     try { (Invoke-WebRequest -Uri "http://localhost:5173" -UseBasicParsing -TimeoutSec 2).StatusCode } catch { 0 }
     ```
   - If not running:
     - Ensure dependencies are present: verify `webapp/node_modules/` exists (run `npm install` in `webapp` if missing).
     - Launch Vite dev server as a background daemon:
       - **CommandLine**: `npm run dev`
       - **Cwd**: `D:\ProiecteRealizate\pyTransDatRO\repo\PyTransDatRO\webapp`
       - **IsDaemon**: `true`
       - **WaitMsBeforeAsync**: `2000`

3. **Open Browser**:
   - Open the application in the default web browser:
     ```powershell
     Start-Process "http://localhost:5173"
     ```

---

### Step 2: Ingest Webapp Status & History

To rapidly get up to speed without reading thousands of lines of code:

1. **Read Latest Changelog**:
   - Read the top section of [`webapp/CHANGELOG.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/CHANGELOG.md) to understand the most recent features, bug fixes, and modifications.
2. **Read Webapp Status**:
   - Review [`webapp/README.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/README.md) to review the component structure and current design system.
3. **Inspect Map State**:
   - Check [`webapp/src/components/MapView.ts`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/src/components/MapView.ts) for:
     - Active layers, basemaps (OSM, CartoDB Positron, CartoDB Dark Matter, ESRI Imagery).
     - Grid extent polygon / bounding box overlay.
     - Marker cluster behavior and point rendering.
     - Event bus dispatchers/listeners (`point:transformed`, `batch:transformed`, `map:click`).
   - Check [`webapp/src/style/map.css`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/src/style/map.css) for custom layer styling, popups, and cluster visuals.

---

### Step 3: Map API Contracts & Capabilities

Understand the backend API bridge behind the scenes:

1. **Frontend REST Client**:
   - Examine [`webapp/src/utils/apiClient.ts`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/src/utils/apiClient.ts) for the data shapes and active operations (`Stereo70ToETRS89`, `ETRS89ToStereo70`).
2. **Backend Endpoints (`api/routes/`)**:
   - Single point: [`api/routes/point.py`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/api/routes/point.py) (`POST /api/v1/transform/point`)
   - Batch transformation: [`api/routes/batch.py`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/api/routes/batch.py) (`POST /api/v1/transform/batch`)
   - File upload: [`api/routes/file.py`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/api/routes/file.py) (`POST /api/v1/transform/file`)
   - Grid & Telemetry info: [`api/routes/info.py`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/api/routes/info.py) (`GET /api/v1/info/grid`, `GET /api/v1/info/stats`)
3. **Live Documentation & Schema Inspection**:
   - OpenAPI is proxied by Vite at `http://localhost:5173/openapi.json` and `http://localhost:5173/docs` (served by FastAPI at `http://127.0.0.1:8000/docs`).
   - When new features require additional data or endpoints (e.g. raster shift layers, contour tiles, polygon bounding checks), they should be added in `api/routes/` and typed in `api/schemas/`.

---

### Step 4: Present Brief Status Summary to User

Provide a clean, bulleted status report to the user:
- **Environment Status**: Backend URL (`http://127.0.0.1:8000`), Frontend URL (`http://localhost:5173`), Browser opened.
- **Web Map State**: Current basemaps, clustering status, grid extent overlay, and guardrails (e.g., 2,000 point warning modal).
- **Recent Evolution**: Brief 2-3 bullet point summary from `webapp/CHANGELOG.md`.
- **API Connectivity**: Available endpoints and schemas ready for integration.
- **Ready**: Ask the user what feature, tweak, or bug fix they would like to tackle.
