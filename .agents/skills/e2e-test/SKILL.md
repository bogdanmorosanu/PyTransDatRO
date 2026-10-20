---
name: e2e-test
description: Executes and debugs the Playwright End-to-End (E2E) browser test suite for PyTransDatRO. Covers headless execution, interactive UI runner mode, targeted spec runs, and diagnostic failure inspection using test-e2e.ps1 and npm scripts. Use whenever validating browser flows, UI regressions, or Playwright specs.
---

# Playwright E2E Browser Testing Skill (`e2e-test`)

## Purpose
This skill governs the execution, debugging, and verification of the Playwright End-to-End (E2E) browser test suite for the PyTransDatRO web application (`webapp/e2e/`). It ensures complete confidence in user-facing flows—single-point conversions, batch text parsing, multi-file drag-and-drop transformations, Web GIS map centering and zoom constraints, bi-directional internationalization (`RO` / `EN`), and API error handling.

---

## Architecture & Test Isolation

The E2E test harness (`webapp/playwright.config.ts`) is designed for 100% deterministic, offline-capable execution:

1. **Automated WebServer Orchestration**:
   - Compiles the frontend into an isolated build folder (`webapp/dist-e2e/`) with `VITE_E2E=true`.
   - Boots an isolated FastAPI backend server on **port 8100** (`http://127.0.0.1:8100`) using the project Conda environment (`D:\anaconda3\envs\pytransdat\python.exe`).
   - Configures temporary SQLite telemetry logs and feedback files in `webapp/e2e/.tmp/` to avoid polluting development logs.
   - Sets rate limits to unlimited (`100000`) to prevent false rate-limit throttling during rapid test execution.
2. **Offline Fixtures & Request Interception (`webapp/e2e/fixtures.ts`)**:
   - Stubs external map tile requests (`tile.openstreetmap.org`, CartoDB, ESRI) with a transparent 1x1 PNG to prevent external network flakiness.
   - Stubs external Google Fonts (`fonts.googleapis.com`, `fonts.gstatic.com`) with empty CSS responses.
3. **Diagnostic Test API (`MapView.ts`)**:
   - In `VITE_E2E=true` builds, `MapView.ts` exposes a read-only `window.__e2e.mapView()` diagnostic snapshot. E2E tests query map center, zoom level, marker counts, and layer state directly rather than relying on fragile DOM scraping.

---

## Test Execution Modes

### 1. Full Headless Suite (Standard Automated Verification)

Run the full suite of all 19+ Chromium browser tests in headless mode:

Using the PowerShell runner from repository root:
```powershell
.\test-e2e.ps1
```

Or directly via npm:
```powershell
npm --prefix webapp run test:e2e
```

---

### 2. Interactive Visual UI Mode (Time-Travel Debugging)

For step-by-step visual debugging, DOM inspection, timeline scrubbing, and locator testing:

Using the PowerShell runner:
```powershell
.\test-e2e.ps1 -Ui
```

Or directly via npm:
```powershell
npm --prefix webapp run test:e2e:ui
```

*Note*: UI mode launches the interactive Playwright Test Runner GUI.

---

### 3. Targeted Spec & Focused Execution

To run a single spec file or a specific test pattern:

1. **Run a Specific Test Spec File**:
   ```powershell
   npm --prefix webapp exec playwright test e2e/smoke.spec.ts
   npm --prefix webapp exec playwright test e2e/point.spec.ts
   npm --prefix webapp exec playwright test e2e/batch.spec.ts
   npm --prefix webapp exec playwright test e2e/file.spec.ts
   npm --prefix webapp exec playwright test e2e/map.spec.ts
   npm --prefix webapp exec playwright test e2e/i18n.spec.ts
   npm --prefix webapp exec playwright test e2e/errors.spec.ts
   ```

2. **Run by Test Name Regex Filter (`-g`)**:
   ```powershell
   npm --prefix webapp exec playwright test -g "handles out of grid coordinates"
   npm --prefix webapp exec playwright test -g "language switch"
   ```

3. **Run with Visible Headed Browser**:
   ```powershell
   npm --prefix webapp exec playwright test e2e/point.spec.ts --headed
   ```

---

### 4. Diagnostic Inspection on Failure

When a test fails:

1. **HTML Test Report**:
   Playwright produces an interactive HTML report. View it with:
   ```powershell
   npm --prefix webapp exec playwright show-report
   ```

2. **Traces & Screenshots**:
   - Check `webapp/test-results/` for failure screenshots and Playwright execution traces (`trace.zip`).
   - Open any trace for complete step replay:
     ```powershell
     npm --prefix webapp exec playwright show-trace webapp/test-results/<test-dir>/trace.zip
     ```

---

## Invariants & Rules

- **Zero Network Flakiness**: All external assets must remain intercepted or mocked in `e2e/fixtures.ts`.
- **E2E Accompaniment Rule**: Any change to covered UI flows (single-point transform, batch paste, file upload/download, map centering/no-auto-pan, language switch, API error display) or `.agents/rules/webapp-design.md` must run or update these Playwright specs.
- **Fast Clean Cleanup**: The runner automatically kills the temporary backend process and cleans up temporary test servers on completion.
