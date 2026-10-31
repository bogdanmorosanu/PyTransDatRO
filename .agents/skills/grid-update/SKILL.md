---
name: grid-update
description: Orchestrates the complete lifecycle of updating or replacing the official Romanian geodetic transformation grid (.spg). Guides the operator through candidate pre-flight inspection, reminds and verifies official ANCPI ground truth test dataset refresh (tests/data/*_expected.csv), executes safe archiving, installs the new grid, regenerates metadata JSON, updates benchmark reports, and verifies quality gates. Use whenever asked to update, replace, or upgrade an ANCPI .spg grid file.
---

# Official Grid Update & Lifecycle Skill (`grid-update`)

## Purpose
This skill orchestrates the end-to-end procedure for updating or replacing the official Romanian 3D binary transformation grid (`.spg` file released by ANCPI / CNC).

It enforces key architectural invariants, ensures zero stale remnants of previous grids remain in the API or web application, and actively reminds the user/operator that **updating the grid file requires refreshing the official ANCPI reference test data (`tests/data/*_expected.csv`)** to maintain valid benchmark reports and passing regression tests.

---

## Operating Principles & Golden Rules

1. **The Dual-Asset Invariant (Grid + Ground Truth Test Data)**:
   > [!IMPORTANT]
   > Replacing an `.spg` grid file alters coordinate transformation outputs. The 9,046 national control points in [`tests/data/`](../../../tests/data/) **will not match** old expected coordinates.
   > The agent must always ask or verify if the operator has generated new expected coordinates using the official ANCPI desktop software (TransDatRO / RomGEO Table Convert) on the 9,046 test points.

2. **Single Active Grid Invariant**:
   [`pytransdatro/grids/`](../../../pytransdatro/grids/) must contain **exactly one** active `.spg` file. The previous grid must be moved to [`pytransdatro/grids/archive/`](../../../pytransdatro/grids/archive/) to prevent [`AmbiguousGridError`](../../../pytransdatro/exceptions.py).

3. **Unified CLI Orchestrator**:
   All operations should be executed through the unified tool [`tools/update_grid.py`](../../../tools/update_grid.py) rather than individual legacy scripts.

4. **In-Memory Server Lifecycle**:
   FastAPI (`api/main.py`) caches the deserialized grid singleton in RAM at startup. Replacing the file on disk requires restarting the backend runtime (`uvicorn` or Docker) to take effect.

---

## Workflow Steps

### Step 1: Pre-Flight Audit & User Engagement

1. **Audit Current Environment**:
   Run the status check to inspect current grid provenance, metadata synchronization, and test data timestamps:
   ```powershell
   D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py status
   ```

2. **Verify Candidate Grid (`inspect`)**:
   Validate the candidate `.spg` file before modifying repository files:
   ```powershell
   D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py inspect <path/to/candidate.spg>
   ```
   - Check that `[Integrity Evaluation]` reports `[PASS]`.
   - Review side-by-side differences with the currently active grid (release version, date, array dimensions, nodata node counts, and Helmert parameter deltas).

3. **Prompt the User Regarding Official ANCPI Test Data**:
   > [!WARNING]
   > Ask the user explicitly:
   > *"Have you generated updated expected reference datasets (`tests/data/etrs89_to_st70_expected.csv` and `tests/data/st70_to_etrs89_expected.csv`) using the official ANCPI application with this new grid?"*

   - **If YES**: Proceed to apply both the grid and the updated expected test files.
   - **If NO**: Explain that the grid can be installed and metadata extracted, but running the benchmark or `pytest tests/test_trans_ro.py` will report failures until the official ANCPI reference coordinates are generated from:
     - Input 1: [`tests/data/etrs89_to_st70_input.csv`](../../../tests/data/etrs89_to_st70_input.csv) $\to$ produce `tests/data/etrs89_to_st70_expected.csv`
     - Input 2: [`tests/data/st70_to_etrs89_input.csv`](../../../tests/data/st70_to_etrs89_input.csv) $\to$ produce `tests/data/st70_to_etrs89_expected.csv`

---

### Step 2: Apply the New Grid

Apply the candidate grid using the unified CLI tool:
```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py apply <path/to/candidate.spg>
```

This automatically:
- Moves the previous `.spg` file into [`pytransdatro/grids/archive/`](../../../pytransdatro/grids/archive/).
- Installs the new `.spg` file into [`pytransdatro/grids/`](../../../pytransdatro/grids/).
- Extracts 100% of embedded parameters, layers, and statistics into [`pytransdatro/grids/spg_metadata.json`](../../../pytransdatro/grids/spg_metadata.json).
- Checks bounding extents against client pre-flight constants in [`webapp/src/utils/gridBounds.ts`](../../../webapp/src/utils/gridBounds.ts) and warns if coverage expands.

If updated test datasets were provided by the user, place them into [`tests/data/`](../../../tests/data/).

---

### Step 3: Run ANCPI Replication Benchmarks

Execute the replication benchmark across the 9,046 national control points:
```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py benchmark
```

- Compares PyTransDatRO forward and inverse calculations against `tests/data/*_expected.csv`.
- Verifies that max residuals do not exceed $0.003\text{ m}$ (or $0.00003''$) across all 9,046 points.
- Writes active grid provenance (`grid_file`, `grid_version`, `grid_sha256`) and overall status (`PASS`) into [`pytransdatro/grids/ancpi_comparison_summary.json`](../../../pytransdatro/grids/ancpi_comparison_summary.json).
- Regenerates summary text files, CSV residual tables, and interactive Leaflet HTML maps in `reports/`.

*(Note: If official test datasets have not yet been refreshed, use `--allow-fail` during interim inspection).*

---

### Step 4: Quality & Regression Verification

Run the full verification suite across Python and TypeScript layers:
```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py verify
```

1. **Python Core & API (`pytest`)**:
   - [`tests/test_trans_ro.py`](../../../tests/test_trans_ro.py) (54 geodetic unit, roundtrip, boundary, and nodata tests).
   - [`tests/test_api_modern.py`](../../../tests/test_api_modern.py) (API endpoints, metadata, bounding boxes, streaming file).
2. **Web Frontend (`vitest`)**:
   - `npm --prefix webapp test` (123 parser, formatting, and UI factory tests).

If web app panels, UI factories, or user flows were modified, execute the Playwright browser E2E test suite:
```powershell
npm --prefix webapp run test:e2e
```

---

### Step 5: Service Restart & Verification

Inform the operator to restart runtime services:
- **Local Dev Server**: Restart Uvicorn on port `8000` so `app.state.trans` instantiates the new grid.
- **Production / Docker**:
  ```bash
  docker compose restart backend
  ```

Verify in the web application:
1. **Grid Metadata Panel**: Opens and displays new grid filename, release version, and layer statistics.
2. **Geodetic Help Panel**: Displays the new baseline grid chip in the benchmark header and confirms `PASS` status.

---

### Step 6: Disaster Recovery & Rollback

If the new grid causes regressions or needs to be rolled back:
```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/update_grid.py rollback
```
This restores the most recently archived `.spg` file from `pytransdatro/grids/archive/`, regenerates `spg_metadata.json`, and returns the system to the previous verified baseline.

---

## Authoritative Documentation Reference
- Living Architecture & Lifecycle Guide: [`docs/grid_lifecycle.md`](../../../docs/grid_lifecycle.md)
- Grid Specification Reference: [`docs/grid.md`](../../../docs/grid.md)
- Testing & Quality Strategy: [`docs/testing.md`](../../../docs/testing.md)
