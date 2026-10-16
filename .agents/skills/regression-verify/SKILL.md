---
name: regression-verify
description: Runs the full geodetic numerical precision suite and API regression tests. Use whenever verifying that code changes or refactoring did not alter coordinate transformation results.
---

# Regression Verification Skill

## Purpose
In geodetic coordinate transformations, "cleaner" or refactored code is useless if it changes mathematical precision. This skill executes the full regression and benchmark test suite to ensure **100% zero-regression** and byte-for-byte fidelity against ANCPI / Romgeo standards.

---

## Verification Pipeline

### Step 1: Core Geodetic Unit & Roundtrip Tests
Run the core transformation test suite:

```powershell
D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_trans_ro.py -v
```

Verifies:
- Stereo 70 to ETRS89 single and batch transformations.
- ETRS89 to Stereo 70 single and batch transformations.
- 2D vs 3D elevation adjustments ($h = H + N$).
- Boundary edge cases (`OutOfGridErr`, `NoDataGridErr`).

---

### Step 2: REST API Contract & Legacy Emulation Tests
Run the API integration tests (or full suite).

> [!IMPORTANT]
> **Windows Sandbox Notice**: Always execute Step 2 (and full `pytest`) with `BypassSandbox: true`. Starlette's `TestClient` uses `anyio` to bridge synchronous test calls to asynchronous FastAPI endpoints using background worker threads and local loopback sockets. Under standard Windows sandbox isolation, local loopback socket communication is blocked and the tests will hang indefinitely.

```powershell
D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_api_modern.py tests/test_api_legacy.py tests/test_api_rate_limiter.py tests/test_telemetry.py -q
```

Verifies:
- Modern JSON endpoints (`/api/v1/transform/*`).
- 100% byte-for-byte backward compatible emulation of the 2012 TransDatOnline Java service (`/transdatonline/cooOpService`).
- Rate limiting and asynchronous SQLite telemetry recording.

---

### Step 3: ANCPI / Romgeo Official Benchmark Comparison
Execute the benchmark comparison tool against reference national control point datasets (can run with `BypassSandbox: false`):

```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/generate_ancpi_comparison_reports.py
```

Criteria for Success:
- Maximum horizontal discrepancy: $< 0.001\text{ mm}$ ($< 10^{-6}\text{ m}$).
- Maximum vertical discrepancy: $< 0.001\text{ mm}$.
- Byte-for-byte conformity on status codes and coordinate strings.

---

### Step 4: Web Application Utility & UI Regression Suite
Execute the frontend unit test suite to verify coordinate formatting, angle parsing, and exporter integrity (can run with `BypassSandbox: false`):

```powershell
npm --prefix webapp test
```

Verifies:
- DMS format conversions and cardinal notations (`dmsParser.test.ts`).
- Single and batch delimiter auto-detection (`pointParser.test.ts`).
- Romania envelope and boundary threshold tests (`gridBounds.test.ts`).
- CSV, GeoJSON, and AutoCAD DXF serialization (`exporters.test.ts`).
- Accessible UI factory functions and keyboard navigation (`ui.test.ts`).

---

## Reporting
If any test fails or any numerical discrepancy exceeds machine epsilon, the change must be reverted or halted immediately until the regression is resolved.

