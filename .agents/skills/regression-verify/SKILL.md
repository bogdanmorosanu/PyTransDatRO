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
Run the API integration tests:

```powershell
D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_api_modern.py tests/test_api_legacy.py tests/test_api_rate_limiter.py tests/test_telemetry.py -q
```

Verifies:
- Modern JSON endpoints (`/api/v1/transform/*`).
- 100% byte-for-byte backward compatible emulation of the 2012 TransDatOnline Java service (`/transdatonline/cooOpService`).
- Rate limiting and asynchronous SQLite telemetry recording.

---

### Step 3: ANCPI / Romgeo Official Benchmark Comparison
Execute the benchmark comparison tool against reference national control point datasets:

```powershell
D:\anaconda3\envs\pytransdat\python.exe tools/generate_ancpi_comparison_reports.py
```

Criteria for Success:
- Maximum horizontal discrepancy: $< 0.001\text{ mm}$ ($< 10^{-6}\text{ m}$).
- Maximum vertical discrepancy: $< 0.001\text{ mm}$.
- Byte-for-byte conformity on status codes and coordinate strings.

---

## Reporting
If any test fails or any numerical discrepancy exceeds machine epsilon, the refactor must be reverted or halted immediately until the mathematical deviation is identified.
