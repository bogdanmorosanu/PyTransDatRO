# Telemetry Out-of-Grid Validation and Logging Timing Fix

## Objectives
- [x] Investigate and resolve telemetry database pollution where out-of-grid coordinates (e.g., false Northing/Easting values of 500,000 passed to ETRS89 endpoints) were recorded as spatial cells in `usage.db`.
- [x] Align execution timing in `TransRO.etrs89_to_st70` with `TransRO.st70_to_etrs89` so telemetry logging executes strictly after 2D grid boundary validation (`_t_gr2d.trans`).
- [x] Preserve pure Python standard library invariants, dynamic grid upgradability (no hardcoded coordinate bounding boxes), and zero latency impact on calculation threads.
- [x] Expand telemetry test suite (`tests/test_telemetry.py`) with comprehensive failure isolation, scalar/batch out-of-grid checks, and mixed-batch regression tests for both transformation directions.

## Implementation Details

1. **Telemetry Logging Timing Relocation (`pytransdatro/trans_ro.py`)**:
   - In `TransRO.etrs89_to_st70`, relocated `self._collector.record(...)` from the method entry point (previously lines 121–126) down to immediately after `self._t_gr2d.trans(r_n, r_e, 1)` and prior to the 1D elevation step.
   - **Architectural Rationale**: In the forward pipeline (`st70_to_etrs89`), telemetry was already recorded downstream of `_t_gr2d.trans`. Moving it in `etrs89_to_st70` achieves complete pipeline symmetry: invalid or out-of-grid coordinates abort execution via `OutOfGridErr` before entering the asynchronous telemetry queue.
   - **Dynamic Grid Coverage**: Relies directly on `Grid2D.covered_by_grid()` dynamically evaluated against the active `.spg` grid file rather than duplicating hardcoded coordinate thresholds.

2. **Expanded Telemetry Unit Testing (`tests/test_telemetry.py`)**:
   - `test_telemetry_etrs89_to_st70_out_of_grid_not_recorded`: Validates that scalar 2D (`0.0, 0.0`), scalar 3D, full out-of-grid batches, and mixed batches (containing both in-grid and out-of-grid coordinates) raise `OutOfGridErr` and record 0 entries to `usage` and `daily_stats`.
   - `test_telemetry_st70_to_etrs89_out_of_grid_not_recorded`: Validates that scalar 2D (`5.0, 1.0`), scalar 3D (`0.0, 0.0`), boundary outliers (`100000.0, 100000.0` / `900000.0, 900000.0`), and mixed batches raise `OutOfGridErr` and record 0 entries to SQLite.
   - `test_telemetry_in_grid_both_directions_recorded`: Verifies that valid in-grid coordinates in both directions (2D, 3D, and batch) properly record entries with `is_st70_to_etrs89` set to 0 and 1.
   - `test_telemetry_failure_isolation`: Verifies that logger failures (simulated disk errors) are caught and suppressed in both `etrs89_to_st70` and `st70_to_etrs89` without failing coordinate mathematics.

## Verification & Testing
- **Test Modules**: `tests/test_telemetry.py`, `tests/test_trans_ro.py`, `tests/test_api_modern.py`, `tests/test_api_legacy.py`, `tests/test_api_rate_limiter.py`
- **Scenarios Covered**:
  - Scalar and batch out-of-grid rejection without database writes.
  - Mixed-batch atomic exception propagation without partial logging.
  - Failure isolation for telemetry collector errors across both transformation directions.
  - Full roundtrip spatial cell and daily volume stats integrity.
- **Results**:
  - `tests/test_telemetry.py`: 5/5 passed in 0.27s.
  - Full test suite: 76/76 passed in 4.09s (100% pass rate).
