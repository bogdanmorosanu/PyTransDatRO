# Task 06: Test Suite Refactoring, Round-Trip Transformations, and No-Data Grid Validation

## Objectives
- [x] Investigate and fix the no-data grid edge-case tests (`test_st70_to_etrs89_no_grid_inner_no_data` and `test_st70_to_etrs89_just_out_ro_border_no_data`) on the unified `.spg` grid.
- [x] Identify root cause of failure: `.spg` grids store unpopulated shift nodes as IEEE `NaN` rather than `999.0`, causing `_grid_vs_no_data()` in `pytransdatro/trans_grid.py` to evaluate bicubic splines on NaNs instead of raising `NoDataGridErr`.
- [x] Update `_grid_vs_no_data()` to check `math.isnan(x)` first (for SPG performance) and `x == self.no_data` second.
- [x] Document the derivation rules for the no-data test coordinates in `tests/test_trans_ro.py`.
- [x] Remove discrete hardcoded benchmark points `P1` through `P9` (and their associated dictionaries, fixtures, and parameterized tests).
- [x] Introduce two high-precision forward-and-back round-trip transformation tests (`Stereo70 -> ETRS89 -> Stereo70` and `ETRS89 -> Stereo70 -> ETRS89`).
- [x] Determine and document numerical tolerances:
  - Planar tolerance: `0.0005 m` (`0.5 mm`) to account for the algebraic linear inverse approximation in 2D Helmert while assuring sub-millimeter precision.
  - Angular tolerance: `1e-10 rad` (matching the 6th decimal of a sexagesimal second / `~2.78e-10 deg`).
  - Vertical tolerance: `1e-9 m` (floating-point epsilon, exact).
- [x] Update documentation in `docs/testing.md`.

---

## Implementation Details

### 1. No-Data Handling & Grid Boundary Analysis
- **Grid Extent & Mask Invariant**: The new `.spg` grid (`rom_grid3d_25.09.spg`) and old `.grd` grid (`ETRS89_KRASOVSCHI42_2DJ.GRD`) have identical bounding boxes and identical dimensions ($53 \times 72 = 3,816$ nodes), with exactly 1,127 unpopulated/no-data nodes at the exact same indices.
- **Bug Resolution**:
  In `pytransdatro/trans_grid.py`, `_grid_vs_no_data()` previously only inspected `if self.no_data in v:` (`self.no_data = 999`). Because SPG loads NaNs, it skipped this check and returned `(nan, nan)` without raising `NoDataGridErr`.
  Changed to:
  ```python
  def _grid_vs_no_data(self, values):
      for v in values:
          if any(math.isnan(x) or x == self.no_data for x in v):
              return True
      return False
  ```
  Checking `math.isnan(x)` first provides optimal execution speed on the default SPG grid.
- **Derivation Rules Documented**:
  - `test_st70_to_etrs89_no_grid_inner_no_data`: Coordinates placed 1 mm inside grid corner nodes $\pm 1 \times \text{step}$ and $\pm 2 \times \text{step}$. The 1 mm inward offset places them inside the coverage bounding box, but the required 4x4 interpolation subgrid includes unpopulated corner nodes, raising `NoDataGridErr`.
  - `test_st70_to_etrs89_just_out_ro_border_no_data`: Representative coordinates immediately across Romania's border (Hungary, Ukraine, Moldova, Bulgaria/Black Sea) that lie within the rectangular bounding box but outside Romania's populated interpolation cells.

### 2. Removal of Hardcoded P1–P9 Tests
- Removed `ST70_POINTS`, `ETRS89_POINTS`, `POINT_IDS`, `st70_pnts`, and `etrs89_pnts` fixtures.
- Removed the 36 parameterized discrete tests:
  - `test_st70_to_etrs89_2D[P1..P9]`
  - `test_st70_to_etrs89_3D[P1..P9]`
  - `test_etrs89_to_st70_2D[P1..P9]`
  - `test_etrs89_to_st70_3D[P1..P9]`

### 3. Round-Trip Tests Added
- **`test_st70_to_etrs89_to_st70_roundtrip`**:
  - Processes all 9,046 points from `tests/data/st70_to_etrs89_input.csv`.
  - Asserts $\Delta N \le 0.0005\text{ m}$, $\Delta E \le 0.0005\text{ m}$, $\Delta Z \le 10^{-9}\text{ m}$.
- **`test_etrs89_to_st70_to_etrs89_roundtrip`**:
  - Processes coordinates from `tests/data/etrs89_to_st70_input.csv`.
  - Correctly catches grid boundary exceptions (`NoDataGridErr`, `OutOfGridErr`) for points on the outermost border whose forward projection places them into margin cells lacking 4x4 interpolation nodes.
  - Asserts $\Delta\text{Lat} \le 10^{-10}\text{ rad}$, $\Delta\text{Lon} \le 10^{-10}\text{ rad}$, $\Delta h \le 10^{-9}\text{ m}$.
  - Confirms $>8,000$ points are validated.

---

## Test Verification Results

Ran `pytest tests/test_trans_ro.py -v`:
- Total tests: 42 test cases
- Passed: 40 passed, 2 xpassed (bulk file tests with release mismatch)
- 100% of active tests passing.
