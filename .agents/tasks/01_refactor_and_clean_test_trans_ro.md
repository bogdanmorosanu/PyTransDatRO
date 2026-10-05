# Task: Refactor and Clean test_trans_ro.py

**Status**: Completed  
**Target File**: [tests/test_trans_ro.py](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/tests/test_trans_ro.py)  
**Created**: 2026-10-05  

---

## Objective

Modernize, clean, and fix identified minor issues in `tests/test_trans_ro.py`, reorganize test coordinate data files into a dedicated directory with a consistent naming convention, and create placeholder files for the reverse transformation without altering core transformation logic.

---

## Detailed Scope of Work

### 1. Test Data Reorganization & File Naming
- **Directory**: Create dedicated directory `tests/data/`.
- **Naming Convention**: Pair files by transformation direction and role (`_input.csv` and `_expected.csv`):
  - **Stereo70 -> ETRS89**:
    - `tests/data/st70_to_etrs89_input.csv` (moved/renamed from `tests/coos_for_testing_st70.txt`)
    - `tests/data/st70_to_etrs89_expected.csv` (moved/renamed from `tests/coos_for_testing_etrs89.txt`)
  - **ETRS89 -> Stereo70** (Placeholders for future test implementation):
    - `tests/data/etrs89_to_st70_input.csv` (empty file)
    - `tests/data/etrs89_to_st70_expected.csv` (empty file)

### 2. Replace Deprecated `pkg_resources` with `pathlib.Path`
- **Issue**: `pkg_resources` is deprecated and produces runtime warnings on modern Python.
- **Resolution**:
  - Remove `import pkg_resources`.
  - Import `from pathlib import Path`.
  - Update `st70_coo_file` and `etrs89_coo_file` fixtures to resolve file paths relative to `Path(__file__).parent / "data"`.

### 3. Fix Vertical Elevation Tolerance Assertion
- **Issue**: In `test_etrs89_to_st70_3D` (line 254), the elevation coordinate ($Z$) comparison uses `coo_plan_tol` instead of `coo_elev_tol`.
- **Resolution**: Update the comparison to `math.isclose(sut[point_id][2], st70[2], abs_tol=coo_elev_tol)`.

### 4. Remove Debug `print()` Calls
- **Issue**: `test_st70_to_etrs89_2D_fromfile` contains leftover debugging statements (`print(sut[i][0], coo_etrs89[0])` and `print(sut[i][1], coo_etrs89[1])`).
- **Resolution**: Remove the print statements to keep pytest output clean.

### 5. Parametrize Discrete Point Tests (P1–P9)
- **Issue**: `test_st70_to_etrs89_2D`, `test_st70_to_etrs89_3D`, `test_etrs89_to_st70_2D`, and `test_etrs89_to_st70_3D` iterate over points using `for` loops. A failure on point P1 stops the test immediately.
- **Resolution**:
  - Parametrize the tests across points `P1` through `P9` using `@pytest.mark.parametrize`.
  - Ensure individual point test IDs appear in pytest reporting (e.g., `test_st70_to_etrs89_2D[P1]`).
  - Maintain coordinate fixtures and docstrings intact.

---

## Acceptance Gates

1. **Gate 1**: `tests/data/` contains all 4 coordinated files (`st70_to_etrs89_input.csv`, `st70_to_etrs89_expected.csv`, and empty placeholders `etrs89_to_st70_input.csv`, `etrs89_to_st70_expected.csv`).
2. **Gate 2**: Zero deprecated `pkg_resources` references in `tests/`. File paths resolved via `pathlib.Path(__file__).parent / "data"`.
3. **Gate 3**: No residual `print()` statements in `tests/test_trans_ro.py`.
4. **Gate 4**: Discrete point tests execute as independent parametrized pytest cases for each point (P1–P9).
5. **Gate 5**: `test_etrs89_to_st70_3D` asserts elevation with `coo_elev_tol`.
6. **Gate 6**: All existing test assertions pass without regression.
