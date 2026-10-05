# Testing in PyTransDatRO

This document summarizes the testing structure, datasets, tolerances, and test execution patterns implemented in `PyTransDatRO` (primarily in `tests/test_trans_ro.py` and the `tests/` directory).

---

## 1. Test Data & Fixtures

The test suite relies on official reference coordinates and verification thresholds based on the official ROMPOS documentation (*Help_TransDatRO_code_source_EN.pdf*) and outputs from the official `TransDatRO v4.08` software.

### A. Reference Points

The test suite defines dictionaries of reference points representing forward and inverse transformations:

- **`st70_pnts` (Stereo70 coordinates)**:
  - **Points P1 to P7**: Official reference points extracted from *Help_TransDatRO_code_source_EN.pdf* (page 3), defined with Northing ($N$), Easting ($E$), and Elevation ($Z$).
  - **Point P8**: Located inside the local 1D quasigeoid grid of Bucharest (`zitaBucx.grd`).
  - **Point P9**: Located immediately outside the Bucharest local grid area, falling back to the national 1D grid (`EGG97_QGRJ.GRD`).
- **`etrs89_pnts` (ETRS89 coordinates)**:
  - Official geographic coordinates in radians (`Latitude`, `Longitude`) and ellipsoidal height $h$ (meters), computed using `TransDatRO v4.08`.
  - Documented sexagesimal DMS equivalents are included in the fixture docstrings.

### B. Verification Tolerances

Tolerances define the maximum acceptable discrepancy against official reference outputs, matching the criteria specified in *Help_TransDatRO_code_source_EN.pdf* (page 3):

| Tolerance Fixture | Value | Physical Meaning | Description |
|-------------------|-------|------------------|-------------|
| `coo_rad_tol` | `0.00000000014544410` rad | $\approx 0.00003''$ | Angular tolerance for geographic latitude and longitude. |
| `coo_plan_tol` | `0.003` m | $3\text{ mm}$ | Planar tolerance for Stereo70 Northing and Easting coordinates. |
| `coo_elev_tol` | `0.003` m | $3\text{ mm}$ | Vertical tolerance for elevation / height components ($Z$ and $h$). |

### C. External Benchmark Files

For extensive surface verification across Romania, benchmark datasets are organized in a dedicated directory under `tests/data/` using a standardized naming convention (`<source>_to_<target>_<role>.csv`):

- **Stereo70 to ETRS89 Transformation**:
  1. **`tests/data/st70_to_etrs89_input.csv`**:
     - Contains 7,866 Stereo70 coordinates formatted as CSV (`id,n,e,z`).
     - Generated across a half-grid (nodes and mid-points) spanning the entire 2D shift grid extent.
  2. **`tests/data/st70_to_etrs89_expected.csv`**:
     - Contains the corresponding 7,866 ETRS89 coordinates processed through official `TransDatRO v4.08`.
     - Encoded in ANSI/Windows-1252, with coordinates in sexagesimal DMS strings (e.g., `47°42'56.40000"N, 22°28'31.99998"E`).
- **ETRS89 to Stereo70 Transformation (Placeholders)**:
  3. **`tests/data/etrs89_to_st70_input.csv`**: Placeholder file for reverse transformation inputs.
  4. **`tests/data/etrs89_to_st70_expected.csv`**: Placeholder file for reverse transformation reference outputs.

File paths are resolved dynamically using `pathlib.Path` through pytest fixtures (`st70_coo_file`, `etrs89_coo_file`).

---

## 2. Test Execution Patterns

All tests are implemented using `pytest` and execute against the top-level coordinator class `pytransdatro.TransRO`.

### A. Discrete Coordinate Transformations (Parametrized)

Discrete reference point tests (points `P1` through `P9`) use `@pytest.mark.parametrize("point_id", POINT_IDS)` where `POINT_IDS = [f"P{i}" for i in range(1, 10)]`. This ensures each point runs as an independent test case, providing isolated failure diagnostics:

- **2D Transformations**:
  - `test_st70_to_etrs89_2D`: Transforms $(N, E) \rightarrow (\text{lat}, \text{lon})$ for each point. Compares calculated angles against `etrs89_pnts` using `math.isclose(..., abs_tol=coo_rad_tol)`.
  - `test_etrs89_to_st70_2D`: Transforms $(\text{lat}, \text{lon}) \rightarrow (N, E)$ for each point. Compares calculated planar coordinates against `st70_pnts` using `math.isclose(..., abs_tol=coo_plan_tol)`.
- **3D Transformations (Elevation & Quasigeoid Shifts)**:
  - `test_st70_to_etrs89_3D`: Transforms $(N, E, Z) \rightarrow (\text{lat}, \text{lon}, h)$ for each point. Asserts angular agreement within `coo_rad_tol` and vertical ellipsoidal height within `coo_elev_tol`.
  - `test_etrs89_to_st70_3D`: Transforms $(\text{lat}, \text{lon}, h) \rightarrow (N, E, Z)$ for each point. Asserts planar coordinates within `coo_plan_tol` and normal elevation $Z$ within `coo_elev_tol`.

### B. Bulk Grid Validation (`test_st70_to_etrs89_2D_fromfile`)

- Ingests the 7,866 points from `tests/data/st70_to_etrs89_input.csv` and `tests/data/st70_to_etrs89_expected.csv`.
- Converts DMS strings to radians using `pytransdatro.utils.sexa_dms_chars_to_rad()`.
- Executes `st70_to_etrs89(n, e)` on all points.
- Validates that every point across the grid matches the official `TransDatRO` output within `coo_rad_tol`.

### C. Boundary & Exception Handling (Parametrized)

Edge cases and out-of-boundary behavior are tested using `@pytest.mark.parametrize`:

- **Out of Grid Boundaries (`test_st70_to_etrs89_extent_inner_out_of_grid`)**:
  - Evaluates a 16-point Cartesian product ($4 \times 4$ of $N$ and $E$) placed at or outside the outer grid bounding box.
  - Asserts that `pytransdatro.exceptions.OutOfGridErr` is raised.
- **No-Data Zones Within Grid Extents (`test_st70_to_etrs89_no_grid_inner_no_data`)**:
  - Evaluates 16 points located inside the grid bounding box, but positioned where the binary grid cells contain no valid interpolation data.
  - Asserts that `pytransdatro.exceptions.NoDataGridErr` is raised.
- **National Border Transition Areas (`test_st70_to_etrs89_just_out_ro_border_no_data`)**:
  - Tests 4 coordinates located immediately outside the Romanian territorial boundary.
  - Asserts that `pytransdatro.exceptions.NoDataGridErr` is raised.

### D. Test Data Generation Helper

- **`create_transdatro_test_file(file, template='st70_half_grid')`**:
  - A utility function contained within `test_trans_ro.py` (not a test itself).
  - Iterates through the active 2D grid (`_t_gr2d`), generating points at grid nodes and half-step intervals.
  - Outputs a CSV file that can be loaded directly into the official TransDatRO desktop software to produce new benchmark files.
