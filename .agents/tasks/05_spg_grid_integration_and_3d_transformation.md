# SPG Grid Integration and 3D Elevation Transformation

## Objectives
- [x] Investigate and understand the new merged 3D grid format (`rom_grid3d_25.09.spg`) from the official Romgeo application.
- [x] Permit `numpy` strictly for loading and deserializing binary grid files while maintaining pure Python standard library for all coordinate processing algorithms.
- [x] Document the Romgeo transformation architecture and elevation algorithm in `docs/romgeo/`.
- [x] Implement `pytransdatro/spg_reader.py` (`SpgReader`) to read, parse, and unpickle `.spg` grids.
- [x] Integrate 2D grid shifts into `Grid2D` in `pytransdatro/trans_grid.py` using bicubic spline interpolation (`BiInterp`).
- [x] Diagnose the 3D elevation discrepancy between `pyTransDatRO` and `Romgeo`.
- [x] Identify grid metadata specifying nearest-neighbor colocation (`INTERP_COLOCATE = 0`) for vertical shifts vs. bicubic spline (`INTERP_BICUBIC = 2`) for horizontal shifts.
- [x] Update `Grid1D.trans()` in `pytransdatro/trans_grid.py` to use nearest-neighbor interpolation for `.spg` vertical grid data.
- [x] Update test expectations in `tests/test_trans_ro.py` for benchmark points P1–P9 to match Romgeo official coordinates and verify all 2D and 3D tests pass.

## Implementation Details

1. **Policy & Dependencies Update**:
   - Updated `GEMINI.md` to allow `numpy` strictly for reading and unpickling `.spg` binary grid files, ensuring all downstream math and interpolation routines remain pure Python.

2. **Algorithm & Architecture Research**:
   - Compared `romgeo_lite` transformation routines with `pytransdatro`.
   - Documented the flow in `docs/romgeo/romgeo_architecture.md` and `docs/romgeo/romgeo_elevation_algorithm.md`.
   - Verified the transformation pipeline sequence:
     - **Stereo70 $\rightarrow$ ETRS89**: 2D inverse projection $\rightarrow$ Helmert (Krassovsky to GRS80) $\rightarrow$ 2D Bicubic Grid Correction $\rightarrow$ Vertical Colocate Grid Correction ($h = H + N$).
     - **ETRS89 $\rightarrow$ Stereo70**: Vertical Colocate Grid Correction ($H = h - N$) $\rightarrow$ 2D Inverse Bicubic Grid Correction $\rightarrow$ Helmert (GRS80 to Krassovsky) $\rightarrow$ 2D Forward Projection.

3. **SPG Reader Implementation (`pytransdatro/spg_reader.py`)**:
   - Implemented `SpgReader` class using `numpy.load(..., allow_pickle=True)`.
   - Extracted and structured horizontal grid data (`shift_n_flat`, `shift_e_flat`, bounds, steps) and vertical geoid undulation data (`height_flat`, bounds, steps).

4. **Planimetric (2D) Integration**:
   - Updated `Grid` and `Grid2D` in `pytransdatro/trans_grid.py` to interface with `SpgReader` when configured with a `.spg` file.
   - Configured `pytransdatro/trans_ro.py` to use `rom_grid3d_25.09.spg` as default.
   - Validated that 2D interpolation matches expected coordinates across points P1–P9.

5. **Vertical (3D) Elevation Resolution**:
   - Discovered that SPG grid metadata explicitly sets `params.interpolation['vertical'] = 0` (`INTERP_COLOCATE`), meaning the official application uses nearest-neighbor lookup rather than bicubic spline interpolation for geoid height shifts.
   - Refactored `Grid1D.trans()` in `pytransdatro/trans_grid.py` to perform nearest-neighbor node lookup when `.spg` grids are active.

6. **Validation & Testing**:
   - Updated benchmark point coordinates (P1 through P9) in `tests/test_trans_ro.py` to align with the Romgeo official output.
   - Enabled and executed parametrized tests:
     - `test_st70_to_etrs89_2D[P1..P9]`
     - `test_etrs89_to_st70_2D[P1..P9]`
     - `test_st70_to_etrs89_3D[P1..P9]`
     - `test_etrs89_to_st70_3D[P1..P9]`
   - Confirmed 100% pass rate for discrete point transformation tests in both 2D and 3D directions.

## Completion Status
The integration of the `.spg` 3D grid format is complete and fully functional. Both planimetric (2D) and vertical (3D) transformations match the official Romgeo application outputs.
