# SPG Grid Coordinates & Shifts CSV Extraction Tool

## Objectives
- [x] Investigate whether 2D planimetric shifts ($dN, dE$) and 1D elevation shifts ($Z / N$) in the unified Romgeo SPG grid file (`romgeo_grid.spg`) are defined at identical spatial locations.
- [x] Determine the appropriate CSV export strategy (single unified file vs. separate files) based on grid structures.
- [x] Implement a CLI tool in `tools/extract_spg_grid_csv.py` to extract grid coordinates and shifts to CSV format for webapp raster basemap generation.
- [x] Generate production extraction datasets and diagnostic report in `reports/`.
- [x] Document native grid properties vs. computed geodetic cross-reference coordinates.

## Architecture & Findings

### 1. Spatial Structure Investigation
An analysis of the binary grid datasets (`geodetic_shifts` and `geoid_heights`) revealed that they reside in fundamentally different coordinate systems and spatial domains:

- **2D Planimetric Distortion Grid (`geodetic_shifts`)**:
  - **Native CRS**: Stereo 70 Projected Coordinates (Metric, meters).
  - **Bounding Extents**: $N \in [213634.564, 785634.564]\,\text{m}$, $E \in [109783.040, 890783.040]\,\text{m}$.
  - **Grid Resolution / Dimensions**: $11,000\,\text{m} \times 11,000\,\text{m}$ step ($53 \times 72$ matrix, 3,816 total nodes).
  - **Active Coverage**: 2,689 valid land nodes (1,127 offshore/border nodes set to `NaN`).
  - **Export Artifact**: `reports/spg_grid_shifts_2d.csv`.

- **1D Elevation / Quasigeoid Undulation Grid (`geoid_heights`)**:
  - **Native CRS**: ETRS89 Ellipsoidal Geographic Coordinates (degrees).
  - **Bounding Extents**: $\text{Lat} \in [43.533333^\circ, 48.400000^\circ]$, $\text{Lon} \in [20.066667^\circ, 29.800000^\circ]$.
  - **Grid Resolution / Dimensions**: $0.033333^\circ \times 0.033333^\circ$ step ($2' \approx 3.7\,\text{km} \times 2.6\,\text{km}$, $147 \times 293$ matrix, 43,071 total nodes).
  - **Active Coverage**: 43,071 valid nodes (100% coverage, 0 NaNs).
  - **Export Artifact**: `reports/spg_grid_shifts_z.csv`.

**Conclusion**: Because the grid points do not coincide spatially, geometrically, or by CRS, two dedicated CSV datasets are required for accurate raster layer construction.

### 2. Extracted vs. Computed Data
- **Native / Extracted Parameters**:
  - `stereo70_n`, `stereo70_e`, `shift_dn_m`, `shift_de_m`, and validity flags in `spg_grid_shifts_2d.csv` are derived directly from the binary grid array and affine metadata.
  - `etrs89_lat_deg`, `etrs89_lon_deg`, `undulation_z_m`, and validity flags in `spg_grid_shifts_z.csv` are derived directly from the geoid grid array and affine metadata.
- **Computed Cross-Reference Coordinates**:
  - `etrs89_lat_deg, etrs89_lon_deg` in `spg_grid_shifts_2d.csv`: Computed via `pytransdatro.trans_ro.TransRO` (reversing 2D distortion, inverse Oblique Stereographic projection, and 7-parameter Helmert transformation to ETRS89) to allow direct ingestion into web map projection pipelines (EPSG:4326/3857).
  - `stereo70_n, stereo70_e` in `spg_grid_shifts_z.csv`: Computed via `TransRO._p_st70` (forward Oblique Stereographic projection) to allow direct rendering in Romanian national grid coordinates (EPSG:31700).

## Files Created & Updated

1. **`tools/extract_spg_grid_csv.py`**:
   - CLI script supporting extraction of 2D and 1D grids.
   - CLI parameters: `--all-nodes` (to include invalid/NaN nodes) and `--output-dir` (custom output location, default `reports/`).
   - Uses `pytransdatro` core coordinate transformer while strictly respecting pure Python invariants.
2. **`reports/spg_grid_shifts_2d.csv`**:
   - 2,689 valid data rows with columns: `row`, `col`, `node_index`, `stereo70_n`, `stereo70_e`, `shift_dn_m`, `shift_de_m`, `etrs89_lat_deg`, `etrs89_lon_deg`, `is_valid`.
3. **`reports/spg_grid_shifts_z.csv`**:
   - 43,071 valid data rows with columns: `row`, `col`, `node_index`, `etrs89_lat_deg`, `etrs89_lon_deg`, `undulation_z_m`, `stereo70_n`, `stereo70_e`, `is_valid`.
4. **`reports/spg_grid_extraction_summary.txt`**:
   - Text report capturing grid matrix dimensions, step sizes, coordinate boundaries, node counts, and CRS specifications.

## Verification & Testing
- **Script Execution**: Executed `tools/extract_spg_grid_csv.py` via `D:\anaconda3\envs\pytransdat\python.exe`.
- **Output Inspection**:
  - Verified row count, schema, and headers for `spg_grid_shifts_2d.csv` (2,689 valid nodes).
  - Verified row count, schema, and headers for `spg_grid_shifts_z.csv` (43,071 valid nodes).
  - Verified summary output in `reports/spg_grid_extraction_summary.txt`.
- **Core Test Suite**:
  - Executed `pytest tests/test_trans_ro.py`: 46 passed.
  - Executed full test suite `pytest`: 71 passed, 0 failed in 5.43s.
