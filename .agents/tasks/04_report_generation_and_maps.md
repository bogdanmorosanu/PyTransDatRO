# Report Generation and Interactive Maps

## Objectives
- [x] Create a tool (`tools/generate_reports.py`) to process the bulk transformation test data and generate comprehensive precision reports for both ETRS89 $\rightarrow$ Stereo70 and Stereo70 $\rightarrow$ ETRS89.
- [x] Calculate precise coordinate differences between the computed and expected values (in meters/seconds) and assign precision flags (0: <0.005m, 1: 0.005-0.015m, 2: 0.015-0.05m, 4: >0.05m).
- [x] Generate CSV files containing detailed point-by-point differences and flags, designed for external software analysis (e.g., Excel, QGIS).
- [x] Generate plain-text summary reports outlining maximum/mean deviations and precision distributions across 2D (planimetric) and 1D (elevation) coordinates.
- [x] Generate interactive HTML maps to visualize the precision of the coordinate transformations across Romania.
- [x] Ensure the maps are strictly open-source, non-commercial, and fully functional locally (via `file:///` protocol) without needing a web server.

## Implementation Details
1. **Script creation**: Implemented `tools/generate_reports.py` which loads the test input/expected data directly and runs the transformations using the core `PyTransDatRO` logic.
2. **Precision Analysis**: Added separate logic for 2D precision (Easting/Northing or Lat/Lon) and 1D precision (Elevation) to reflect that the app uses two separate grid files for planar and vertical conversions.
3. **Artifact generation**: The script outputs files to the `reports/` folder (`.csv`, `.txt`, `.html`).
4. **Interactive Maps**: Built using `Leaflet` and `OpenTopoMap` raster tiles. This specific stack was chosen to avoid CORS / `fetch()` restrictions that block modern vector tiles (like MapLibre/OpenFreeMap) when opened directly from the filesystem.
5. **Map features**: Includes separate toggleable layers for 2D and 1D precisions (since the grids update independently) and detailed point popups indicating differences and flag values.

## Completion Status
The script runs successfully, generating accurate statistical reports and fully rendering local HTML maps with working interactive layers and open-source topographic backgrounds. The task is considered complete.
