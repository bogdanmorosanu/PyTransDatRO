# PyTransDatRO Architecture

## Overview
`pytransdatro` provides bidirectional coordinate transformations between **Stereo70** (Romanian projected CRS) and **ETRS89** (geographic CRS). The library is designed to match the output of the official CNC / Romgeo application byte-for-byte, supporting unified 2D/3D transformations via modern `.spg` grid packages.

## Class Structure
- **`TransRO`** (`trans_ro.py`): The main coordinator class. It composes instances of grid, Helmert, and projection classes to perform the full mathematical pipeline.
- **`SpgReader`** (`spg_reader.py`): Singleton reader that auto-discovers and deserializes the active `.spg` grid file, caching flat Python coordinate arrays and metadata parameters.
- **`Grid` (Abstract)** (`trans_grid.py`): Base class for grid coordinate validation and interpolation caching.
  - **`BiInterp`**: Nested class for evaluating the 16-coefficient bicubic spline polynomial on a 4x4 subgrid.
  - **`Grid1D`**: Concrete class for height/elevation corrections (Z-shift / quasigeoid undulation) using nearest-neighbor colocate indexing.
  - **`Grid2D`**: Concrete class for northing/easting planimetric distortion corrections (N, E shifts) using bicubic spline interpolation.
- **`Helmert2D`** (`trans_helmert2d.py`): Encapsulates the 4-parameter 2D Helmert transformation constants (Translations, PPM/Scale, Rotation) between Stereo70 and StereoGRS80, loaded dynamically from the `.spg` metadata.
- **`StereoProj`** (`proj_stereo.py`): Implements Oblique Stereographic projection on the WGS84 ellipsoid.
  - **`_Ellipsoid`**: Handles basic ellipsoid element calculations.
  - **`_ConfSph`**: Defines the conformal sphere and computes conformal latitudes and longitudes.

## Pipeline Flow
### Stereo70 -> ETRS89 (`st70_to_etrs89`)
1. **2D Grid Transformation**: Applies `geodetic_shifts` planar distortion shifts from the SPG grid (subtraction).
2. **2D Helmert Transformation**: Converts Stereo70 to StereoGRS80 using dynamic SPG parameters.
3. **Stereographic Projection**: Converts grid coordinates to geographic coordinates `(Lat, Long)`.
4. **Height Correction (Optional)**: If a `Z` coordinate is provided, the 1D geoid layer (`geoid_heights`) provides the quasigeoid undulation $N$ via nearest-neighbor lookup: $h = H + N$.

### ETRS89 -> Stereo70 (`etrs89_to_st70`)
1. **Stereographic Projection**: Transforms geographic coordinates to StereoGRS80 grid coordinates.
2. **2D Helmert Transformation (Inverse)**: Converts StereoGRS80 to Stereo70.
3. **2D Grid Transformation**: Applies `geodetic_shifts` planar distortion shifts (addition).
4. **Height Correction (Optional)**: Subtracts height correction via the 1D geoid layer: $H = h - N$.

## Grid Management
The system automatically discovers and loads the single `.spg` grid located in `pytransdatro/grids/`. If multiple grids or no grids are found, explicit descriptive exceptions (`AmbiguousGridError`, `MissingGridError`) are raised to prevent transformation ambiguity. Legacy binary `.grd` files have been retired.
