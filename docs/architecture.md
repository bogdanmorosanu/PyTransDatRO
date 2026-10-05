# PyTransDatRO Architecture

## Overview
`pytransdatro` provides bidirectional coordinate transformations between **Stereo70** (Romanian projected CRS) and **ETRS89** (geographic CRS). The library is designed to match the output of the official TransDatRO application byte-for-byte, including undocumented edge cases.

## Class Structure
- **`TransRO`** (`trans_ro.py`): The main coordinator class. It composes instances of grid, Helmert, and projection classes to perform the full mathematical pipeline.
- **`Grid` (Abstract)** (`trans_grid.py`): Base class for binary grid parsing, coordinate checking, and bicubic interpolation caching.
  - **`BiInterp`**: Nested class for evaluating the 16-coefficient bicubic spline polynomial on a 4x4 subgrid.
  - **`Grid1D`**: Concrete class for height/elevation corrections (Z-shift).
  - **`Grid2D`**: Concrete class for northing/easting corrections (N, E shifts).
- **`Helmert2D`** (`trans_helmert2d.py`): Encapsulates the 4-parameter 2D Helmert transformation constants (Translations, PPM/Scale, Rotation) between Stereo70 and StereoGRS80.
- **`StereoProj`** (`proj_stereo.py`): Implements Oblique Stereographic projection on the WGS84 ellipsoid.
  - **`_Ellipsoid`**: Handles basic ellipsoid element calculations.
  - **`_ConfSph`**: Defines the conformal sphere and computes conformal latitudes and longitudes.

## Pipeline Flow
### Stereo70 -> ETRS89 (`st70_to_etrs89`)
1. **2D Grid Transformation**: Applies `ETRS89_KRASOVSCHI42_2DJ.GRD` shifts (subtraction).
2. **2D Helmert Transformation**: Converts Stereo70 to StereoGRS80.
3. **Stereographic Projection**: Converts grid coordinates to geographic coordinates `(Lat, Long)`.
4. **Height Correction (Optional)**: If a `Z` coordinate is provided, a 1D grid is used for the height correction. (Defaults to `EGG97_QGRJ.GRD`, but switches to `zitaBucx.grd` in the Bucharest bounding box).

### ETRS89 -> Stereo70 (`etrs89_to_st70`)
1. **Stereographic Projection**: Transforms geographic coordinates to StereoGRS80 grid coordinates.
2. **2D Helmert Transformation (Inverse)**: Converts StereoGRS80 to Stereo70.
3. **2D Grid Transformation**: Applies the `ETRS89_KRASOVSCHI42_2DJ.GRD` shifts (addition).
4. **Height Correction (Optional)**: Subtracts height correction via the appropriate 1D grid.

## Current Status
The mathematical pipeline is complete. The system correctly parses binary grid files, caches lookups, and employs bicubic interpolation. Notably, the library correctly handles the undocumented Bucharest 1D grid fallback (`zitaBucx.grd`) mirroring the behavior of the official software.
