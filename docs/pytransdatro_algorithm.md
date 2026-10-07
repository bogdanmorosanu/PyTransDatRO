# PyTransDatRO Transformation Algorithms

This document outlines the sequence of operations for 2D and 3D coordinate transformations in `pytransdatro` (found in `pytransdatro/trans_ro.py`).

## 1. Stereo70 to ETRS89 (`st70_to_etrs89`)

**Input:** `n` (Northing), `e` (Easting) in Stereo70, and optionally `z` (Elevation/Orthometric height) in MN75.

1.  **Grid Shift (Inverse):** 
    Applies a correction using the distortion grid to account for local geodetic anomalies.
    *   **Operation:** `r_n, r_e = _t_gr2d.trans(n, e, -1)`
    *   **Logic:** Uses bicubic interpolation (`Grid2D.interp`) to find shifts `shift_n`, `shift_e` at the provided `n, e`. Because direction is `-1`, it subtracts the shift: `r_n = n - shift_n`, `r_e = e - shift_e`.
    *   **Note:** The 2D grid shift is applied *first*.

2.  **Helmert 2D Transformation (Forward):**
    Applies the Helmert 7-parameter (simplified to 2D) conformal transformation.
    *   **Operation:** `r_n, r_e = _t_h2d.trans(r_n, r_e, 1)`
    *   **Logic:** Transforms from Dealul Piscului (Stereo70) system to the global ETRS89 system.

3.  **Stereographic Projection to Geodetic:**
    Converts from grid coordinates (Northing, Easting) to geographic coordinates (Latitude, Longitude).
    *   **Operation:** `lat, lon = _p_st70.to_geo(r_n, r_e)`
    *   **Logic:** Uses the oblique stereographic projection formulas for the WGS84 ellipsoid.

4.  **Elevation / 1D Grid Transformation (Optional):**
    If the `z` parameter is provided, a 1D grid interpolation calculates the geoid height (undulation).
    *   **Operation:** `r_z, = grid_1d.trans(lat_deg, lon_deg, z, 1)`
    *   **Logic:** Extracts the interpolated height shift (`h_shift`) based on the output ETRS89 `lat` and `lon` (converted to degrees). Because direction is `1`, it calculates ellipsoidal height as: `z_etrs89 = z + h_shift`.

**Output:** `lat` (Latitude), `lon` (Longitude), and optionally `z_etrs89` (Ellipsoidal Height) in ETRS89.

---

## 2. ETRS89 to Stereo70 (`etrs89_to_st70`)

**Input:** `lat` (Latitude), `lon` (Longitude) in ETRS89, and optionally `h` (Ellipsoidal Height).

1.  **Geodetic to Stereographic Projection:**
    Converts geographic coordinates (Latitude, Longitude) to grid coordinates (Northing, Easting) on the projection plane.
    *   **Operation:** `r_n, r_e = _p_st70.to_grid(lat, lon)`
    *   **Logic:** Uses oblique stereographic projection.

2.  **Helmert 2D Transformation (Inverse):**
    Reverses the Helmert transformation to convert back to the local Dealul Piscului system.
    *   **Operation:** `r_n, r_e = _t_h2d.trans(r_n, r_e, -1)`

3.  **Grid Shift (Forward):**
    Applies the distortion grid shift to obtain the final localized coordinates.
    *   **Operation:** `r_n, r_e = _t_gr2d.trans(r_n, r_e, 1)`
    *   **Logic:** Uses bicubic interpolation (`Grid2D.interp`) on the intermediate `r_n, r_e` to get `shift_n, shift_e`. Because direction is `1`, it adds the shift: `final_n = r_n + shift_n`, `final_e = r_e + shift_e`.
    *   **Note:** The 2D grid shift is applied *at the end*.

4.  **Elevation / 1D Grid Transformation (Optional):**
    If the `h` parameter is provided, a 1D grid interpolation calculates the orthometric height.
    *   **Operation:** `r_z, = grid_1d.trans(lat_deg, lon_deg, h, -1)`
    *   **Logic:** Extracts the interpolated height shift (`h_shift`) based on the *input* ETRS89 `lat` and `lon` (converted to degrees). Because direction is `-1`, it calculates orthometric height as: `H_stereo70 = h - h_shift`.

**Output:** `n` (Northing), `e` (Easting), and optionally `z` (Orthometric height MN75) in Stereo70.
