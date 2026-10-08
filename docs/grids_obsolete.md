> [!CAUTION]
> **OBSOLETE DOCUMENTATION**
> 
> This document describes the legacy binary `.GRD` files (`ETRS89_KRASOVSCHI42_2DJ.GRD`, `EGG97_QGRJ.GRD`, and `zitaBucx.grd`) used by previous versions of TransDatRO / pyTransDatRO.
> 
> As of the 2025.09 release, **pyTransDatRO uses a single unified 3D SPG grid file (`rom_grid3d_25.09.spg`)**.
> 
> For active grid documentation, please refer to [grid.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/grid.md).

---

# Legacy Grid Files Reference (Obsolete)

This document describes the legacy binary grid files (`.GRD` or `.grd`) previously utilized by `pytransdatro` for coordinate transformations between the Stereo70 and ETRS89 coordinate reference systems.

## Legacy Grid Files & Their Purpose

1. **`ETRS89_KRASOVSCHI42_2DJ.GRD`**
   - **Type**: 2D Grid
   - **Purpose**: Provided horizontal coordinate shifts (Northing and Easting). It was used during both the forward and reverse transformations between ETRS89 and Stereo70.

2. **`EGG97_QGRJ.GRD`**
   - **Type**: 1D Grid
   - **Purpose**: Provided vertical (elevation/height) corrections based on the European Gravimetric (Quasi)Geoid 1997.

3. **`zitaBucx.grd`**
   - **Type**: 1D Grid
   - **Purpose**: Localized grid providing elevation corrections strictly for the Bucharest metropolitan area.

---

## Legacy Binary Grid Format Specification

The legacy grid files were written in a strictly structured binary format with byte offsets.

### 1. Header Section (Bytes 0 - 47)
Every grid file began with a 48-byte header consisting of six consecutive 64-bit (8-byte) double-precision floating-point numbers stored in Little-Endian format (`<d` in Python's `struct`).

| Byte Range | Type   | Name     | Description |
|------------|--------|----------|-------------|
| `0 - 7`    | Double | `e_min`  | The minimum Easting (or Longitude) boundary of the grid. |
| `8 - 15`   | Double | `e_max`  | The maximum Easting (or Longitude) boundary of the grid. |
| `16 - 23`  | Double | `n_min`  | The minimum Northing (or Latitude) boundary of the grid. |
| `24 - 31`  | Double | `n_max`  | The maximum Northing (or Latitude) boundary of the grid. |
| `32 - 39`  | Double | `e_step` | The grid resolution (distance between nodes) along the Easting axis. |
| `40 - 47`  | Double | `n_step` | The grid resolution (distance between nodes) along the Northing axis. |

*Dimensions of the grid were calculated as:*
- **Columns** (`c_count`): `round((e_max - e_min) / e_step) + 1`
- **Rows** (`r_count`): `round((n_max - n_min) / n_step) + 1`

### 2. Data Section (Byte 48 Onwards)
Immediately following the header, the file contained the raw grid node values ordered row by row (South to North), moving left to right (West to East):

- **1D Grids** (`v_size = 1`): Each node contains 1 Double (8 bytes) representing the elevation shift.
- **2D Grids** (`v_size = 2`): Each node contains 2 Doubles (16 bytes) representing Easting and Northing shifts.

### 3. No Data Values
Legacy grids used `999.0` as a sentinel "No Data" value.
