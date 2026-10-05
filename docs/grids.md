# Grid Files Reference

This document describes the binary grid files (`.GRD` or `.grd`) utilized by `pytransdatro` for coordinate transformations between the Stereo70 and ETRS89 coordinate reference systems. These grids are supplied by the official ROMPOS/TransDatRO distribution.

## Grid Files & Their Purpose

1. **`ETRS89_KRASOVSCHI42_2DJ.GRD`**
   - **Type**: 2D Grid
   - **Purpose**: Provides horizontal coordinate shifts (Northing and Easting). It is used during both the forward and reverse transformations between ETRS89 and Stereo70 (specifically acting as the distortion surface over the Krasovsky 1942 basis).

2. **`EGG97_QGRJ.GRD`**
   - **Type**: 1D Grid
   - **Purpose**: Provides vertical (elevation/height) corrections based on the European Gravimetric (Quasi)Geoid 1997. It is applied when transforming a 3D coordinate (i.e., when a Z/Height value is provided).

3. **`zitaBucx.grd`**
   - **Type**: 1D Grid
   - **Purpose**: A highly specific, localized grid providing elevation corrections strictly for the Bucharest metropolitan area. The library switches to this grid automatically when coordinates fall within the Bucharest bounding box, replicating undocumented TransDatRO behavior.

---

## Binary Grid Format Specification

The grid files are written in a strictly structured binary format. There is no string-based metadata or XML; everything is parsed mathematically based on byte offsets. 

### 1. Header Section (Bytes 0 - 47)
Every grid file begins with a 48-byte header consisting of six consecutive 64-bit (8-byte) double-precision floating-point numbers. They are stored in Little-Endian format (`<d` in Python's `struct`).

| Byte Range | Type   | Name     | Description |
|------------|--------|----------|-------------|
| `0 - 7`    | Double | `e_min`  | The minimum Easting (or Longitude) boundary of the grid. |
| `8 - 15`   | Double | `e_max`  | The maximum Easting (or Longitude) boundary of the grid. |
| `16 - 23`  | Double | `n_min`  | The minimum Northing (or Latitude) boundary of the grid. |
| `24 - 31`  | Double | `n_max`  | The maximum Northing (or Latitude) boundary of the grid. |
| `32 - 39`  | Double | `e_step` | The grid resolution (distance between nodes) along the Easting axis. |
| `40 - 47`  | Double | `n_step` | The grid resolution (distance between nodes) along the Northing axis. |

*From this header, the dimensions of the grid are calculated as:*
- **Columns** (`c_count`): `round((e_max - e_min) / e_step) + 1`
- **Rows** (`r_count`): `round((n_max - n_min) / n_step) + 1`

### 2. Data Section (Byte 48 Onwards)
Immediately following the header, the file contains the raw grid node values. The nodes are ordered row by row (from South to North), moving left to right (West to East) across columns.

- **1D Grids** (`v_size = 1`):
  Each node contains **1 Double** (8 bytes) representing the elevation shift.
- **2D Grids** (`v_size = 2`):
  Each node contains **2 Doubles** (16 bytes) representing the Easting and Northing shifts respectively.

#### Locating Data
To find the grid values for a specific node at Row `r` and Column `c`:
1. **Calculate the flat index**: `idx = r * c_count + c`
2. **Calculate the byte offset**: `offset = 48 + (idx * v_size * 8)`
3. **Read the values**: Read the next `v_size * 8` bytes.

### 3. No Data Values
To mark areas outside the Romanian border (but still inside the rectangular grid bounding box), the grids use **`999.0`** as a sentinel "No Data" value. If any of the 16 nodes required for a bicubic spline interpolation contains `999.0`, the library raises a `NoDataGridErr`.
