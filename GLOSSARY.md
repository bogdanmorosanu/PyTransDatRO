# PyTransDatRO

Bidirectional transformation of Romanian coordinates between the Stereo 70 national projection
and ETRS89, following the official ANCPI algorithm.

## Grid coverage

**Grid envelope**:
The rectangular bounding box of a transformation grid. A coordinate outside it cannot be
transformed at all.
_Avoid_: Grid extent, coverage area, outside grid (for no-data cases)

**No-data zone**:
An area inside the grid envelope where the grid holds no valid values, typically just beyond
the national border. A coordinate here is inside the envelope but still cannot be transformed.
_Avoid_: Out of grid, outside grid

**Untransformable point**:
A point that falls outside the grid envelope or in a no-data zone; its output carries a warning
naming which of the two applies, never a coordinate.
_Avoid_: Failed point, invalid point
