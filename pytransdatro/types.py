"""Type definitions and aliases for PyTransDatRO.

Provides standardized type aliases for 2D and 3D coordinate tuples,
polymorphic input sequences, and batch transformation structures.
"""

from typing import Sequence, Tuple, Union

# Scalar coordinate values
CoordinateScalar = float

# Sequence of coordinate values (lists, tuples)
CoordinateSequence = Sequence[float]

# 2D coordinate pair: (Northing, Easting) or (Latitude, Longitude) in radians
Coord2D = Tuple[float, float]

# 3D coordinate triplet: (Northing, Easting, Normal Elevation) or (Latitude, Longitude, Ellipsoidal Height)
Coord3D = Tuple[float, float, float]

# Polymorphic input allowing a single float coordinate or an iterable sequence of floats
ScalarOrSeq = Union[float, Sequence[float]]

# Batch transformation return types: parallel lists of coordinates
Sequence2D = Tuple[list[float], list[float]]
Sequence3D = Tuple[list[float], list[float], list[float]]
