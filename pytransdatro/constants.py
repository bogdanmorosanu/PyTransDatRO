"""Geodetic and mathematical constants for PyTransDatRO.

This module houses all physical, ellipsoid, projection, and geographic envelope
constants used across coordinate transformation pipelines. All values are pure
Python constants adhering to the standard library invariant.
"""

import math

# Angular conversion factors
DEG_TO_RAD_FACTOR: float = math.pi / 180.0
RAD_TO_DEG_FACTOR: float = 180.0 / math.pi

# WGS84 Reference Ellipsoid parameters
WGS84_ELL_SMAJ_AXIS: float = 6378137.0  # Semi-major axis 'a' (meters)
WGS84_ELL_INV_FLAT: float = 298.257223563  # Inverse flattening '1/f'

# Stereographic Oblique 1970 (Stereo 70) Projection parameters
# Natural origin in radians (sexagesimal equivalent: 46° 00' 00.0" N and 25° 00' 00.0" E)
STEREO70_ORIG_LAT: float = 46.0 * DEG_TO_RAD_FACTOR  # 46° 00' 00.0" N (~0.8028514559173916 rad)
STEREO70_ORIG_LON: float = 25.0 * DEG_TO_RAD_FACTOR  # 25° 00' 00.0" E (~0.4363323129985824 rad)
STEREO70_FALSE_N: float = 500000.0  # False northing (meters)
STEREO70_FALSE_E: float = 500000.0  # False easting (meters)
STEREO70_SCALE_FACTOR: float = 0.99975  # Scale factor at natural origin (k_0)

# Projection algorithm iteration safeguards
STEREO70_CONV_TOLERANCE: float = 0.0000000000484814  # Convergence threshold in radians
STEREO70_MAX_ITERATIONS: int = 50  # Iteration cap for conformal latitude resolution

# Romania Territorial Geographic Envelope (WGS84 Bounding Box)
ROMANIA_BBOX_LAT_MIN: float = 43.5  # Southern boundary limit (degrees)
ROMANIA_BBOX_LAT_MAX: float = 48.3  # Northern boundary limit (degrees)
ROMANIA_BBOX_LON_MIN: float = 20.2  # Western boundary limit (degrees)
ROMANIA_BBOX_LON_MAX: float = 29.8  # Eastern boundary limit (degrees)
