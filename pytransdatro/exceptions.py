"""Defines domain exceptions that are raised by PyTransDatRO."""

from typing import Any


class PyTransDatError(Exception):
    """Base exception class for all PyTransDatRO errors."""

    pass


class OutOfGridErr(PyTransDatError):
    """Signifies that the input location (N,E) is not covered by a grid
    or is outside the grid's bounding envelope.
    """

    def __init__(self, n: float, e: float, grid: Any) -> None:
        """Initialize OutOfGridErr with out-of-bounds coordinates and grid metadata.

        :param n: Northing coordinate in meters or latitude in degrees
        :type n: float
        :param e: Easting coordinate in meters or longitude in degrees
        :type e: float
        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid
        """
        message = (
            f"({n}, {e}) location is out of the grid area defined by file {grid.source}. "
            f"Location must be inside the grid defined by "
            f"Min N = {grid.n_min}, Min E = {grid.e_min}, "
            f"Max N = {grid.n_max}, Max E  = {grid.e_max}."
        )
        super().__init__(message)


class NoDataGridErr(PyTransDatError):
    """Signifies that grid shift values required for interpolation at input location (N,E)
    contain No Data (NaN or 999) values, preventing computation.
    """

    def __init__(self, n: float, e: float, grid: Any) -> None:
        """Initialize NoDataGridErr with location and grid metadata.

        :param n: Northing coordinate in meters or latitude in degrees
        :type n: float
        :param e: Easting coordinate in meters or longitude in degrees
        :type e: float
        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid
        """
        message = (
            f"(The subgrid of size {grid.sg_size}x{grid.sg_size} nodes required "
            f"for the computation of interpolated value at {n}, {e} location "
            f"has No Data value(s)."
        )
        super().__init__(message)


class OutOfRangeIndexGridErr(PyTransDatError):
    """Signifies that the input index is out of the range of available grid indexes."""

    def __init__(self, idx: int, grid: Any) -> None:
        """Initialize OutOfRangeIndexGridErr with invalid index and grid metadata.

        :param idx: Linear node index queried
        :type idx: int
        :param grid: Grid against which index coordinates are queried
        :type grid: Grid
        """
        message = (
            f"(Grid {grid.file_name} has index values ranging from 0 to "
            f"{grid.c_count * grid.r_count - 1}. The requested index ({idx}) "
            f"is out of range."
        )
        super().__init__(message)


class MissingGridError(FileNotFoundError, PyTransDatError):
    """Signifies that no valid grid file was found in the grids directory."""

    pass


class AmbiguousGridError(ValueError, PyTransDatError):
    """Signifies that multiple grid files were found when exactly one was expected."""

    pass
