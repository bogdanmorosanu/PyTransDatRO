"""Stereo70 and ETRS89 Geodetic Transformation Coordinator.

This module provides the TransRO class to coordinate 2D and 3D coordinate transformations
between Romania's official national cartographic projection (Stereo 70) and the
European Terrestrial Reference System 1989 (ETRS89).

Transformation Pipelines:
    - Forward (Stereo 70 -> ETRS89):
        1. 2D grid inverse shift correction (:class:`~pytransdatro.trans_grid.Grid2D`)
        2. 2D Helmert conformal transformation (:class:`~pytransdatro.trans_helmert2d.Helmert2D`)
        3. Oblique Stereographic projection to ellipsoidal coordinates (:class:`~pytransdatro.proj_stereo.StereoProj`)
        4. Optional 1D geoid height correction: :math:`h = H + N` (:class:`~pytransdatro.trans_grid.Grid1D`)

    - Inverse (ETRS89 -> Stereo 70):
        1. Oblique Stereographic conversion to projected coordinates (:class:`~pytransdatro.proj_stereo.StereoProj`)
        2. 2D Helmert inverse transformation (:class:`~pytransdatro.trans_helmert2d.Helmert2D`)
        3. 2D grid forward shift correction (:class:`~pytransdatro.trans_grid.Grid2D`)
        4. Optional 1D geoid height correction: :math:`H = h - N` (:class:`~pytransdatro.trans_grid.Grid1D`)
"""

import logging
import math
from typing import Any, Optional, Sequence, Union, overload

from pytransdatro.proj_stereo import StereoProj
from pytransdatro.trans_grid import Grid1D, Grid2D
from pytransdatro.trans_helmert2d import Helmert2D
from pytransdatro.types import Coord2D, Coord3D, ScalarOrSeq, Sequence2D, Sequence3D

logger = logging.getLogger(__name__)


class TransRO:
    """Coordinates geodetic transformations between Stereo 70 and ETRS89 coordinate systems."""

    def __init__(
        self,
        grid_filename: Optional[str] = None,
        telemetry_logger: Optional[Any] = None,
    ) -> None:
        """Initialize transformation engines with SPG binary grid and optional telemetry logger.

        :param grid_filename: Optional filename or path of the .spg grid file.
            If None, the grid is auto-discovered from pytransdatro/grids/.
        :type grid_filename: Optional[str]
        :param telemetry_logger: Optional asynchronous telemetry collector instance
        :type telemetry_logger: Optional[Any]
        """
        self._t_gr2d: Grid2D = Grid2D(grid_filename)
        self._t_h2d: Helmert2D = Helmert2D()
        self._p_st70: StereoProj = StereoProj()
        self._t_gr1d: Grid1D = Grid1D(grid_filename)
        self._collector: Optional[Any] = telemetry_logger

    @overload
    def st70_to_etrs89(
        self,
        n: float,
        e: float,
        z: None = None,
        source: int = 0,
    ) -> Coord2D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: float,
        e: float,
        z: float,
        source: int = 0,
    ) -> Coord3D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: None = None,
        source: int = 0,
    ) -> Sequence2D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: Sequence[float],
        source: int = 0,
    ) -> Sequence3D: ...

    def st70_to_etrs89(
        self,
        n: ScalarOrSeq,
        e: ScalarOrSeq,
        z: Optional[ScalarOrSeq] = None,
        source: int = 0,
    ) -> Union[Coord2D, Coord3D, Sequence2D, Sequence3D]:
        """Transforms Stereo70 grid coordinates to ETRS89 geographic coordinates.

        Supports both single coordinates (floats) and coordinate sequences (lists/tuples).
        :math:`(N, E[, H]) \\to (\\text{Lat}, \\text{Long}[, h])`

        :param n: Northing coordinate in meters (or sequence of northings)
        :type n: float or Sequence[float]
        :param e: Easting coordinate in meters (or sequence of eastings)
        :type e: float or Sequence[float]
        :param z: Normal elevation in meters (or sequence of elevations), optional
        :type z: float or Sequence[float], optional
        :param source: Integer source identifier (e.g. 0=unknown, 1=web_map, 2=rest_api)
        :type source: int, optional
        :return: (lat, lon[, h]) in ETRS89 coordinate reference system (radians, meters)
        :rtype: Tuple[float, float] or Tuple[float, float, float] or Sequence2D or Sequence3D
        :raises OutOfGridErr: If (n, e) falls outside the active grid boundaries
        :raises NoDataGridErr: If subgrid interpolation encounters No Data (NaN or 999) values
        """
        is_scalar = not isinstance(n, (list, tuple))
        n_seq: Sequence[float]
        e_seq: Sequence[float]
        z_seq: Optional[Sequence[float]]

        if is_scalar:
            n_seq, e_seq = (float(n),), (float(e),)  # type: ignore[arg-type]
            z_seq = (float(z),) if z is not None else None  # type: ignore[arg-type]
        else:
            n_seq, e_seq = n, e  # type: ignore[assignment]
            z_seq = z  # type: ignore[assignment]

        r_n, r_e = self._t_gr2d.trans(n_seq, e_seq, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, 1)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)

        if self._collector is not None:
            try:
                is_3d = 1 if z is not None else 0
                self._collector.record(1, r_n, r_e, is_3d=is_3d, source=source)
            except Exception as exc:
                logger.debug("Telemetry logging skipped during st70_to_etrs89: %s", exc)

        # If z value provided, apply 1D grid elevation correction: h = H + N
        if z_seq is not None:
            _degrees = math.degrees
            n_deg = [_degrees(x) for x in r_n]
            e_deg = [_degrees(y) for y in r_e]
            (r_z,) = self._t_gr1d.trans(n_deg, e_deg, z_seq, 1)
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)

    @overload
    def etrs89_to_st70(
        self,
        lat: float,
        lon: float,
        h: None = None,
        source: int = 0,
    ) -> Coord2D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: float,
        lon: float,
        h: float,
        source: int = 0,
    ) -> Coord3D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: Sequence[float],
        lon: Sequence[float],
        h: None = None,
        source: int = 0,
    ) -> Sequence2D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: Sequence[float],
        lon: Sequence[float],
        h: Sequence[float],
        source: int = 0,
    ) -> Sequence3D: ...

    def etrs89_to_st70(
        self,
        lat: ScalarOrSeq,
        lon: ScalarOrSeq,
        h: Optional[ScalarOrSeq] = None,
        source: int = 0,
    ) -> Union[Coord2D, Coord3D, Sequence2D, Sequence3D]:
        """Transforms ETRS89 geographic coordinates to Stereo70 grid coordinates.

        Supports both single coordinates (floats) and coordinate sequences (lists/tuples).
        :math:`(\\text{Lat}, \\text{Long}[, h]) \\to (N, E[, H])`

        :param lat: Latitude in radians (or sequence of latitudes)
        :type lat: float or Sequence[float]
        :param lon: Longitude in radians (or sequence of longitudes)
        :type lon: float or Sequence[float]
        :param h: Ellipsoidal elevation in meters (or sequence of elevations), optional
        :type h: float or Sequence[float], optional
        :param source: Integer source identifier (e.g. 0=unknown, 1=web_map, 2=rest_api)
        :type source: int, optional
        :return: (n, e[, z]) in Stereo70 coordinate reference system (meters)
        :rtype: Tuple[float, float] or Tuple[float, float, float] or Sequence2D or Sequence3D
        :raises OutOfGridErr: If (lat, lon) falls outside the active grid boundaries
        :raises NoDataGridErr: If subgrid interpolation encounters No Data (NaN or 999) values
        """
        is_scalar = not isinstance(lat, (list, tuple))
        lat_seq: Sequence[float]
        lon_seq: Sequence[float]
        h_seq: Optional[Sequence[float]]

        if is_scalar:
            lat_seq, lon_seq = (float(lat),), (float(lon),)  # type: ignore[arg-type]
            h_seq = (float(h),) if h is not None else None  # type: ignore[arg-type]
        else:
            lat_seq, lon_seq = lat, lon  # type: ignore[assignment]
            h_seq = h  # type: ignore[assignment]

        r_n, r_e = self._p_st70.to_grid(lat_seq, lon_seq)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, -1)
        r_n, r_e = self._t_gr2d.trans(r_n, r_e, 1)

        if self._collector is not None:
            try:
                is_3d = 1 if h is not None else 0
                self._collector.record(2, lat_seq, lon_seq, is_3d=is_3d, source=source)
            except Exception as exc:
                logger.debug("Telemetry logging skipped during etrs89_to_st70: %s", exc)

        # If h value provided, apply 1D grid elevation correction: H = h - N
        if h_seq is not None:
            _degrees = math.degrees
            lat_deg = [_degrees(x) for x in lat_seq]
            lon_deg = [_degrees(y) for y in lon_seq]
            (r_z,) = self._t_gr1d.trans(lat_deg, lon_deg, h_seq, -1)
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)
