"""Module which stores functions to transform coordinates between Stereo70 and
ETRS89 coordinate reference systems.
Stereo70 is the official cartographic projection used in Romania.

Functions list:
    - st70_to_etrs89: (N, E[, H]) -> (Lat, Long[, h])
    - etrs89_to_st70: (Lat, Long[, h]) -> (N, E[, H])
"""
import math
import pytransdatro.trans_grid 
import pytransdatro.trans_helmert2d 
import pytransdatro.proj_stereo 
import pytransdatro.exceptions
import pytransdatro.utils

class TransRO():
    """Class which stores functions to transform coordinates between 
    Stereo70 and ETRS89 coordinate reference systems.
    """
    def __init__(self, grid_filename=None, telemetry_logger=None):
        """Initialize the trans_ro coordinate transformation.
        
        :param grid_filename: Optional filename or path of the .spg grid.
            If None, the grid is auto-discovered from pytransdatro/grids/.
        :ivar _t_gr2d: 2D grid transformation
        :ivar _t_h2d: 2D Helmert transformation Stereo70 to StereoGRS80
        :ivar _p_st70: stereographic oblique projection on WGS84 ellipsoid
        :ivar _t_gr1d: 1D grid transformation
        """
        self._t_gr2d = pytransdatro.trans_grid.Grid2D(grid_filename)
        self._t_h2d = pytransdatro.trans_helmert2d.Helmert2D() 
        self._p_st70 = pytransdatro.proj_stereo.StereoProj()
        self._t_gr1d = pytransdatro.trans_grid.Grid1D(grid_filename)
        self._collector = telemetry_logger
        
    def st70_to_etrs89(self, n, e, z=None, source=0):
        """Transforms Stereo70 grid coordinates to ETRS89 geographic coordinates.
        Supports both single coordinates (floats) and sequences (lists/tuples).
        (N, E[, H]) -> (Lat, Long[, h])

        :param n: northing (meters) or sequence of northings
        :type n: float or sequence of floats

        :param e: easting (meters) or sequence of eastings
        :type e: float or sequence of floats

        :param z: normal elevation (meters) or sequence of elevations, optional
        :type z: float or sequence of floats, optional

        :param source: integer source ID (e.g. 0=unknown, 1=web_map, 2=rest_api), optional
        :type source: int, optional

        :return: (lat, lon[, h]) in ETRS89 coordinate reference system
        :rtype: tuple of floats or tuple of lists

        :raises OutOfGridErr: if (n,e) is out of the grid / not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s).
        """
        is_scalar = not isinstance(n, (list, tuple))
        if is_scalar:
            n, e = (n,), (e,)
            if z is not None:
                z = (z,)

        r_n, r_e = self._t_gr2d.trans(n, e, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, 1)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)

        if self._collector is not None:
            try:
                is_3d = 1 if z is not None else 0
                self._collector.record(1, r_n, r_e, is_3d=is_3d, source=source)
            except Exception:
                pass

        # If z value provided, apply 1D grid elevation correction: h = H + N
        if z is not None:
            _degrees = math.degrees
            n_deg = [_degrees(x) for x in r_n]
            e_deg = [_degrees(y) for y in r_e]
            r_z, = self._t_gr1d.trans(n_deg, e_deg, z, 1)
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)
    
    def etrs89_to_st70(self, lat, lon, h=None, source=0):
        """Transforms ETRS89 geographic coordinates to Stereo70 grid coordinates.
        Supports both single coordinates (floats) and sequences (lists/tuples).
        (Lat, Long[, h]) -> (N, E[, H])

        :param lat: latitude (radians) or sequence of latitudes
        :type lat: float or sequence of floats

        :param lon: longitude (radians) or sequence of longitudes
        :type lon: float or sequence of floats

        :param h: ellipsoidal elevation (meters) or sequence of elevations, optional
        :type h: float or sequence of floats, optional

        :param source: integer source ID (e.g. 0=unknown, 1=web_map, 2=rest_api), optional
        :type source: int, optional

        :return: (n, e[, z]) in Stereo70 coordinate reference system
        :rtype: tuple of floats or tuple of lists

        :raises OutOfGridErr: if (n,e) is out of the grid / not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s).
        """        
        is_scalar = not isinstance(lat, (list, tuple))
        if is_scalar:
            lat, lon = (lat,), (lon,)
            if h is not None:
                h = (h,)

        if self._collector is not None:
            try:
                is_3d = 1 if h is not None else 0
                self._collector.record(0, lat, lon, is_3d=is_3d, source=source)
            except Exception:
                pass

        r_n, r_e = self._p_st70.to_grid(lat, lon)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, -1)
        r_n, r_e = self._t_gr2d.trans(r_n, r_e, 1)

        # If h value provided, apply 1D grid elevation correction: H = h - N
        if h is not None:
            _degrees = math.degrees
            lat_deg = [_degrees(x) for x in lat]
            lon_deg = [_degrees(y) for y in lon]
            r_z, = self._t_gr1d.trans(lat_deg, lon_deg, h, -1)
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)
