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
    def __init__(self, grid_filename=None):
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
        
    def st70_to_etrs89(self, n, e, z=None):
        """Transforms Stereo70 grid coordinates to ETRS89 geographic coordinates:
        (N, E[, H]) -> (Lat, Long[, h])

        :param n: northing (meters)
        :type n: float

        :param e: easting (meters)
        :type e: float 

        :param z: normal elevation (meters), optional
        :type z: float, optional

        :return: (lat, lon[, h]) in ETRS89 coordinate reference system
        :rtype: tuple of floats

        :raises OutOfGridErr: if (n,e) is out of the grid / not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s).
        """
        r_n, r_e = self._t_gr2d.trans(n, e, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, 1)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)

        # If z value provided, apply 1D grid elevation correction: h = H + N
        if z is not None:
            n_deg = math.degrees(r_n)
            e_deg = math.degrees(r_e)
            r_z, = self._t_gr1d.trans(n_deg, e_deg, z, 1)
            return (r_n, r_e, r_z)
        else:
            return (r_n, r_e)
    
    def etrs89_to_st70(self, lat, lon, h=None):
        """Transforms ETRS89 geographic coordinates to Stereo70 grid coordinates:
        (Lat, Long[, h]) -> (N, E[, H])

        :param lat: latitude (radians)
        :type lat: float

        :param lon: longitude (radians)
        :type lon: float 

        :param h: ellipsoidal elevation (meters), optional
        :type h: float, optional

        :return: (n, e[, z]) in Stereo70 coordinate reference system
        :rtype: tuple of floats

        :raises OutOfGridErr: if (n,e) is out of the grid / not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s).
        """        
        r_n, r_e = self._p_st70.to_grid(lat, lon)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, -1)
        r_n, r_e = self._t_gr2d.trans(r_n, r_e, 1)

        # If h value provided, apply 1D grid elevation correction: H = h - N
        if h is not None:
            lat_deg = math.degrees(lat)
            lon_deg = math.degrees(lon)
            r_z, = self._t_gr1d.trans(lat_deg, lon_deg, h, -1)
            return (r_n, r_e, r_z)
        else:
            return (r_n, r_e)
