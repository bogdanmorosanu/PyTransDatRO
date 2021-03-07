"""Module which stores functions to transform coordinates between Stereo70 and
ETRS89 coordinate reference systems.
Stere70 is the official cartographic projection used in Romania

Functions list:
    - st70_to_etrs89: (N, E[, H]) -> (Lat, Long[, h])
    - etrs89_to_st70: (Lat, Long[, h]) -> (N, E[, H])
"""
import math
import pytransdatro.trans_grid 
import pytransdatro.trans_helmert2d 
import pytransdatro.proj_stereo 


import pytransdatro.utils

class TransRO():
    """Class which stores functions to transform coordinates between 
    Stereo70 and ETRS89 coordinate reference systems
    """
    def __init__(self):
        """Initialize the trans_ro coordinate transformation
        
        :ivar _t_gr2d: 2D grid transformation using grid ETRS89_KRASOVSCHI42_2DJ
        :ivar _t_h2d: 2D Helmert transformation Stereo70 to StereoGRS80
        :ivar _p_st70: stereographic oblique projection on WGS84 ellipsoid
        :ivar _t_gr1d: 1D grid transformation using grid EGG97_QGRJ      
        """
        self._t_gr2d = pytransdatro.trans_grid.Grid2D()
        self._t_h2d = pytransdatro.trans_helmert2d.Helmert2D() 
        self._p_st70 = pytransdatro.proj_stereo.StereoProj()
        self._t_gr1d = pytransdatro.trans_grid .Grid1D()

    def st70_to_etrs89(self, n, e, z = None):
        """Transforms Stereo70 grid coordinates to ETRS89 geopgraphic coordinates
        (N, E[, H]) -> (Lat, Long[, h])

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param z: elevation
        :type z: float, optional

        :return: (lat, lon[, h]) in ETRS89 coordinate reference system
        :rtype: tuple of floats

        :raises OutOfGridErr: if (n,e) is out of the grid/not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s). In general, this would apply to locations which 
            are outside of Romania's border but still covered by the grid      
        """
        r_n, r_e = self._t_gr2d.trans(n, e, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, 1)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)
        # IF z value provided, apply 1D grid transformation (height correction)
        if z is not None:
            r_z, = self._t_gr1d.trans(math.degrees(r_n), 
                                          math.degrees(r_e),
                                          z, 1)
            return (r_n, r_e, r_z)
        else:
            return (r_n, r_e)
    
    def etrs89_to_st70(self, lat, lon, h = None):
        """Transforms ETRS89 geopgraphic coordinates to Stereo70 grid coordinates
        (Lat, Long[, h]) -> (N, E[, H])

        :param lat: latitude
        :type lat: float

        :param lon: longitude
        :type lon: float 

        :param h: elevation
        :type h: float, optional

        :return: (n, e[, z]) in Stereo70 coordinate reference system
        :rtype: tuple of floats

        :raises OutOfGridErr: if (n,e) is out of the grid/not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (n,e) has
            No Data value(s). In general, this would apply to locations which
            are outside of Romania's border but still covered by the grid        
        """        
        r_n, r_e = self._p_st70.to_grid(lat, lon)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, -1)
        r_n, r_e = self._t_gr2d.trans(r_n, r_e, 1)

        # IF h value provided, apply 1D grid transformation (height correction)
        if h is not None:
            r_z, = self._t_gr1d.trans(math.degrees(lat), 
                                          math.degrees(lon),
                                          h, -1)
            return (r_n, r_e, r_z)
        else:
            return (r_n, r_e)

