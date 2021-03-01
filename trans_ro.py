"""Module which stores functions to transform coordinates between Stereo70 and
ETRS89 coordinate reference systems (CRS).
Stere70 is the official cartographic projection used in Romania

Functions list:
    - st70_to_etrs89: (N, E) -> (Lat, Long)
    - etrs89_to_st70: (Lat, Long) -> (N, E)
"""
import math
import trans_grid as t_gr
import trans_helmert2d as t_h2d
import proj_stereo as p_st

import utils

class TransRO():
    """Class which stores functions to transform coordinates between 
    Stereo70 and ETRS89 coordinate reference systems
    """
    def __init__(self):
        self._t_gr2d = t_gr.Grid2D()  # 2D grid transformation ETRS89_KRASOVSCHI42_2DJ
        # 2D Helmert transformation Stereo70 to StereoGRS80
        self._t_h2d = t_h2d.Helmert2D(False)  # TODO: init somehow helmert to work both ways
        self._p_st70 = p_st.StereoProj()
        self._t_gr1d = t_gr.Grid1D()  # 1D grid transformation EGG97_QGRJ

    def st70_to_etrs89(self, n, e, z = None):
        """Transforms Stereo70 grid coordinates to ETRS89 geopgraphic coordinates
        (N, E) -> (Lat, Long)

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :return: (Latitude,Longitude) in ETRS89 coordinate reference
            system
        :rtype: tuple of float

        :raises OutOfGridErr: if (N, E) is out of the grid/not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (N, E) has
            No Data value(s). In general, this would apply to locations which are
            outside of Romania's border        
        """
        r_n, r_e = self._t_gr2d.transform(n, e, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)
        # IF z value provided, apply 1D grid transformation (height correction)
        if z is not None:
            r_z, = self._t_gr1d.transform(math.degrees(r_n), 
                                          math.degrees(r_e),
                                          z, 1)
            return (r_n, r_e, r_z)
        return (r_n, r_e)
    
    def etrs89_to_st70(self, lat, lon, h = None):
        r_n, r_e = self._p_st70.to_grid(lat, lon)
    
        h2d = t_h2d.Helmert2D(True)   # TODO: asta trebe stearsa trebe rezolva helmert reveresed
        r_n, r_e = h2d.trans(r_n, r_e)
        r_n, r_e = self._t_gr2d.transform(r_n, r_e, 1)

        # IF h value provided, apply 1D grid transformation (height correction)
        if h is not None:
            r_z, = self._t_gr1d.transform(math.degrees(lat), 
                                          math.degrees(lon),
                                          h, -1)
            return (r_n, r_e, r_z)
        return (r_n, r_e)
        
t = TransRO()
n = 693771.731
e = 310723.518
z = 122.714
lat, lon, h = t.st70_to_etrs89(n, e, z)
print(utils.rad_to_sexa(lat), utils.rad_to_sexa(lon), h)
print(t.etrs89_to_st70(lat, lon, h))
