"""Module which stores functions to transform coordinates between Stereo70 and
ETRS89 coordinate reference systems (CRS).
Stere70 is the official cartographic projection used in Romania

Functions list:
    - st70_to_etrs89: (N, E) -> (Lat, Long)
    - etrs89_to_st70: (Lat, Long) -> (N, E)
"""

import trans_grid as t_gr
import trans_helmert2d as t_h2d
import proj_stereo as p_st

import utils

class TransRO():
    """Class which stores functions to transform coordinates between 
    Stereo70 and ETRS89 coordinate reference systems
    """
    def __init__(self):
        self._t_gr = t_gr.Grid2D()  # 2D grid transformation ETRS89_KRASOVSCHI42_2DJ
        
        # 2D Helmert transformation Stereo70 to StereoGRS80
        self._t_h2d = t_h2d.Helmert2D(False)  # TODO: init somehow helmert to work both ways

        self._p_st70 = p_st.StereoProj()



    def st70_to_etrs89(self, n, e):
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
        # apply 2D grid transformation
        r_n, r_e = self._t_gr.trans(n, e, -1)

        # apply 2D Helmert transformation
        r_n, r_e = self._t_h2d.trans(r_n, r_e)

        # apply stereografic projection
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)

        return (r_n, r_e)


t = TransRO()
n = 693771.731
e = 310723.518
lat, lon = t.st70_to_etrs89(n, e)
print(utils.rad_to_sexa(lat), utils.rad_to_sexa(lon))