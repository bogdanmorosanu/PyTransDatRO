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
import pytransdatro.exceptions


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
        :ivar _t_gr1d_buc: 1D grid transformation using grid zitaBucx 
        """
        self._t_gr2d = pytransdatro.trans_grid.Grid2D('ETRS89_KRASOVSCHI42_2DJ.GRD')
        self._t_h2d = pytransdatro.trans_helmert2d.Helmert2D() 
        self._p_st70 = pytransdatro.proj_stereo.StereoProj()
        self._t_gr1d = pytransdatro.trans_grid.Grid1D('EGG97_QGRJ.GRD')
        
        # grid added to match the transformation done by TransdatRO. Normaly,
        # there should be only one 1D grid. But TransdatRO, intorduced another
        # 1D grid for Bucharest area. This change in the computation algorithm 
        # is not documented anywhere, so it was introduced in the py_transdat
        # implementation in the same way so future grid updates can be done
        # just by using copy and paste (another approach would have been to 
        # merge the two grids into one so the algorithm stays consistent)
        # This grid is used in an bounding box defined by min/max lat/lon
        # extracted by checking where the results of using the EGG97_QGRJ.GRD
        # started to be different from the one obtained by using TransDatRO. 
        self._t_gr1d_buc = pytransdatro.trans_grid.Grid1D('zitaBucx.grd')
    
    def _grid1d_sel(self, n_deg, e_deg):
        """Returns the 1D Grid which will be used for the h correction.
        See comments for class member self._t_gr1d_buc for more details.

        :param n_deg: lat (degrees)
        :type n: float

        :param e_deg: lon (degrees)
        :type e: float         

        :return: a reference to the 1D grid relevant for given location
        :rtype: Grid1D       
        """

        if (e_deg >= 25.981719000000002 and e_deg <= 26.223385860000004
            and n_deg >= 44.35033746000002 and n_deg <= 44.517004260000164):
            return self._t_gr1d_buc
        else:
            return self._t_gr1d
       
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
            n_deg = math.degrees(r_n)
            e_deg = math.degrees(r_e)
            grid_1d = self._grid1d_sel(n_deg, e_deg)
            r_z, = grid_1d.trans(n_deg, e_deg, z, 1)
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
            lat_deg = math.degrees(lat)
            lon_deg = math.degrees(lon)
            grid_1d = self._grid1d_sel(lat_deg, lon_deg)
            r_z, = grid_1d.trans(lat_deg, lon_deg, h, -1)
            return (r_n, r_e, r_z)
        else:
            return (r_n, r_e)

