"""This module stores classes for the coordinate conversion of the 
Oblique Stereographic projection used by PyTransDatRO. The projection's 
parameters values were extracted from the document 
Help_TransDatRO_code_source_EN.pdf and are assigned when a StereoProj object is
instantiated 

Usage:
The purpose of proj_stereo is to be  referenced in the trans_ro module as part 
of the Stere70 <-> ETRS89 transformation.
Code example:
    st_proj = StereoProj()
    lat, lon = st_proj.to_geo(n, e)   # (Lat, Long) -> (N, E)
    n, e = st_proj.to_grid(lat, lon)  # (N, E) -> (Lat, Long) 

Notes:
    - geographic cooordinates (lat, lon) are expressed in radians
    - grid (projected) coordinates (n, e) are expressed in meters

Classes:
    - StereoProj
        - _Ellipsoid
        - _ConfSph
 
"""
import math
from pytransdatro.utils import sexa_to_rad

class StereoProj():
    """Class which provides functions for Oblique Stereographic projection
    coordinate conversion.
    (Lat, Long) -> (N, E)
    (N, E) -> (Lat, Long)
    """
    # constants which define the WGS84 ellipsoid used by this projection
    WGS84_ELL_SMAJ_AXIS = 6378137   # ellipsoid's semi-major axis 
    WGS84_ELL_INV_FLAT = 298.257223563   # ellipsoid's inverse flattening
    
    # constants which define the projection parameters
    ORIG_LAT = '46 0 0.0'   # latitude of natural origin
    ORIG_LON = '25 0 0.0'   # longitude of natural origin
    FALSE_N = 500000   # false northing
    FALSE_E = 500000   # false easting
    PROJ_SCALE = 0.99975   # scale factor

    def __init__(self):
        """Contructor method which sets the parameters of this projection.
        """
        self._ell = self._Ellipsoid(StereoProj.WGS84_ELL_SMAJ_AXIS,
                                    StereoProj.WGS84_ELL_INV_FLAT)
        self._orig_lat = sexa_to_rad(StereoProj.ORIG_LAT)
        self._orig_lon = sexa_to_rad(StereoProj.ORIG_LON)
        self._f_n = StereoProj.FALSE_N
        self._f_e = StereoProj.FALSE_E
        self._s = StereoProj.PROJ_SCALE
        self.__conf_sphere = self._ConfSph(self)

    def to_grid(self, lat, lon):
        """Returns the grid coordinates (projected) of the input geodetic
        (geographic) coordinates (sequences).
        (Lat, Long) -> (N, E)

        :param lat: sequence of latitudes
        :type lat: sequence of floats

        :param lon: sequence of longitudes
        :type lon: sequence of floats

        :return: (n_out, e_out) lists of northing and easting coordinates
        :rtype: tuple of lists
        """
        lat_c_arr, lon_c_arr = self.__conf_sphere.get_latlon_from_geo(lat, lon)
        lat_c0 = self.__conf_sphere.orig_lat
        lon_c0 = self.__conf_sphere.orig_lon
        sin_lat_c0 = math.sin(lat_c0)
        cos_lat_c0 = math.cos(lat_c0)
        f_n = self._f_n
        f_e = self._f_e
        r2s = 2 * self.__conf_sphere.r * self._s
        _sin = math.sin
        _cos = math.cos

        n_out = []
        e_out = []
        for lat_c, lon_c in zip(lat_c_arr, lon_c_arr):
            sin_lat_c = _sin(lat_c)
            cos_lat_c = _cos(lat_c)
            dlon = lon_c - lon_c0
            b = 1 + sin_lat_c * sin_lat_c0 + cos_lat_c * cos_lat_c0 * _cos(dlon)
            n_out.append(f_n + r2s * (sin_lat_c * cos_lat_c0 - cos_lat_c * sin_lat_c0 * _cos(dlon)) / b)
            e_out.append(f_e + r2s * cos_lat_c * _sin(dlon) / b)
        return n_out, e_out

    def to_geo(self, n, e):
        """Returns the geodetic (geographic) coordinates of the input grid
        (projected) coordinates (sequences).
        (N, E) -> (Lat, Long)

        :param n: sequence of northings
        :type n: sequence of floats

        :param e: sequence of eastings
        :type e: sequence of floats

        :return: (lat_out, lon_out) lists of latitudes and longitudes
        :rtype: tuple of lists
        """
        lat_c_arr, lon_c_arr = self.__conf_sphere.get_latlon_from_grid(n, e)

        e_first = self._ell.first_ecc
        e2 = self._ell.first_ecc2
        c_conf = self.__conf_sphere.c
        n_conf = self.__conf_sphere.n
        orig_lon = self._orig_lon

        _sin = math.sin
        _cos = math.cos
        _log = math.log
        _exp = math.exp
        _atan = math.atan
        _tan = math.tan
        _pow = math.pow
        pi_2 = math.pi / 2
        pi_4 = math.pi / 4
        tol = 0.0000000000484814

        lat_out = []
        lon_out = []
        for lat_c, lon_c in zip(lat_c_arr, lon_c_arr):
            sin_lat_c = _sin(lat_c)
            iso_lat = (0.5 * _log((1 + sin_lat_c) / (c_conf * (1 - sin_lat_c)))) / n_conf
            r_lat = 2 * _atan(_exp(iso_lat)) - pi_2
            diff = tol + 1
            i = 1
            while diff >= tol and i < 50:
                r_lat_before_next_iter = r_lat
                sin_rlat = _sin(r_lat)
                iso_lat_i = _log(_tan(r_lat / 2 + pi_4) * _pow((1 - e_first * sin_rlat) / (1 + e_first * sin_rlat), e_first / 2))
                r_lat = r_lat - ((iso_lat_i - iso_lat) * _cos(r_lat) * (1 - e2 * sin_rlat * sin_rlat) / (1 - e2))
                diff = abs(r_lat_before_next_iter - r_lat)
                i += 1
            lat_out.append(r_lat)
            lon_out.append(orig_lon + (lon_c - orig_lon) / n_conf)
        return lat_out, lon_out

    class _Ellipsoid:
        """Class which provides functionality for computation of basic 
        ellipsoid's elements.
        """
        def __init__(self, smaj_axis, inv_flat):
            """
            :param smaj_axis: semi-major axis of the ellipsoid
            :type smaj_axis: float

            :param inv_flat: inverse flattening of the ellipsoid (1/flattening)
            :type inv_flat: float

            :ivar flat: flattening
            :ivar smin_axis: flattening
            :ivar first_ecc2: # first eccentricity squared
            :ivar first_ecc: first eccentricity
            :ivar sec_ecc2: second eccentricity squared
            :ivar sec_ecc: second eccentricity
            """
            self.smaj_axis = smaj_axis
            self.inv_flat = inv_flat
            self.flat = 1 / inv_flat
            self.smin_axis = smaj_axis * (1 - 1 / inv_flat)
            self.first_ecc2 = 2 * self.flat - self.flat * self.flat
            self.first_ecc = math.sqrt(self.first_ecc2)
            self.sec_ecc2 = self.first_ecc2 / (1 - self.first_ecc2)
            self.sec_ecc = math.sqrt(self.sec_ecc2)
        
        def get_rad_M_N(self, lat):
            """Calculates the radius of curvature in meridian - north-south 
            direction (M) and the radius of curvature in prime vertical - 
            east-west direction (N)

            :param lat: the latitude for which M and N values are computed
            :type lat: float

            :return: The values of M and N - (M,N)
            :rtype: tuple of floats
            """
            # temporary value used for computation of both M and N
            # sqrt(1 - e^2 * sin^2(lat)), where e = first_ecc
            tmp_v = math.sqrt((1 - self.first_ecc2 * math.pow(math.sin(lat), 2)))
            m = self.smaj_axis * (1 - self.first_ecc2) / math.pow(tmp_v, 3)
            n = self.smaj_axis / tmp_v
            return m, n

    class _ConfSph():
        """Class defining a conformal sphere.
        Provides methods to compute the equivalent conformal latitude and 
        longitude of a given point
        """
        def __init__(self, proj):
            """
            :param proj: the projection associated with the conformal sphere
            :type proj: StereoOblProj 
            """

            # calculate parameters defining the conformal sphere (tagged CSP)
            ## some vars used in calculus (precompute values to optimize
            ## computation and improve readability)
            r_m, r_n = proj._ell.get_rad_M_N(proj._orig_lat)
            e2 = proj._ell.first_ecc2
            e = proj._ell.first_ecc
            sin_orig_lat = math.sin(proj._orig_lat)
            s1 = (1 + sin_orig_lat) / (1 - sin_orig_lat)
            s2 = (1 - e * sin_orig_lat) / (1 + e * sin_orig_lat)
            r = math.sqrt(r_m * r_n)   # CSP
            n = math.sqrt(1 + (e2 * math.pow(math.cos(proj._orig_lat), 4) 
                               / (1 - e2)))   # CSP
            w1 = math.pow(s1, n) * math.pow(s2, e * n)
            sin_chi0 = (w1 - 1) / (w1 + 1)
            c = (n + sin_orig_lat) * (1 - sin_chi0) / ((n - sin_orig_lat) 
                                                       * (1 + sin_chi0))   # CSP
            w2 = c * w1
            
            self.__proj = proj
            self.r = r
            self.n = n
            self.c = c
            self.orig_lat = math.asin((w2 - 1) / (w2 + 1))
            self.orig_lon = proj._orig_lon

        def get_latlon_from_grid(self, n, e):
            """Computes the equivalent conformal latitude and longitude for a 
            given point with stereographic grid coordinates (N,E)

            :param n: northing
            :type n: float

            :param e: easting
            :type e: float 

            :return: the values of latitude and longitude (lat,lon)
            :rtype: tuple of floats               
            """

            # some vars used in calculus (precompute values to optimize 
            # computation and improve readability) 
            sf = self.__proj._s   # projection scale factor
            r = self.r
            g = 2 * r * sf * math.tan((math.pi / 4) - (self.orig_lat / 2))
            h = 4 * r * sf * math.tan(self.orig_lat) + g
            f_n = self.__proj._f_n
            f_e = self.__proj._f_e
            orig_lat = self.orig_lat
            orig_lon = self.orig_lon
            r2s = 2 * r * sf
            _atan = math.atan
            _tan = math.tan

            lat_out = []
            lon_out = []
            for n_val, e_val in zip(n, e):
                dn = n_val - f_n
                de = e_val - f_e
                i = _atan(de / (h + dn))
                j = _atan(de / (g - dn)) - i
                lat_out.append(orig_lat + 2 * _atan((dn - de * _tan(j / 2)) / r2s))
                lon_out.append(j + 2 * i + orig_lon)
            return lat_out, lon_out

        def get_latlon_from_geo(self, lat, lon):
            """Computes the equivalent conformal latitude and longitude for 
            sequences of geodetic geographic coordinates (latitudes, longitudes).
                
            :param lat: sequence of latitudes in radians
            :type lat: sequence of floats

            :param lon: sequence of longitudes in radians
            :type lon: sequence of floats

            :return: (lat_out, lon_out) lists of conformal latitudes and longitudes
            :rtype: tuple of lists
            """
            e = self.__proj._ell.first_ecc
            n = self.n
            orig_lon = self.orig_lon
            c = self.c
            _sin = math.sin
            _asin = math.asin
            _pow = math.pow

            lat_out = []
            lon_out = []
            for lat_val, lon_val in zip(lat, lon):
                sin_lat = _sin(lat_val)
                sa = (1 + sin_lat) / (1 - sin_lat)
                sb = (1 - e * sin_lat) / (1 + e * sin_lat)
                w = c * _pow((sa * _pow(sb, e)), n)
                lat_out.append(_asin((w - 1) / (w + 1)))
                lon_out.append(n * (lon_val - orig_lon) + orig_lon)
            return lat_out, lon_out




