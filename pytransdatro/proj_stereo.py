"""This module stores classes for the coordinate conversion of the 
Oblique Stereographic projection used by PyTransDatRO. The projection's 
parameters values were extracted from the document 
Help_TransDatRO_code_source_EN.pdf and are assigned when a StereoProj object is
instantiated 

Usage:
The purpose of StereoProj is to be part of the Stere70 <-> ETRS89 transformation
defined in the module trans_ro. There should be no purpose of using it in
another context.
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
        (geographic) coordinates 
        (Lat, Long) -> (N, E)

        :param lat: latitude
        :type lat: float

        :param lon: longitude
        :type lon: float 

        :return: The values of northing and easting (n,e)
        :rtype: tuple of floats           
        """
        lat_c, lon_c = self.__conf_sphere.get_latlon_from_geo(lat, lon)
        lat_c0 = self.__conf_sphere.orig_lat
        lon_c0 = self.__conf_sphere.orig_lon

        # some vars used in calculus (precompute values to optimize computation
        # and improve readability)        
        sin_lat_c = math.sin(lat_c)
        cos_lat_c = math.cos(lat_c)
        sin_lat_c0 = math.sin(lat_c0)
        cos_lat_c0 = math.cos(lat_c0)
        dlon = lon_c - lon_c0
        b = (1 + sin_lat_c * sin_lat_c0 
             + cos_lat_c * cos_lat_c0 * math.cos(dlon))
        
        return ((self._f_n + 2 * self.__conf_sphere.r * self._s 
                               * (sin_lat_c * cos_lat_c0 
                                  - cos_lat_c * sin_lat_c0 * math.cos(dlon)) 
                                / b),
                (self._f_e + 2 * self.__conf_sphere.r 
                               * self._s * cos_lat_c * math.sin(dlon) / b))

    def to_geo(self, n, e):
        """Returns the geodetic (geographic) coordinates of the input grid
        (projected) coordinates 
        (N, E) -> (Lat, Long)

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :return: The values of latitude and longitude  (lat,lon)
        :rtype: tuple of floats           
        """

        lat_c, lon_c = self.__conf_sphere.get_latlon_from_grid(n, e)

        iso_lat = ((0.5 * math.log((1 + math.sin(lat_c))
                   / (self.__conf_sphere.c * (1 - math.sin(lat_c))))) 
                   / self.__conf_sphere.n)
        e = self._ell.first_ecc
        e2 = self._ell.first_ecc2

        # next computation involves an iterative process until 
        # precision is achieved (tolerance and counter used to control the loop)
        tol = 0.0000000000484814   # radians of DMS value: 0 0 0.00001
        i = 0   # counter to avoid infinite loop
        diff = tol + 1
        r_lat = 0.0   # assign a default value for result latitude
        while diff >= tol and i < 50:
            if i == 0:   # if first iteration compute first approximation
                r_lat = 2 * math.atan(math.exp(iso_lat)) - math.pi / 2
            else:
                r_lat_before_next_iter = r_lat
                iso_lat_i = math.log(math.tan(r_lat / 2 + math.pi / 4)
                                     * math.pow((1 - e * math.sin(r_lat))
                                     / (1 + e * math.sin(r_lat)), e / 2))

                r_lat = r_lat - ((iso_lat_i - iso_lat) * math.cos(r_lat)
                                 * (1 - e2 * math.pow(math.sin(r_lat), 2)) 
                                 / (1 - e2))
                diff = abs(r_lat_before_next_iter - r_lat)
            i+=1
        
        return (r_lat,
                self._orig_lon + (lon_c - self._orig_lon) / self.__conf_sphere.n)

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
            g = 2 * self.r * sf * math.tan((math.pi / 4) - (self.orig_lat / 2))
            h = 4 * self.r * sf * math.tan(self.orig_lat) + g
            dn = n - self.__proj._f_n
            de = e - self.__proj._f_e
            i = math.atan(de / (h + dn))
            j = (math.atan(de / (g - dn))) - i  

            return (self.orig_lat + 2 * math.atan((dn - de * math.tan(j / 2))
                                                    / (2 * self.r * sf)), 
                    j + 2 * i + self.orig_lon)

        def get_latlon_from_geo(self, lat, lon):
            """Computes the equivalent conformal latitude and longitude for a 
            given point with geodetic geographic coordinates (latitude, 
            longitude)
                
            :param lat: latitude in radians
            :type lat: float

            :param lon: longitude in radians
            :type lon: float 

            :return: the values of latitude and longitude (lat,lon)
            :rtype: tuple of floats         
            """

            # some vars used in calculus (precompute values to optimize 
            # computation and improve readability)
            e = self.__proj._ell.first_ecc
            sin_lat = math.sin(lat)
            sa = (1 + sin_lat) / (1 - sin_lat)
            sb = (1 - e * sin_lat) / (1 + e * sin_lat)
            w = self.c * math.pow((sa * math.pow(sb, e)), self.n)

            return (math.asin((w - 1) / (w + 1)),
                    self.n * (lon - self.orig_lon) + self.orig_lon)       




