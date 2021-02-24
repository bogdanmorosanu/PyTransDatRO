import math

class Ellipsoid:
    """
    Class which provides basic ellipsoid element computation.

    """
    def __init__(self, smaj_axis, inv_flat):
        """
        :param smaj_axis: semi-major axis of the ellipsoid
        :type smaj_axis: float

        :param inv_flat: inverse flattening of the ellipsoid (1 / flattening)
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
        """
        Calculates the radius of curvature in meridian - north-south direction (M)
        and the radius of curvature in prime vertical - east-west direction (N)

        :param lat: the latitude for which M and N values are computed
        :type lat: float

        :return: The values of M and N
        :rtype: tuple of floats
        """
        # temporary value used for computation of both M and N
        # sqrt(1 - e^2 * sin^2(lat)), where e - first_ecc2
        tmp_v = math.sqrt((1 - self.first_ecc2 * math.pow(math.sin(lat), 2)))
        m = self.smaj_axis * (1 - self.first_ecc2) / math.pow(tmp_v, 3)
        n = self.smaj_axis / tmp_v
        return m, n

class StereoOblProj():
    """
    Class which provides functions for Oblique Stereographic projection
    coordinate conversion.
    (Lat, Long) -> (N, E)
    (N, E) -> (Lat, Long)
    """
    def __init__(self, ell, orig_lat, orig_lon, f_n, f_e, s):
        """
        :param ell: projection's ellipsoid
        :type ell: Ellipsoid

        :param orig_lat: projection's latitude of natural origin
        :type orig_lat: float

        :param orig_lon: projection's longitude of natural origin
        :type orig_lon: float

        :param f_n: projection's false northing
        :type f_n: float

        :param f_e: projection's false easting
        :type f_e: float

        :param s: projection's scale factor
        :type s: float    
        """

        self.__ell = ell
        self.__orig_lat = orig_lat
        self.__orig_lon = orig_lon
        self.__f_n = f_n
        self.__f_e = f_e
        self.__s = s


class ConfSph():
    """
    Class defining a conformal sphere
    """
    def __init__(self, ell, proj_orig_lat, proj_orig_lon,
                 proj_f_n, proj_f_e, proj_s):
        """
        :param ell: projection's ellipsoid
        :type ell: Ellipsoid

        :param proj_orig_lat: projection's geodetic origin latitude
        :type proj_orig_lat: float

        :param proj_orig_lon: projection's geodetic origin longitude
        :type proj_orig_lon: float 

        :param f_n: projection's false northing
        :type f_n: float

        :param f_e: projection's false easting
        :type f_e: float

        :param s: projection's scale factor
        :type s: float   
        """

        # calculate parameters defining the conformal sphere (tagged CSP)
        ## some vars used in calculus (precompute values to optimeze and 
        ## improve readability)
        r_m, r_n = ell.get_rad_M_N(proj_orig_lat)
        e2 = ell.first_ecc2
        e = ell.first_ecc
        sin_orig_lat = math.sin(proj_orig_lat)
        s1 = (1 + sin_orig_lat) / (1 - sin_orig_lat)
        s2 = (1 - e * sin_orig_lat) / (1 + e * sin_orig_lat)
        r = math.sqrt(r_m * r_n)   # CSP
        n = math.sqrt(1 + (e2 * math.pow(math.cos(proj_orig_lat), 4) / (1 - e2)))   # CSP
        w1 = math.pow(s1, n) * math.pow(s2, e * n)
        sin_chi0 = (w1 - 1) / (w1 + 1)
        c = (n + sin_orig_lat) * (1 - sin_chi0) / ((n - sin_orig_lat) * (1 + sin_chi0))   # CSP
        w2 = c * w1
        
        self.__f_n = proj_f_n
        self.__f_e = proj_f_e
        self.__s = proj_s
        self.__ell_first_ecc = ell.first_ecc
        self.__r = r
        self.__n = n
        self.__c = c
        self.__orig_lat = math.asin((w2 - 1) / (w2 + 1))
        self.__orig_lon = proj_orig_lon

    def get_latlon_from_grid(self, n, e):
        """
        Computes the equivalent conformal latitude and longitude for a given 
        point with Stereographic grid coordinates (N,E)

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :return: The values of Latitude and Longitude
        :rtype: tuple of floats               
        """

        # some vars used in calculus (precompute values to optimeze and 
        # improve readability)
        g = 2 * self.__r * self.__s * math.tan((math.pi / 4) - (self.__orig_lat / 2))
        h = 4 * self.__r * self.__s * math.tan(self.__orig_lat) + g
        dn = n - self.__f_n
        de = e - self.__f_e
        i = math.atan(de / (h + dn))
        j = (math.atan(de / (g - dn))) - i  

        return (self.__orig_lat + 2 * math.atan((dn - de * math.tan(j / 2))
                                                / (2 * self.__r * self.__s)), 
                j + 2 * i + self.__orig_lon)

    def get_latlon_from_geo(self, lat, lon):
        """
        Computes the equivalent conformal latitude and longitude for a given
        point with geodetic geographic coordinates (latitude, longitude)
              
        :param lat: latitude
        :type lat: float

        :param lon: longitude
        :type lon: float 

        :return: The values of Latitude and Longitude
        :rtype: tuple of floats         
        """

        # some vars used in calculus (precompute values to optimeze and 
        # improve readability)
        sin_lat = math.sin(lat)
        sa = (1 + sin_lat) / (1 - sin_lat)
        sb = (1 - self.__ell_first_ecc * sin_lat) / (1 + self.__ell_first_ecc * sin_lat)
        w = self.__c * math.pow((sa * math.pow(sb, self.__ell_first_ecc)), self.__n)

        return (self.__n * (lat - self.__orig_lon) + self.__orig_lon,
                math.asin((w - 1) / (w + 1)))       


class Helmert2D():
    """
    Class which does exposes a 2D Helmert transformation
    """

    def __init__(self, tx, ty, s, r):
        """
        Define the parameters of the transformation

        :param tx: translation X
        :type tx: float

        :param ty: translation Y
        :type ty: float   

        :param s: unity scale, not ppm
        :type s: float

        :param r: rotation value in radians
        :type r: float                             
        """

        self.__tx = tx
        self.__ty = ty
        self.__sxsinr = math.sin(r)   # s * sin(r)
        self.__sxcosr = math.cos(r)   # s * cos(r)
    
    def trans(self, x, y):
        """
        Calculates the new x', y' values of x and y input 
        by using the transfomration's parameters

        :param x: x
        :type x: float

        :param y: Y
        :type y: float 

        :return: transformed x' and y' values
        :rtype: tuple of floats                 
        """
        return (x * self.__sxcosr + y * self.__sxsinr + self.__ty,
                y * self.__sxcosr - x * self.__sxsinr + self.__tx)

