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
    
    def rad_M_N(self, lat):
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

