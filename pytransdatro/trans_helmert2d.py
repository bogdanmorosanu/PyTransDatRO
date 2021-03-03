import math
from pytransdatro import utils

class Helmert2D():
    """
    Class which calculates the 2D Helmert transformation between Stereo70 and 
    Stereo on GRS80
    """

    # constants to define the Stereo70 to Stereo on GRS80 Helmert transformation
    T_N = -31.8051
    T_E = -119.7358
    PPM = -0.11559991   # ppm value
    ROTATION = '0 0 0.22739706'   # DMS value
    
    def __init__(self):
        """Constructor sets parameters for Stereo70 to StereoGRS80 transformation        
        """
        self.__tn = Helmert2D.T_N   # translation on North
        self.__te = Helmert2D.T_E   # translation on East
        self.__ppm = Helmert2D.PPM  # scale (ppm)
        self.__r = utils.sexa_to_rad(Helmert2D.ROTATION)   # rotation (radians)
        
        # some precomputed values to optimize speed
        self.__sinr = math.sin(self.__r)   # sin(r)
        self.__cosr = math.cos(self.__r)   # cos(r)

    def trans(self, n, e, sign):
        """
        Transforms the values of n and e
        (N,E) -> (N',E')

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param sign: Value of 1 to transform using defined parameters or -1 for 
            the inverse transformation (params x -1)
        :type sign: int     

        :return: the new values of n and e (N',E')
        :rtype: tuple of floats                 
        """       
        return (
            (n * (1 + sign * self.__ppm * 1E-6) * self.__cosr 
            + e * (1 + sign *  self.__ppm * 1E-6) * sign * self.__sinr 
            + sign * self.__tn),
            (e * (1 + sign * self.__ppm * 1E-6) * self.__cosr 
            - n * (1 + sign * self.__ppm * 1E-6) * sign * self.__sinr
            + sign * self.__te)
        )

