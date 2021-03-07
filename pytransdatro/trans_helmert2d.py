"""This module stores a class for the Helmert 2D (4 parameters) coordinate
tranformation. The transformation's parameters values were extracted from the
document Help_TransDatRO_code_source_EN.pdf and are assigned when a Helmert2D 
object is instantiated 

Usage:
The purpose of Helmert2D is to be part of the Stere70 <-> ETRS89 transformation
defined in the module trans_ro. There should be no purpose of using it in
another context.
Code example:
    h2D = Helmert2D()
    n, e = h2D.trans(n, e) 

Notes:
    - in and out coordinates (n, e) are expressed in meters

Classes:
    - Helmert2D
"""
import math
from pytransdatro import utils

class Helmert2D():
    """Class which calculates the 2D Helmert transformation between Stereo70 and 
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
        self.__tn = Helmert2D.T_N   # translation on north
        self.__te = Helmert2D.T_E   # translation on east
        self.__ppm = Helmert2D.PPM  # scale (ppm)
        self.__r = utils.sexa_to_rad(Helmert2D.ROTATION)   # rotation (radians)
        
        # some precomputed values to optimize speed
        self.__sinr = math.sin(self.__r)   # sin(r)
        self.__cosr = math.cos(self.__r)   # cos(r)

    def trans(self, n, e, sign):
        """Transforms the values of n and e
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

