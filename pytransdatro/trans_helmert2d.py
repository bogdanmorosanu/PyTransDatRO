"""This module stores a class for the Helmert 2D (4 parameters) coordinate
transformation. The transformation parameter values are dynamically read from 
the SPG grid metadata via SpgReader when a Helmert2D object is instantiated.

Usage:
The purpose of Helmert2D is to be part of the Stereo70 <-> ETRS89 transformation
defined in the module trans_ro.
Code example:
    h2D = Helmert2D()
    n, e = h2D.trans(n, e, 1) 

Notes:
    - in and out coordinates (n, e) are expressed in meters

Classes:
    - Helmert2D
"""
import math
from pytransdatro.spg_reader import SpgReader

class Helmert2D():
    """Class which calculates the 2D Helmert transformation between Stereo70 and 
    Stereo on GRS80
    """
    
    def __init__(self, tn=None, te=None, ppm=None, rot_rad=None):
        """Constructor sets parameters for Stereo70 to StereoGRS80 transformation.
        If parameters are omitted, they are dynamically loaded from the active SPG grid.
        """
        reader = SpgReader()
        self.__tn = reader.helmert_tn if tn is None else float(tn)
        self.__te = reader.helmert_te if te is None else float(te)
        self.__ppm = reader.helmert_ppm if ppm is None else float(ppm)
        self.__r = reader.helmert_rot_rad if rot_rad is None else float(rot_rad)
        
        # Precompute constants to optimize speed for batch calculations
        self.__sinr = math.sin(self.__r)
        self.__cosr = math.cos(self.__r)
        self.__scale_factor = self.__ppm * 1E-6

    @property
    def tn(self):
        return self.__tn

    @property
    def te(self):
        return self.__te

    @property
    def ppm(self):
        return self.__ppm

    @property
    def r(self):
        return self.__r

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
        scale = 1.0 + sign * self.__scale_factor
        return (
            (n * scale * self.__cosr 
            + e * scale * sign * self.__sinr 
            + sign * self.__tn),
            (e * scale * self.__cosr 
            - n * scale * sign * self.__sinr
            + sign * self.__te)
        )
