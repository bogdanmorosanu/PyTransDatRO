import math
from utils import ppm_to_unity
from utils import sexa_to_rad

class Helmert2D():
    """
    Class which defines a 2D Helmert transformation
    """

    # constants to define the Stereo70 to Stereo on GRS80 Helmert transformation
    _T_N = -31.8051   # translation on North
    _T_E = -119.7358   # translation on East
    _SCALE = -0.11559991   # ppm scale factor
    _ROTATION = '0 0 0.22739706'   # latitude of natural origin
    

    def __init__(self, reversed):
        """
        :param reversed: initialize the reversed transformation
        :type reversed: bool 
        """
        self.__t_e = Helmert2D._T_E if not reversed else -Helmert2D._T_E
        self.__t_n = Helmert2D._T_N if not reversed else -Helmert2D._T_N
        scale_unity = ppm_to_unity(Helmert2D._SCALE)
        tmp_s = scale_unity if not reversed else 2 - scale_unity
        rotation_rads = sexa_to_rad(Helmert2D._ROTATION)
        tmp_r = rotation_rads if not reversed else - rotation_rads
        self.__sxsinr = tmp_s * math.sin(tmp_r)   # s * sin(r)
        self.__sxcosr = tmp_s * math.cos(tmp_r)   # s * cos(r)
    
    def trans(self, n, e):
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
        return (n * self.__sxcosr + e * self.__sxsinr + self.__t_n,
                e * self.__sxcosr - n * self.__sxsinr + self.__t_e)
