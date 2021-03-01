import math
from pytransdatro.utils import ppm_to_unity
from pytransdatro.utils import sexa_to_rad

class Helmert2D():
    """
    Class which defines a 2D Helmert transformation
    """

    # constants to define the Stereo70 to Stereo on GRS80 Helmert transformation
    T_N = -31.8051   # translation on North
    T_E = -119.7358   # translation on East
    PPM_SCALE = -0.11559991   # ppm scale factor
    ROTATION = '0 0 0.22739706'   # latitude of natural origin
    

    def __init__(self, reversed):
        """
        :param reversed: initialize the reversed transformation
        :type reversed: bool 
        """
        if reversed:
            self.__t_e = -Helmert2D.T_E
            self.__t_n = -Helmert2D.T_N
            tmp_s = 2 - ppm_to_unity(Helmert2D.PPM_SCALE)
            tmp_r = - sexa_to_rad(Helmert2D.ROTATION)
        else:
            self.__t_e = Helmert2D.T_E 
            self.__t_n = Helmert2D.T_N 
            tmp_s = ppm_to_unity(Helmert2D.PPM_SCALE)
            tmp_r = sexa_to_rad(Helmert2D.ROTATION)

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
