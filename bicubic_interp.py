import math

class BiInterp():
    """
    Class which does a bicubic spline interpolation.
    Must be first initialized with grid values (values used for interpolation)
  
    """
    def __init__(self, g):
        """
        :param g: grid values
        :type g: Ellipsoid  
        """ 
        self.__a = {}      
        self.init_coefs(g)

    def init_coefs(self, g):

        df = {}
        df[0] = g[5]
        df[1] = g[6]
        df[2] = g[9]
        df[3] = g[10]

        # Derivatives in the East direction and the North direction
        df[4] = (-g[7] + 4 * g[6] - 3 * g[5]) / 2
        df[5] = (3 * g[6] - 4 * g[5] + g[4]) / 2
        df[6] = (-g[11] + 4 * g[10] - 3 * g[9]) / 2
        df[7] = (3 * g[10] - 4 * g[9] + g[8]) / 2
        df[8] = (-g[13] + 4 * g[9] - 3 * g[5]) / 2
        df[9] = (-g[14] + 4 * g[10] - 3 * g[6]) / 2
        df[10] = (3 * g[9] - 4 * g[5] + g[1]) / 2
        df[11] = (3 * g[10] - 4 * g[6] + g[2]) / 2

        # Equations for the cross derivative
        df[12] = ((g[0] + g[10]) - (g[2] + g[8])) / 4
        df[13] = ((g[1] + g[11]) - (g[3] + g[9])) / 4
        df[14] = ((g[4] + g[14]) - (g[6] + g[12])) / 4
        df[15] = ((g[5] + g[15]) - (g[7] + g[13])) / 4

        self.__a[0] = df[0]
        self.__a[1] = df[4]
        self.__a[2] = -3 * df[0] + 3 * df[1] - 2 * df[4] - df[5]
        self.__a[3] = 2 * df[0] - 2 * df[1] + df[4] + df[5]
        self.__a[4] = df[8]
        self.__a[5] = df[12]
        self.__a[6] = -3 * df[8] + 3 * df[9] - 2 * df[12] - df[13]
        self.__a[7] = 2 * df[8] - 2 * df[9] + df[12] + df[13]
        self.__a[8] = -3 * df[0] + 3 * df[2] - 2 * df[8] - df[10]
        self.__a[9] = -3 * df[4] + 3 * df[6] - 2 * df[12] - df[14]
        self.__a[10] = (9 * df[0] - 9 * df[1] - 9 * df[2] + 9 * df[3] + 6 * df[4] + 3
                        * df[5] - 6 * df[6] - 3 * df[7] + 6 * df[8] - 6 * df[9] + 3
                        * df[10] - 3 * df[11] + 4 * df[12] + 2 * df[13] + 2 * df[14]
                        + df[15])
        self.__a[11] = (-6 * df[0] + 6 * df[1] + 6 * df[2] - 6 * df[3] - 3 * df[4] - 3
                        * df[5] + 3 * df[6] + 3 * df[7] - 4 * df[8] + 4 * df[9] - 2
                        * df[10] + 2 * df[11] - 2 * df[12] - 2 * df[13] - df[14]
                        - df[15])
        self.__a[12] = 2 * df[0] - 2 * df[2] + df[8] + df[10]
        self.__a[13] = 2 * df[4] - 2 * df[6] + df[12] + df[14]
        self.__a[14] = (-6 * df[0] + 6 * df[1] + 6 * df[2] - 6 * df[3] - 4 * df[4] - 2
                        * df[5] + 4 * df[6] + 2 * df[7] - 3 * df[8] + 3 * df[9] - 3
                        * df[10] + 3 * df[11] - 2 * df[12] - df[13] - 2 * df[14]
                        - df[15])
        self.__a[15] = (4 * df[0] - 4 * df[1] - 4 * df[2] + 4 * df[3] + 2 * df[4] + 2
                        * df[5] - 2 * df[6] - 2 * df[7] + 2 * df[8] - 2 * df[9] + 2
                        * df[10] - 2 * df[11] + df[12] + df[13] + df[14] + df[15])
    
    def interp(self, x, y):
        """
        Interpolates value at given x and y
        
        :param g: grid values
        :type g: Ellipsoid  
        """         
        p={}
        for i in range(4):
            for j in range(4):
                p[i * 4 + j] = math.pow(x, j) * math.pow(y, i)
        
        r = 0.0   # result value
        for k in range(16):
            r += self.__a[k] * p[k]
        return r

