import struct
import time
import math
import exceptions as ex_td  # Transdat exceptions

class Grid2D():
    """
    TODO: describe the grid, how is indexed (positions), step, min, max
     - mention the neighbours used for interpolation - used by point is covered by grid
    :ivar source: the grid source file which stores the shift values     
    """
    GRID2D_SRC = r'ETRS89_KRASOVSCHI42_2DJ.GRD'
    def __init__(self):
        self.source = Grid2D.GRID2D_SRC
        try:
            # read grid header
            with open(self.source, 'rb') as f:
                self.e_min, self.e_max = struct.unpack('<dd', f.read(16))
                self.n_min, self.n_max = struct.unpack('<dd', f.read(16))
                self.e_step, self.n_step = struct.unpack('<dd', f.read(16))
        except IOError:
            raise IOError(f'Failed to open file: {self.source}')
        except Exception:
            raise Exception(f'Failed to read the content of binary grid file: {self.source}')
        self.c_count = round((self.e_max - self.e_min) 
                            / self.e_step) + 1   # grid columns count
        self.r_count = round((self.n_max - self.n_min) 
                            / self.n_step) + 1   # grid rows count
        self.sg_size = 4   # interpolation subgrid size (4 columns by 4 rows)

        self.no_data = 999.0

    def is_inside_grid(self, n, e):
        """
        Returns True if the input location (N,E) is inside (covered) by the grid,
        False otherwise. It also checks the neighbouring grid nodes 
        (subgrid 4x4 centered in the input N,E) which is required for 
        computing the interpolated grid shift values.
        
        NE ± d must be inside the grid mins and maxs values to return True
        ---------
        |   d   |
        |   |   |
        |d-N,E-d|
        |   |   |
        |   d   |
        ---------

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :return: (N,E) location and required subgrid is inside grid
        :rtype: bool

        """
        # distances (d) from input N, E at defining a location which should
        # still be covered by the grid so interpolated values can be computed
        n_sg_dist = (self.sg_size / 2 - 1) * self.n_step
        e_sg_dist = (self.sg_size / 2 - 1) * self.e_step

        return (n - n_sg_dist >= self.n_min 
                and n + n_sg_dist <= self.n_max
                and e - e_sg_dist >= self.e_min
                and e + e_sg_dist <= self.e_max)

    def get_inter_sgrid_idxs(self, n, e):
        """
        Returns the indexes of a subgrid required for bicupic interpolation
        The returned subgrid size is 4x4, centered on the input N and E

        13 14 15 16
        9  10 11 12
            N,E - input location 
        4  5  7  8
        1  2  3  4

        Note:
        The function doesn't check if the (N,E) is covered by the grid. Wrong
        result is returned if used in this case
        
        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :return: The indexes of the subgrid
        :rtype: tuple of indexes values
        """
        # calculate indexes of the bottom left corner (position 1)
        r_idx = int((n - self.n_min - self.n_step * (self.sg_size / 2 - 1)) 
                    / self.n_step)
        c_idx = int((e - self.e_min - self.e_step * (self.sg_size / 2 - 1)) 
                    / self.e_step)

        # result ( len = subgrid nodes count)
        r = [None] * (self.sg_size * self.sg_size)   
        for i in range(4):
            for j in range(4):
                r[i * self.sg_size + j] = (r_idx + i) * self.c_count + c_idx + j
        return tuple(r)

    def get_shifts_at_idxs(self, idxs):
        """
        Returns values from grid at specified indexes. Each index returns two
        values, i.e. values are stored in pairs of 2 (8 bits for each value),
        N shift and E shift

        :param idxs: indexes from which grid values will be returned.
            Indexes are 0 index based and start with the first set 
            (header ignored) of shift values (shift E and shif N)
        :type idxs: set

        :return: Grid (shift) values for the input indexes
        :rtype: dict (key = pos index, value = tuple of shift values for E and N)         
        """
        r = {}   # result dict
        with open(self.source, 'rb') as f:
            # 0, 1, 2 indexes are used by the grid header
            # 48 = 8(bits) * 2(values) * 3(sets of values)
            start_bit = 48   
            for idx in (idxs):
                f.seek(start_bit + idx * 16)
                r[idx] = struct.unpack('<dd', f.read(16))
        return r

    def values_have_no_data(self, values):
        """
        Returns True if a value from the grid data stores a No Data value,
        False otherwise

        :param values: values which will be checked for No Data value
        :type values: list of tuple(float, float)

        :return: True if a No Data value is found
        :rtype: boolean
        """
        for v in values:
            if self.no_data in v:
                return True
        return False

    def reduce_to_unity(self, n, e):
        """Reduce input (N,E) to a location relative to the unity cell of 
        the grid.
        Result coordinates are within the interval [0, 1]
        TODO: add preciese description of what this does

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float

        :return: The unity coordinates of input (N,E)
        :rtype: tuple of floats
        """
        return (((n - self.n_min) % self.n_step) / self.n_step,
                ((e - self.e_min) % self.e_step) / self.e_step)

    # TODO: cache: cred ca aici e (o sa mai vad)
    def init_interp(self, values):
        bi = self.BiInterp(values)
        return bi

    def trans (self, n, e, corr_sgn):
        """
        Transforms input (N,E) based on grid shift values 
        and bicubic interpolation

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param corr_sgn: Determines how the shift grid corrections are applied.
            Can have 2 values only. 1, corrections will be added and -1,
            corrections will be substracted
        :type corr_sgn: int 

        :return: The transformed values of northing and easting (N, E)
        :rtype: tuple of floats

        :raises OutOfGridErr: if (N, E) is out of the grid/not covered by grid
        :raises NoDataGridErr: if subgrid used for interpolation at (N, E) has
            No Data value(s)
        """
        if not self.is_inside_grid(n, e):
            raise  ex_td.OutOfGridErr(n, e, self)

        sg_idxs = self.get_inter_sgrid_idxs(n, e)
        sg_v_dict = self.get_shifts_at_idxs(sg_idxs)    # todo: cache

        if self.values_have_no_data(sg_v_dict.values()):
            raise ex_td.NoDataGridErr(n, e, self)

        n_unity, e_unity = self.reduce_to_unity(n, e)
        sg_v_n_list = []
        sg_v_e_list = []
        for v in sg_v_dict.values():
            sg_v_n_list.append(v[1])    
            sg_v_e_list.append(v[0])
        
        bi = self.init_interp(sg_v_n_list)
        shift_n = bi.interp(n_unity, e_unity)

        bi = self.init_interp(sg_v_e_list)
        shift_e = bi.interp(n_unity, e_unity)  

        return (n + corr_sgn * shift_n,
                e + corr_sgn * shift_e)     

        



        

        
            


    def __str__(self):
        r = (
            f"e_min={self.e_min}\n"
            f"e_max={self.e_max}\n"
            f"n_min={self.n_min}\n"
            f"n_max={self.n_max}\n"
            f"e_step={self.e_step}\n"
            f"n_step={self.n_step}\n"
            f"r_count={self.r_count}\n"
            f"c_count={self.c_count}\n"
            )
        return r

    class BiInterp():
        """
        Class which does a bicubic spline interpolation.
        Must be first initialized with grid values (values used for interpolation)
    
        """
        def __init__(self, sg_values):
            """
            :param sg_values: interpolation subgrid values
            :type sg_values: list of floats  
            """ 
            self.__a = {}      
            self.__init_coefs(sg_values)

        def __init_coefs(self, g):
            """
            Initialize the coeficients of the bicubic interpolation.

            :param sg_values: interpolation subgrid values
            :type sg_values: list of floats          
            """
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
                    p[i * 4 + j] = math.pow(y, j) * math.pow(x, i)
            
            r = 0.0   # result value
            for k in range(16):
                r += self.__a[k] * p[k]
            return r        

# file_content = np.fromfile('ETRS89_KRASOVSCHI42_2DJ.GRD', np.float64)
# print(file_content)
my_grid = Grid2D()
#  print(my_grid)

# test point in grid
# n = 219135
# e = 115173
# print(my_grid.is_ne_covered_by_grid(n, e))



# test point in grid
# n = 239135
# e = 128173
# print(my_grid.get_inter_sgrid_idxs(n, e))


# test get shifts
# n = 757856
# e = 556349
# idxs = my_grid.get_inter_sgrid_idxs(n, e)
# print(idxs)
# print(my_grid.get_shifts_at_idxs(idxs))


# test no data
# n = 768856
# e = 632349
# idxs = my_grid.get_inter_sgrid_idxs(n, e)
# sg_values = my_grid.get_shifts_at_idxs(idxs)
# print(my_grid.values_have_no_data(sg_values.values()))

# test transformation
# n_in = 500000
# e_in = 500000

# n_out, e_out = my_grid.trans(n_in, e_in, 1)
# print(n_out, e_out)






# idxs = set(range(3800))
# # print(idxs)
# start_time = time.time()
# values = my_grid.get_shifts_at_idxs(idxs)
# print(len(values))
# print("Process finished --- %s seconds ---" % (time.time() - start_time))



