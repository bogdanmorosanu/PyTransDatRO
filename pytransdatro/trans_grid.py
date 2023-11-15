"""This module stores classes for coordinate transformation using grids.
The grid files (files which store the grid values) were copied from the
official distribtion of the TransDatRO application (download link available
here: https://rompos.ro/index.php/download/category/2-software).
The grid files are stored in the /grid folder and referenced when Grid1D/Grid2D
objects are instantiated.

Usage:
The purpose of trans_grid is to be  referenced in the trans_ro module as part 
of the Stere70 <-> ETRS89 transformation.
Code example:
    t_grid1D = Grid1D()
    z_out = t_grid1D.trans(lat_in, lon_in, z_in, corr_sgn)

    t_grid2D = Grid2D()
    # z value not included
    n_out, e_out = t_grid1D.trans(n_in, e_in, corr_sgn)
    # z value included
    n_out, e_out, z_out = t_grid1D.trans(n_in, e_in, z_in, corr_sgn)

Notes:
    - coordinates used as input and output by functions of this module are named
    n and e even though they don't necessarily describe easting and northing 
    (for example, they could be latitudes and longitudes). The only meaning 
    which can be associated with them is the axis direction:
        - n: S-N/vertical direction
        - e: W-E/horizontal direction

Classes:
    - Grid  # abstract class
        - BiInterp
    - Grid1D(Grid)
    - Grid2D(Grid)
"""

import struct
import math
import abc
import pkg_resources
import functools
from pytransdatro import exceptions

_GRID_DIR = 'grids'   # grid files directory name

class Grid(abc.ABC):
    """An abstract class to be used as parent class for Grid1D and Grid2D
    concrete classes.
    This class defines a grid as a collection of nodes (locations) which store
    values (grid values). It can store a single value (Grid1D concrete class) or two values
    (concrete class Grid2D). Each node is described by position and location.
    Its position is described by an index as shown below:

    (n-1)m   (n-1)m+1  (n-1)m+2   ...   nm-1
       ...       ...       ...    ...    ...
        2m       2m+1      2m+2   ...   3m-1
         m       m+1        m+2   ...   2m-1
         0       1          2     ...    m-1
    
    The above grid has m x n nodes, the first being at index 0 and the last one
    at index n x m - 1.
    Its location is computed based on the location of the node at index 0. This
    can be done as distances between nodes in a row/column are equal. 
    The distance between them is called Step E, for rows or W to E direction, 
    and Step N, for columns or S to N direction.
    So, the nodes of a grid cover an area which is defined by the most 
    southwestern (SW) node and the most northeastern (NE) node.
    The grid transformation consists of values being added/substracted to the
    input coordinates. The values applied are computed based on a bicubic 
    interpolation which uses a subgrid of 4 by 4 nodes centered in the input
    location.
    A grid cell is defined by 4 nodes which are next to each other
    """

    def __init__(self, file_name):
        """Contructor method which sets the grid file reference and reads
        grid definition data from its header.

        :ivar source: the grid's source file. This is the binary file which is 
        distributed with every official realease of TransDatRO application 
        (https://rompos.ro/index.php/download/category/2-software) 
        :ivar v_size: the grid's values count (1 for Grid1D and 2 for Grid2D)
        :ivar n_min, e_min, n_max, e_max: min and max ccordinates for grid area
        :ivar n_step, e_step: Step N and Step E values
        :ivar c_count: grid's number of columns (m-1)
        :ivar r_count: grid's number of rows (n)
        :ivar sg_size: interpolation subgrid size (sg_size x sg_size nodes)
        :ivar no_data: grid value which indicates that there is no data
            available for the corresponding node
        
        :raises FileNotFoundError: if grid file not found
        :raises IOError: if it fails to read the header content of the grid file  
        """
        self.file_name = file_name
        self.source = pkg_resources.resource_filename(__name__, 
                                               f'{_GRID_DIR}/{file_name}')
        
        self.v_size = self._get_v_size
        try:
            # read grid header
            with open(self.source, 'rb') as f:
                self.e_min, self.e_max = struct.unpack('<dd', f.read(16))
                self.n_min, self.n_max = struct.unpack('<dd', f.read(16))
                self.e_step, self.n_step = struct.unpack('<dd', f.read(16))
        
        except FileNotFoundError:
            raise FileNotFoundError(f'Failed to open grid file: {self.source}')
        except (OSError, IOError):
            raise IOError(f'Failed to read the grid file content: {self.source}')
        self.c_count = round((self.e_max - self.e_min) 
                            / self.e_step) + 1   # grid columns count
        self.r_count = round((self.n_max - self.n_min) 
                            / self.n_step) + 1   # grid rows count
        self.sg_size = 4   # interpolation subgrid size (4 columns by 4 rows)

    @abc.abstractproperty
    def _get_v_size(self):
        """Abstract property to get the grid's number of values for each node
        """
        pass

    def covered_by_grid(self, n, e):
        """
        Returns true if the interpolation subgrid (4x4 nodes centered in the 
        input location) is covered by the grid, false otherwise
        ne ± d must be inside the grid mins and maxs values to return True
        ---------
        |   d   |
        |   |   |
        |d-n,e-d|
        |   |   |
        |   d   |
        ---------

        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float 

        :return: interpolation subgrid is covered by the grid area
        :rtype: bool
        """
        # distances (d) from n, e input at which a location should still be
        #  covered by the grid such that interpolated values can be computed
        n_sg_dist = (self.sg_size / 2 - 1) * self.n_step
        e_sg_dist = (self.sg_size / 2 - 1) * self.e_step

        return (n - n_sg_dist > self.n_min 
                and n + n_sg_dist < self.n_max
                and e - e_sg_dist > self.e_min
                and e + e_sg_dist < self.e_max)    

    def _sgrid_idxs(self, n, e):
        """Returns the indexes of the subgrid required for bicupic interpolation
        The returned subgrid size is 4x4, centered on the input n and e

        13 14 15 16
        9  10 11 12
            n,e - input location 
        4  5  7  8
        1  2  3  4

        .. warning:: the function doesn't check if the input location is
        covered by the grid. Wrong result returned if used in this scenario
        
        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float 

        :return: the indexes of the subgrid
        :rtype: tuple of int
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

    def _grid_vs_at_idxs(self, idxs):
        """Returns grid values at specified indexes.

        :param idxs: indexes from which grid values will be returned.
        :type idxs: tuple of int

        :return: grid values at input indexes
        :rtype: dict (key = position/index, value = tuple of float(s))         
        """
        r = {}   # result dict
        tmp_format = f'<{"d" * self.v_size}'
        tmp_buffer = self.v_size * 8
        with open(self.source, 'rb') as f:
            # 0, 1, 2 indexes are used by the grid header
            # 48 = 8(bits) * 2(values) * 3(sets of values)
            start_bit = 48   
            for idx in idxs:
                f.seek(start_bit + idx * tmp_buffer)
                r[idx] = struct.unpack(tmp_format, f.read(tmp_buffer))
        return r

    def _grid_vs_no_data(self, values):
        """Returns true if a value from values stores a No Data value,
        false otherwise

        :param values: values which will be checked for No Data value
        :type values: list of tuple(float[, float])

        :return: True if input contains No Data value, false otherwise
        :rtype: bool
        """
        for v in values:
            if self.no_data in v:
                return True
        return False

    def _reduce_to_unity(self, n, e):
        """Reduce input location to coordinates relative to the unity cell.
        The origin of the returned coordinates is defined by the SV node of the
        cell in which the input location resides. The unit of measure is 
        relative to the Step N and Step E values so that the result coordinates 
        are within the interval [0, 1]

        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float

        :return: The unity coordinates of input location
        :rtype: tuple of floats
        """
        return (((n - self.n_min) % self.n_step) / self.n_step,
                ((e - self.e_min) % self.e_step) / self.e_step)

    @functools.lru_cache(maxsize = 3072) 
    def _init_interp(self, values):
        """Function used to cache the BiInterp class instances. This is usefull
        for locations which are within the same grid cell. In these cases the
        bicubic interpolation coefficients are the same

        :param values: 4x4 grid values used for interpolation
        :type values: tuple

        :return: an instance of BiInterp initialized with coefficients 
            derived from the input values
        :rtype: BiInterp
        """
        bi = self.BiInterp(values)
        return bi

    def interp (self, n, e):
        """Returns interpolated grid values at input location

        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float 

        :return: The interpolated value(s) (one or two if Grid1D or Grid2D 
            respectively)
        :rtype: tuple of floats

        :raises OutOfGridErr: if interpolation subgrid not covered by the grid
        :raises NoDataGridErr: if interpolation subgrid has No Data values
        """
        if not self.covered_by_grid(n, e):
            raise  exceptions.OutOfGridErr(n, e, self)

        sg_idxs = self._sgrid_idxs(n, e)
        sg_v_dict = self.get_vs_at_idxs_cached(sg_idxs)
       
        if self._grid_vs_no_data(sg_v_dict.values()):
            raise exceptions.NoDataGridErr(n, e, self)

        n_unity, e_unity = self._reduce_to_unity(n, e)

        # if instance of Grid1D
        if self.v_size == 1:
            sg_v_z_list = tuple(v[0] for v in sg_v_dict.values())  
            bi = self._init_interp(sg_v_z_list)
            interp_v_z = bi.interp(n_unity, e_unity)
            return (interp_v_z, )        

        # if instance of Grid2D
        if self.v_size == 2:
            sg_v_e_list, sg_v_n_list = zip(*sg_v_dict.values())

            bi = self._init_interp(sg_v_n_list)
            shift_n = bi.interp(n_unity, e_unity)

            bi = self._init_interp(sg_v_e_list)
            shift_e = bi.interp(n_unity, e_unity)  
            return (shift_n, shift_e)             
    
    @abc.abstractmethod
    def trans(self):
        """Abstrasct function to define the transformation method of the grid
        """
        pass

    @abc.abstractmethod
    def get_vs_at_idxs_cached(self):
        """Function used to cache values read from the grid.
        """
        pass

    class BiInterp():
        """Class which has a function to compute a bicubic spline interpolation.
        """
        def __init__(self, sg_values):
            """Constructor method which initializes the coefficients used in
            the interpolation calculus

            :param sg_values: interpolation subgrid values
            :type sg_values: list of floats  
            """ 
            self.__a = {}      
            self.__init_coefs(sg_values)

        def __init_coefs(self, g):
            """Initialize the coeficients of the bicubic interpolation.

            :param g: interpolation subgrid values
            :type g: list of floats          
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
        
        def interp(self, n_unity, e_unity):
            """ Interpolates value at given location (unity coordinates)
            
            :param n_unity: unity coordinate on north direction
            :type n_unity: float 

            :param e_unity: unity coordinate on east direction
            :type e_unity: float              
            """         
            p={}
            for i in range(4):
                for j in range(4):
                    p[i * 4 + j] = math.pow(e_unity, j) * math.pow(n_unity, i)
            
            r = 0.0   # result value
            for k in range(16):
                r += self.__a[k] * p[k]
            return r  

    def __str__(self):
        """Returns a string representation of a Grid instance
        """
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

class Grid1D(Grid):       

    @property
    def _get_v_size(self):
        """Returns the number of grid values stored by each node 
        """        
        return 1

    @functools.lru_cache(maxsize = 2048)
    def get_vs_at_idxs_cached(self, idxs):
        """Cached method used to return grid nodes values. Returned grid values
        are the same for coordinates within the same grid cell.

        :param idxs: indexes from which grid values will be returned.
        :type idxs: tuple of int

        :return: grid values at input indexes
        :rtype: dict (key = position/index, value = tuple of float) 
        """              
        return super()._grid_vs_at_idxs(idxs)

    def trans(self, n, e, z, corr_sgn):
        """Transforms z by adding/substracting interpolated value 
        (correction) obtained from the grid at the input location
        
        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float

        :param corr_sgn: 1 to add the grid correction, -1 to substract it  
        :type corr_sgn: int

        :return: the input z coordinate with the grid correction applied
        :rtype: tuple of float
        """        
        corrs = self.interp(n, e)
        return (z + corr_sgn * corrs[0],)          
         
class Grid2D(Grid):

    @property
    def _get_v_size(self):
        """Returns the number of grid values stored by each node 
        """
        return 2

    @functools.lru_cache(maxsize = 1024)    
    def get_vs_at_idxs_cached(self, idxs):
        """Cached method used to return grid nodes values. Returned grid values
        are the same for coordinates within the same grid cell.

        :param idxs: indexes from which grid values will be returned.
        :type idxs: tuple of int

        :return: grid values at input indexes
        :rtype: dict (key = position/index, value = tuple of floats)   
        """
        return super()._grid_vs_at_idxs(idxs)        

    def trans(self, n, e, corr_sgn):
        """Transforms n, e by adding/substracting interpolated values 
        (corrections) obtained from the grid at the input location
        
        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float

        :param corr_sgn: 1 to add the grid corrections, -1 to substract them  
        :type corr_sgn: int

        :return: the input coordinates with the grid corrections applied
        :rtype: tuple of floats
        """
        corrs = self.interp(n, e)
        return (n + corr_sgn * corrs[0], e + corr_sgn * corrs[1])         

