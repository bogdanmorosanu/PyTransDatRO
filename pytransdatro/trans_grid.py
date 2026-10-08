"""This module stores classes for coordinate transformation using grids.
Grid values and bounds are dynamically loaded from the active .spg binary grid 
via SpgReader.

Usage:
The purpose of trans_grid is to be referenced in the trans_ro module as part 
of the Stereo70 <-> ETRS89 transformation.
Code example:
    t_grid1D = Grid1D()
    z_out, = t_grid1D.trans(lat_in, lon_in, z_in, corr_sgn)

    t_grid2D = Grid2D()
    n_out, e_out = t_grid2D.trans(n_in, e_in, corr_sgn)

Classes:
    - Grid  # abstract class
        - BiInterp
    - Grid1D(Grid)
    - Grid2D(Grid)
"""

import math
import abc
import functools
from pytransdatro import exceptions
from pytransdatro.spg_reader import SpgReader


class Grid(abc.ABC):
    """An abstract class to be used as parent class for Grid1D and Grid2D
    concrete classes.
    """

    def __init__(self, filename=None):
        """Initializes grid data and bounds from SpgReader."""
        self.reader = SpgReader(filename)
        self.file_name = self.reader.grid_filename
        self.source = self.reader.grid_path
        
        self.v_size = self._get_v_size
        self._set_bounds()
        
        self.c_count = round((self.e_max - self.e_min) / self.e_step) + 1
        self.r_count = round((self.n_max - self.n_min) / self.n_step) + 1
        self.nodes_count = self.c_count * self.r_count
        self.sg_size = 4   # interpolation subgrid size (4 columns by 4 rows)
        self.no_data = 999

    @abc.abstractmethod
    def _set_bounds(self):
        """Method to be implemented by subclasses to set bounds from the SPG reader."""
        pass

    def coo_at_idx(self, idx):
        """Returns node's coordinates at specified grid index (0 based index)
                
        :param idx: grid index at which the coordinates will be returned
        :type idx: int

        :return: the coordinates at the specified grid index
        :rtype: tuple of float
        """
        if 0 <= idx < self.nodes_count:
            return (
                self.n_min + idx // self.c_count % self.r_count * self.n_step,
                self.e_min + idx % self.c_count * self.e_step
            )
        else:
            raise exceptions.OutOfRangeIndexGridErr(idx, self)

    def idx_at_coo(self, n, e):
        """Returns the index of the node which is at the specified location.
        The location must coincide with the node location, not be inside a 
        grid cell.

        :param n: coordinate on north direction
        :type n: float

        :param e: coordinate on east direction
        :type e: float 

        :return: grid index at input location
        :rtype: int
        """
        if (n >= self.n_min and n <= self.n_max 
            and e >= self.e_min and e <= self.e_max):
            
            c_idx = (e - self.e_min) / self.e_step
            r_idx = (n - self.n_min) / self.n_step
            
            if c_idx.is_integer() and r_idx.is_integer():
                return int(r_idx) * self.c_count + int(c_idx)
            else:
                raise exceptions.OutOfGridErr(n, e, self)
        else:
            raise exceptions.OutOfGridErr(n, e, self)

    @property
    @abc.abstractmethod
    def _get_v_size(self):
        pass

    def covered_by_grid(self, n, e):
        """Returns true if the interpolation subgrid (4x4 nodes centered in the 
        input location) is covered by the grid, false otherwise.
        """
        n_sg_dist = (self.sg_size / 2 - 1) * self.n_step
        e_sg_dist = (self.sg_size / 2 - 1) * self.e_step

        return (n - n_sg_dist > self.n_min 
                and n + n_sg_dist < self.n_max
                and e - e_sg_dist > self.e_min
                and e + e_sg_dist < self.e_max)    

    def _sgrid_idxs(self, n, e):
        """Returns the indexes of the subgrid required for bicubic interpolation."""
        r_idx = int((n - self.n_min - self.n_step * (self.sg_size / 2 - 1)) / self.n_step)
        c_idx = int((e - self.e_min - self.e_step * (self.sg_size / 2 - 1)) / self.e_step)
        c_count = self.c_count

        return (
            r_idx * c_count + c_idx,
            r_idx * c_count + c_idx + 1,
            r_idx * c_count + c_idx + 2,
            r_idx * c_count + c_idx + 3,
            (r_idx + 1) * c_count + c_idx,
            (r_idx + 1) * c_count + c_idx + 1,
            (r_idx + 1) * c_count + c_idx + 2,
            (r_idx + 1) * c_count + c_idx + 3,
            (r_idx + 2) * c_count + c_idx,
            (r_idx + 2) * c_count + c_idx + 1,
            (r_idx + 2) * c_count + c_idx + 2,
            (r_idx + 2) * c_count + c_idx + 3,
            (r_idx + 3) * c_count + c_idx,
            (r_idx + 3) * c_count + c_idx + 1,
            (r_idx + 3) * c_count + c_idx + 2,
            (r_idx + 3) * c_count + c_idx + 3
        )

    @abc.abstractmethod
    def _grid_vs_at_idxs(self, idxs):
        """Returns grid values from SPG memory for given indices."""
        pass

    def _grid_vs_no_data(self, values):
        """Returns true if a value from values stores a No Data value (NaN or 999)."""
        for v in values:
            if any(math.isnan(x) or x == self.no_data for x in v):
                return True
        return False

    def _reduce_to_unity(self, n, e):
        """Reduce input location to coordinates relative to the unity cell."""
        return (((n - self.n_min) % self.n_step) / self.n_step,
                ((e - self.e_min) % self.e_step) / self.e_step)

    @functools.lru_cache(maxsize=8192) 
    def _init_interp(self, values):
        """Function used to cache BiInterp class instances."""
        return self.BiInterp(values)

    def interp(self, n, e):
        """Returns arrays of interpolated grid values at input locations (sequences)."""
        if self.v_size == 1:
            z_out = []
            for n_val, e_val in zip(n, e):
                if not self.covered_by_grid(n_val, e_val):
                    raise exceptions.OutOfGridErr(n_val, e_val, self)

                sg_idxs = self._sgrid_idxs(n_val, e_val)
                sg_v_dict = self.get_vs_at_idxs_cached(sg_idxs)
            
                if self._grid_vs_no_data(sg_v_dict.values()):
                    raise exceptions.NoDataGridErr(n_val, e_val, self)

                n_unity, e_unity = self._reduce_to_unity(n_val, e_val)

                sg_v_z_list = tuple(v[0] for v in sg_v_dict.values())  
                bi = self._init_interp(sg_v_z_list)
                interp_v_z = bi.interp(n_unity, e_unity)
                z_out.append(interp_v_z)
            return (z_out,)        
        
        elif self.v_size == 2:
            n_out = []
            e_out = []
            for n_val, e_val in zip(n, e):
                if not self.covered_by_grid(n_val, e_val):
                    raise exceptions.OutOfGridErr(n_val, e_val, self)

                sg_idxs = self._sgrid_idxs(n_val, e_val)
                sg_v_dict = self.get_vs_at_idxs_cached(sg_idxs)
            
                if self._grid_vs_no_data(sg_v_dict.values()):
                    raise exceptions.NoDataGridErr(n_val, e_val, self)

                n_unity, e_unity = self._reduce_to_unity(n_val, e_val)

                sg_v_e_list, sg_v_n_list = zip(*sg_v_dict.values())

                bi_n = self._init_interp(sg_v_n_list)
                shift_n = bi_n.interp(n_unity, e_unity)

                bi_e = self._init_interp(sg_v_e_list)
                shift_e = bi_e.interp(n_unity, e_unity)  

                n_out.append(shift_n)
                e_out.append(shift_e)

            return (n_out, e_out)
    class BiInterp():
        """Class which computes the 16 coefficients of the bicubic interpolation
        polynomial and performs interpolation on a 4x4 subgrid.
        """
        def __init__(self, g):
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

            self.__a = (
                df[0],
                df[4],
                -3 * df[0] + 3 * df[1] - 2 * df[4] - df[5],
                2 * df[0] - 2 * df[1] + df[4] + df[5],
                df[8],
                df[12],
                -3 * df[8] + 3 * df[9] - 2 * df[12] - df[13],
                2 * df[8] - 2 * df[9] + df[12] + df[13],
                -3 * df[0] + 3 * df[2] - 2 * df[8] - df[10],
                -3 * df[4] + 3 * df[6] - 2 * df[12] - df[14],
                (9 * df[0] - 9 * df[1] - 9 * df[2] + 9 * df[3] + 6 * df[4] + 3
                 * df[5] - 6 * df[6] - 3 * df[7] + 6 * df[8] - 6 * df[9] + 3
                 * df[10] - 3 * df[11] + 4 * df[12] + 2 * df[13] + 2 * df[14]
                 + df[15]),
                (-6 * df[0] + 6 * df[1] + 6 * df[2] - 6 * df[3] - 3 * df[4] - 3
                 * df[5] + 3 * df[6] + 3 * df[7] - 4 * df[8] + 4 * df[9] - 2
                 * df[10] + 2 * df[11] - 2 * df[12] - 2 * df[13] - df[14]
                 - df[15]),
                2 * df[0] - 2 * df[2] + df[8] + df[10],
                2 * df[4] - 2 * df[6] + df[12] + df[14],
                (-6 * df[0] + 6 * df[1] + 6 * df[2] - 6 * df[3] - 4 * df[4] - 2
                 * df[5] + 4 * df[6] + 2 * df[7] - 3 * df[8] + 3 * df[9] - 3
                 * df[10] + 3 * df[11] - 2 * df[12] - df[13] - 2 * df[14]
                 - df[15]),
                (4 * df[0] - 4 * df[1] - 4 * df[2] + 4 * df[3] + 2 * df[4] + 2
                 * df[5] - 2 * df[6] - 2 * df[7] + 2 * df[8] - 2 * df[9] + 2
                 * df[10] - 2 * df[11] + df[12] + df[13] + df[14] + df[15])
            )

        def interp(self, n_unity, e_unity):
            """Interpolates value at given location (unity coordinates)."""
            r = 0.0
            for i in range(4):
                for j in range(4):
                    r += self.__a[i * 4 + j] * (e_unity ** j) * (n_unity ** i)
            return r



class Grid1D(Grid):       

    def __init__(self, filename=None):
        super().__init__(filename)
        
        # Dynamically bind the transformation method once based on grid metadata
        method = self.reader.interp_vertical
        if method == 0:  # INTERP_COLOCATE (Nearest-Neighbor)
            self.trans = self._trans_colocate
        elif method == 2:  # INTERP_BICUBIC (Bicubic Spline)
            self.trans = self._trans_bicubic
        else:
            raise NotImplementedError(
                f"Unsupported vertical interpolation strategy code: {method}"
            )

    def _set_bounds(self):
        self.e_min, self.e_max, self.n_min, self.n_max, self.e_step, self.n_step = self.reader.height_bounds

    def _grid_vs_at_idxs(self, idxs):
        return {idx: (self.reader.height_flat[idx],) for idx in idxs}

    @property
    def _get_v_size(self):
        return 1

    @functools.lru_cache(maxsize=4096)
    def get_vs_at_idxs_cached(self, idxs):
        return self._grid_vs_at_idxs(idxs)

    def _trans_colocate(self, n, e, z, corr_sgn):
        """Transforms z by adding/subtracting grid height anomaly using
        nearest-neighbor (colocate) node lookup for sequences.
        """        
        c_count = self.c_count
        r_count = self.r_count
        e_min = self.e_min
        n_min = self.n_min
        e_step = self.e_step
        n_step = self.n_step
        height_flat = self.reader.height_flat

        z_out = []
        for n_val, e_val, z_val in zip(n, e, z):
            c_idx = int(round((e_val - e_min) / e_step))
            r_idx = int(round((n_val - n_min) / n_step))
            
            if c_idx < 0: c_idx = 0
            elif c_idx >= c_count: c_idx = c_count - 1
            
            if r_idx < 0: r_idx = 0
            elif r_idx >= r_count: r_idx = r_count - 1
            
            idx = r_idx * c_count + c_idx
            corr = height_flat[idx]
            z_out.append(z_val + corr_sgn * corr)
        return (z_out,)

    def _trans_bicubic(self, n, e, z, corr_sgn):
        """Transforms z by adding/subtracting grid height anomaly using
        bicubic spline polynomial interpolation on a 4x4 subgrid for sequences.
        """
        corrs = self.interp(n, e)
        return ([z_val + corr_sgn * corr for z_val, corr in zip(z, corrs[0])],)


class Grid2D(Grid):

    def __init__(self, filename=None):
        super().__init__(filename)

    def _set_bounds(self):
        self.e_min, self.e_max, self.n_min, self.n_max, self.e_step, self.n_step = self.reader.shift_bounds

    def _grid_vs_at_idxs(self, idxs):
        return {idx: (self.reader.shift_e_flat[idx], self.reader.shift_n_flat[idx]) for idx in idxs}

    @property
    def _get_v_size(self):
        return 2

    @functools.lru_cache(maxsize=4096)    
    def get_vs_at_idxs_cached(self, idxs):
        return self._grid_vs_at_idxs(idxs)        

    def trans(self, n, e, corr_sgn):
        """Transforms n, e by adding/subtracting interpolated values for sequences."""
        corrs = self.interp(n, e)
        return ([n_val + corr_sgn * corr for n_val, corr in zip(n, corrs[0])], 
                [e_val + corr_sgn * corr for e_val, corr in zip(e, corrs[1])])
