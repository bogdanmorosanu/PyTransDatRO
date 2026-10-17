"""Grid-based Coordinate Transformations and Bicubic Interpolation.

This module provides classes for coordinate transformations using binary SPG grids:
    - :class:`Grid2D`: 2D planar distortion grid interpolating horizontal shifts
      :math:`(\\delta N, \\delta E)` between Romanian national networks.
    - :class:`Grid1D`: 1D vertical grid applying quasi-geoid undulation / height anomaly
      corrections :math:`N` (:math:`h = H + N` and :math:`H = h - N`).
    - :class:`Grid.BiInterp`: 16-coefficient bicubic polynomial interpolation engine
      evaluating continuous surfaces across 4x4 node subgrids.

Grid values and dimensions are loaded dynamically from the active `.spg` binary file
via :class:`~pytransdatro.spg_reader.SpgReader`.
"""

import abc
import functools
import math
from typing import Callable, Iterable, Optional, Sequence, Tuple

from pytransdatro import exceptions
from pytransdatro.spg_reader import SpgReader


class Grid(abc.ABC):
    """Abstract base class for 1D and 2D geodetic grid transformation engines."""

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initializes grid boundaries, node counts, and data accessors from SpgReader.

        :param filename: Optional grid filename or path. If None, auto-discovered.
        :type filename: Optional[str]
        """
        self.reader: SpgReader = SpgReader(filename)
        self.file_name: str = self.reader.grid_filename
        self.source: str = self.reader.grid_path

        self.e_min: float = 0.0
        self.e_max: float = 0.0
        self.n_min: float = 0.0
        self.n_max: float = 0.0
        self.e_step: float = 0.0
        self.n_step: float = 0.0

        self.v_size: int = self._get_v_size
        self._set_bounds()

        self.c_count: int = round((self.e_max - self.e_min) / self.e_step) + 1
        self.r_count: int = round((self.n_max - self.n_min) / self.n_step) + 1
        self.nodes_count: int = self.c_count * self.r_count
        self.sg_size: int = 4  # Interpolation subgrid size (4 columns by 4 rows)
        self.no_data: float = 999.0

    @abc.abstractmethod
    def _set_bounds(self) -> None:
        """Sets grid coordinate boundaries and steps from SPG reader metadata."""
        pass

    def coo_at_idx(self, idx: int) -> Tuple[float, float]:
        """Returns node coordinates at specified grid linear index (0-based index).

        :param idx: Linear node index
        :type idx: int
        :return: Tuple of ``(northing, easting)`` node coordinates
        :rtype: Tuple[float, float]
        :raises OutOfRangeIndexGridErr: If idx is negative or exceeds available nodes
        """
        if 0 <= idx < self.nodes_count:
            return (
                self.n_min + idx // self.c_count % self.r_count * self.n_step,
                self.e_min + idx % self.c_count * self.e_step,
            )
        else:
            raise exceptions.OutOfRangeIndexGridErr(idx, self)

    _out_of_grid_exc: type[exceptions.OutOfGridErr] = exceptions.OutOfGridErr

    def idx_at_coo(self, n: float, e: float) -> int:
        """Returns the linear index of the node that coincides with specified coordinates.

        The coordinates must fall precisely on a grid node, not inside a grid cell.

        :param n: Northing coordinate (or latitude in degrees)
        :type n: float
        :param e: Easting coordinate (or longitude in degrees)
        :type e: float
        :return: Linear grid node index
        :rtype: int
        :raises OutOfGridErr: If location is outside grid or does not coincide with a node
        """
        if self.n_min <= n <= self.n_max and self.e_min <= e <= self.e_max:
            c_idx = (e - self.e_min) / self.e_step
            r_idx = (n - self.n_min) / self.n_step

            if c_idx.is_integer() and r_idx.is_integer():
                return int(r_idx) * self.c_count + int(c_idx)
            else:
                raise self._out_of_grid_exc(n, e, self)
        else:
            raise self._out_of_grid_exc(n, e, self)

    @property
    @abc.abstractmethod
    def _get_v_size(self) -> int:
        """Number of components per grid node (1 for elevation, 2 for horizontal shifts)."""
        pass

    def covered_by_grid(self, n: float, e: float) -> bool:
        """Checks if the 4x4 interpolation subgrid centered at (n, e) is fully within grid bounds.

        :param n: Northing coordinate (or latitude in degrees)
        :type n: float
        :param e: Easting coordinate (or longitude in degrees)
        :type e: float
        :return: True if full 4x4 subgrid is covered by grid bounds, False otherwise
        :rtype: bool
        """
        n_sg_dist = (self.sg_size / 2 - 1) * self.n_step
        e_sg_dist = (self.sg_size / 2 - 1) * self.e_step

        return (
            n - n_sg_dist > self.n_min
            and n + n_sg_dist < self.n_max
            and e - e_sg_dist > self.e_min
            and e + e_sg_dist < self.e_max
        )

    def _sgrid_idxs(self, n: float, e: float) -> Tuple[int, ...]:
        """Returns the 16 node indices of the 4x4 subgrid required for bicubic interpolation.

        :param n: Northing coordinate
        :type n: float
        :param e: Easting coordinate
        :type e: float
        :return: 16-tuple of node linear indices
        :rtype: Tuple[int, ...]
        """
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
            (r_idx + 3) * c_count + c_idx + 3,
        )

    @abc.abstractmethod
    def _grid_vs_at_idxs(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Returns grid component values from SPG memory for given indices."""
        pass

    def _grid_vs_no_data(self, values: Iterable[Tuple[float, ...]]) -> bool:
        """Checks if any node value contains No Data (NaN or 999).

        :param values: Iterable of node component value tuples
        :type values: Iterable[Tuple[float, ...]]
        :return: True if No Data is detected, False otherwise
        :rtype: bool
        """
        for v in values:
            if any(math.isnan(x) or x == self.no_data for x in v):
                return True
        return False

    def _reduce_to_unity(self, n: float, e: float) -> Tuple[float, float]:
        """Reduces input coordinates to normalized unit-cell coordinates in :math:`[0, 1] \\times [0, 1]`.

        :param n: Northing coordinate
        :type n: float
        :param e: Easting coordinate
        :type e: float
        :return: Normalized ``(n_unity, e_unity)`` coordinates
        :rtype: Tuple[float, float]
        """
        return (
            ((n - self.n_min) % self.n_step) / self.n_step,
            ((e - self.e_min) % self.e_step) / self.e_step,
        )

    @abc.abstractmethod
    def get_vs_at_idxs_cached(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Cached retrieval of grid node values."""
        pass

    @functools.lru_cache(maxsize=8192)
    def _init_interp(self, values: Tuple[float, ...]) -> "Grid.BiInterp":
        """Factory method with LRU caching for 16-coefficient BiInterp instances."""
        return self.BiInterp(values)

    def interp(
        self,
        n: Sequence[float],
        e: Sequence[float],
    ) -> Tuple[list[float], ...]:
        """Interpolates grid values at input coordinate sequences using bicubic spline interpolation.

        :param n: Sequence of northing coordinates (or latitudes in degrees)
        :type n: Sequence[float]
        :param e: Sequence of easting coordinates (or longitudes in degrees)
        :type e: Sequence[float]
        :return: Tuple containing interpolated lists (1 list for 1D, 2 lists for 2D)
        :rtype: Tuple[list[float], ...]
        :raises OutOfGridErr: If any point does not have a complete 4x4 subgrid
        :raises NoDataGridErr: If any node within the 4x4 subgrid carries a No Data value
        """
        if self.v_size == 1:
            z_out: list[float] = []
            for n_val, e_val in zip(n, e):
                if not self.covered_by_grid(n_val, e_val):
                    raise self._out_of_grid_exc(n_val, e_val, self)

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
            n_out: list[float] = []
            e_out: list[float] = []
            for n_val, e_val in zip(n, e):
                if not self.covered_by_grid(n_val, e_val):
                    raise self._out_of_grid_exc(n_val, e_val, self)

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

        raise ValueError(f"Unexpected grid component size: {self.v_size}")

    class BiInterp:
        """Bicubic spline polynomial interpolation engine on a 4x4 node subgrid.

        Mathematical Formulation:
            Evaluates the continuous 2D surface within the unit cell
            :math:`(u, v) \\in [0, 1] \\times [0, 1]` via:
            .. math::
                P(u, v) = \\sum_{i=0}^3 \\sum_{j=0}^3 a_{4i+j} \\, v^j \\, u^i

            where :math:`u = n_\\text{unity}` and :math:`v = e_\\text{unity}`.
            The 16 polynomial coefficients :math:`a` are computed from the 16 nodal
            values :math:`g_0, \\dots, g_{15}` using central finite-difference partial
            derivatives along East, North, and cross-derivatives.
        """

        def __init__(self, g: Sequence[float]) -> None:
            """Solves the 16 bicubic coefficients from 16 subgrid node values.

            :param g: Sequence of 16 float values ordered row-wise across 4x4 subgrid
            :type g: Sequence[float]
            """
            df: dict[int, float] = {}
            df[0] = g[5]
            df[1] = g[6]
            df[2] = g[9]
            df[3] = g[10]

            # Derivatives in the East direction and the North direction
            df[4] = (-g[7] + 4 * g[6] - 3 * g[5]) / 2.0
            df[5] = (3 * g[6] - 4 * g[5] + g[4]) / 2.0
            df[6] = (-g[11] + 4 * g[10] - 3 * g[9]) / 2.0
            df[7] = (3 * g[10] - 4 * g[9] + g[8]) / 2.0
            df[8] = (-g[13] + 4 * g[9] - 3 * g[5]) / 2.0
            df[9] = (-g[14] + 4 * g[10] - 3 * g[6]) / 2.0
            df[10] = (3 * g[9] - 4 * g[5] + g[1]) / 2.0
            df[11] = (3 * g[10] - 4 * g[6] + g[2]) / 2.0

            # Equations for the cross derivative
            df[12] = ((g[0] + g[10]) - (g[2] + g[8])) / 4.0
            df[13] = ((g[1] + g[11]) - (g[3] + g[9])) / 4.0
            df[14] = ((g[4] + g[14]) - (g[6] + g[12])) / 4.0
            df[15] = ((g[5] + g[15]) - (g[7] + g[13])) / 4.0

            self.__a: Tuple[float, ...] = (
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
                (
                    9 * df[0]
                    - 9 * df[1]
                    - 9 * df[2]
                    + 9 * df[3]
                    + 6 * df[4]
                    + 3 * df[5]
                    - 6 * df[6]
                    - 3 * df[7]
                    + 6 * df[8]
                    - 6 * df[9]
                    + 3 * df[10]
                    - 3 * df[11]
                    + 4 * df[12]
                    + 2 * df[13]
                    + 2 * df[14]
                    + df[15]
                ),
                (
                    -6 * df[0]
                    + 6 * df[1]
                    + 6 * df[2]
                    - 6 * df[3]
                    - 3 * df[4]
                    - 3 * df[5]
                    + 3 * df[6]
                    + 3 * df[7]
                    - 4 * df[8]
                    + 4 * df[9]
                    - 2 * df[10]
                    + 2 * df[11]
                    - 2 * df[12]
                    - 2 * df[13]
                    - df[14]
                    - df[15]
                ),
                2 * df[0] - 2 * df[2] + df[8] + df[10],
                2 * df[4] - 2 * df[6] + df[12] + df[14],
                (
                    -6 * df[0]
                    + 6 * df[1]
                    + 6 * df[2]
                    - 6 * df[3]
                    - 4 * df[4]
                    - 2 * df[5]
                    + 4 * df[6]
                    + 2 * df[7]
                    - 3 * df[8]
                    + 3 * df[9]
                    - 3 * df[10]
                    + 3 * df[11]
                    - 2 * df[12]
                    - df[13]
                    - 2 * df[14]
                    - df[15]
                ),
                (
                    4 * df[0]
                    - 4 * df[1]
                    - 4 * df[2]
                    + 4 * df[3]
                    + 2 * df[4]
                    + 2 * df[5]
                    - 2 * df[6]
                    - 2 * df[7]
                    + 2 * df[8]
                    - 2 * df[9]
                    + 2 * df[10]
                    - 2 * df[11]
                    + df[12]
                    + df[13]
                    + df[14]
                    + df[15]
                ),
            )

        def interp(self, n_unity: float, e_unity: float) -> float:
            """Interpolates continuous value at unit-cell normalized coordinates.

            :param n_unity: Normalized Northing coordinate in :math:`[0, 1]`
            :type n_unity: float
            :param e_unity: Normalized Easting coordinate in :math:`[0, 1]`
            :type e_unity: float
            :return: Interpolated surface value
            :rtype: float
            """
            r = 0.0
            for i in range(4):
                for j in range(4):
                    r += self.__a[i * 4 + j] * (e_unity**j) * (n_unity**i)
            return r


class Grid1D(Grid):
    """1D Vertical quasi-geoid anomaly grid transformation engine."""

    _out_of_grid_exc = exceptions.OutOfGrid1DErr

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initializes 1D grid and binds transformation worker dynamically based on metadata.

        :param filename: Optional grid filename or path
        :type filename: Optional[str]
        """
        super().__init__(filename)

        # Dynamically bind the transformation method once based on grid metadata
        method = self.reader.interp_vertical
        self.trans: Callable[
            [Sequence[float], Sequence[float], Sequence[float], int],
            Tuple[list[float]],
        ]
        if method == 0:  # INTERP_COLOCATE (Nearest-Neighbor)
            self.trans = self._trans_colocate
        elif method == 2:  # INTERP_BICUBIC (Bicubic Spline)
            self.trans = self._trans_bicubic
        else:
            raise NotImplementedError(f"Unsupported vertical interpolation strategy code: {method}")

    def _set_bounds(self) -> None:
        """Sets bounds for 1D geoid heights from SPG reader."""
        (
            self.e_min,
            self.e_max,
            self.n_min,
            self.n_max,
            self.e_step,
            self.n_step,
        ) = self.reader.height_bounds

    def _grid_vs_at_idxs(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Returns height values for requested linear indices."""
        return {idx: (self.reader.height_flat[idx],) for idx in idxs}

    @property
    def _get_v_size(self) -> int:
        """Returns 1 component (geoid undulation)."""
        return 1

    @functools.lru_cache(maxsize=4096)
    def get_vs_at_idxs_cached(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Cached retrieval of 1D height node values."""
        return self._grid_vs_at_idxs(idxs)

    def _trans_colocate(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: Sequence[float],
        corr_sgn: int,
    ) -> Tuple[list[float]]:
        """Transforms elevations using nearest-neighbor (colocate) node lookup.

        :param n: Sequence of latitudes in degrees
        :type n: Sequence[float]
        :param e: Sequence of longitudes in degrees
        :type e: Sequence[float]
        :param z: Sequence of elevation values in meters
        :type z: Sequence[float]
        :param corr_sgn: Direction sign: ``+1`` for :math:`h = H + N`, ``-1`` for :math:`H = h - N`
        :type corr_sgn: int
        :return: 1-element tuple containing list of transformed elevations
        :rtype: Tuple[list[float]]
        """
        c_count = self.c_count
        r_count = self.r_count
        e_min = self.e_min
        e_max = self.e_max
        n_min = self.n_min
        n_max = self.n_max
        e_step = self.e_step
        n_step = self.n_step
        no_data = self.no_data
        height_flat = self.reader.height_flat

        z_out: list[float] = []
        for n_val, e_val, z_val in zip(n, e, z):
            if not (n_min <= n_val <= n_max and e_min <= e_val <= e_max):
                raise self._out_of_grid_exc(n_val, e_val, self)

            c_idx = int(round((e_val - e_min) / e_step))
            r_idx = int(round((n_val - n_min) / n_step))

            if c_idx < 0:
                c_idx = 0
            elif c_idx >= c_count:
                c_idx = c_count - 1

            if r_idx < 0:
                r_idx = 0
            elif r_idx >= r_count:
                r_idx = r_count - 1

            idx = r_idx * c_count + c_idx
            corr = height_flat[idx]
            if math.isnan(corr) or corr == no_data:
                raise exceptions.NoDataGridErr(n_val, e_val, self)

            z_out.append(z_val + corr_sgn * corr)
        return (z_out,)

    def _trans_bicubic(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: Sequence[float],
        corr_sgn: int,
    ) -> Tuple[list[float]]:
        """Transforms elevations using 4x4 bicubic spline polynomial interpolation.

        :param n: Sequence of latitudes in degrees
        :type n: Sequence[float]
        :param e: Sequence of longitudes in degrees
        :type e: Sequence[float]
        :param z: Sequence of elevation values in meters
        :type z: Sequence[float]
        :param corr_sgn: Direction sign: ``+1`` for :math:`h = H + N`, ``-1`` for :math:`H = h - N`
        :type corr_sgn: int
        :return: 1-element tuple containing list of transformed elevations
        :rtype: Tuple[list[float]]
        """
        corrs = self.interp(n, e)
        return ([z_val + corr_sgn * corr for z_val, corr in zip(z, corrs[0])],)


class Grid2D(Grid):
    """2D Horizontal distortion grid transformation engine."""

    _out_of_grid_exc = exceptions.OutOfGrid2DErr

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initializes 2D grid and bounds.

        :param filename: Optional grid filename or path
        :type filename: Optional[str]
        """
        super().__init__(filename)

    def _set_bounds(self) -> None:
        """Sets bounds for 2D horizontal shifts from SPG reader."""
        (
            self.e_min,
            self.e_max,
            self.n_min,
            self.n_max,
            self.e_step,
            self.n_step,
        ) = self.reader.shift_bounds

    def _grid_vs_at_idxs(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Returns 2D shift values (Easting shift, Northing shift) for requested indices."""
        return {idx: (self.reader.shift_e_flat[idx], self.reader.shift_n_flat[idx]) for idx in idxs}

    @property
    def _get_v_size(self) -> int:
        """Returns 2 components (Easting and Northing shifts)."""
        return 2

    @functools.lru_cache(maxsize=4096)
    def get_vs_at_idxs_cached(self, idxs: Tuple[int, ...]) -> dict[int, Tuple[float, ...]]:
        """Cached retrieval of 2D shift node values."""
        return self._grid_vs_at_idxs(idxs)

    def trans(
        self,
        n: Sequence[float],
        e: Sequence[float],
        corr_sgn: int,
    ) -> Tuple[list[float], list[float]]:
        """Transforms Northing and Easting coordinates by adding/subtracting interpolated shifts.

        :param n: Sequence of northing coordinates in meters
        :type n: Sequence[float]
        :param e: Sequence of easting coordinates in meters
        :type e: Sequence[float]
        :param corr_sgn: Direction sign: ``+1`` or ``-1``
        :type corr_sgn: int
        :return: Tuple of ``(n_out, e_out)`` lists with applied shifts
        :rtype: Tuple[list[float], list[float]]
        """
        corrs = self.interp(n, e)
        return (
            [n_val + corr_sgn * corr for n_val, corr in zip(n, corrs[0])],
            [e_val + corr_sgn * corr for e_val, corr in zip(e, corrs[1])],
        )
