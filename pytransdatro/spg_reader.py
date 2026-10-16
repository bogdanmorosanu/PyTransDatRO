"""Binary SPG Grid Deserializer and Cache.

This module provides the SpgReader singleton class that loads and deserializes
binary `.spg` transformation grids, extracting interpolation parameters,
2D geodetic shifts, 1D geoid heights, and Helmert transformation values into
in-memory pure Python tuples.
"""

import glob
import os
import pickle
from typing import Any, Optional, Tuple

from pytransdatro import utils
from pytransdatro.exceptions import AmbiguousGridError, MissingGridError


class SpgReader:
    """Singleton class to load the .spg binary file exactly once and keep
    its flat tuple data and transformation metadata in memory for pure-Python processing.
    """

    _instance: Optional["SpgReader"] = None

    def __new__(cls, grid_filename: Optional[str] = None) -> "SpgReader":
        """Instantiate or return the existing singleton SpgReader instance.

        :param grid_filename: Optional filename or path of the .spg grid file.
            If None, the grid is auto-discovered from pytransdatro/grids/.
        :type grid_filename: Optional[str]
        :return: Singleton SpgReader instance
        :rtype: SpgReader
        """
        if cls._instance is None:
            cls._instance = super(SpgReader, cls).__new__(cls)
            cls._instance._load(grid_filename)
        return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset the singleton instance (useful for unit testing)."""
        cls._instance = None

    def _resolve_grid_file(self, grid_filename: Optional[str]) -> str:
        """Resolves the absolute path to the .spg grid file.

        :param grid_filename: Optional filename or explicit absolute path
        :type grid_filename: Optional[str]
        :return: Absolute path to the verified grid file
        :rtype: str
        :raises MissingGridError: If no grid file is found
        :raises AmbiguousGridError: If multiple grid files exist without explicit selection
        """
        grids_dir = os.path.join(os.path.dirname(__file__), "grids")

        if grid_filename is not None:
            if os.path.isabs(grid_filename):
                grid_path = grid_filename
            else:
                grid_path = os.path.join(grids_dir, grid_filename)
            if not os.path.isfile(grid_path):
                raise MissingGridError(f"Grid file '{grid_path}' not found.")
            return grid_path

        spg_files = glob.glob(os.path.join(grids_dir, "*.spg"))
        if len(spg_files) == 0:
            raise MissingGridError(f"No .spg grid file found in grids directory: '{grids_dir}'.")
        elif len(spg_files) > 1:
            file_names = [os.path.basename(f) for f in spg_files]
            raise AmbiguousGridError(
                f"Multiple .spg grid files found in '{grids_dir}': {file_names}. "
                "Please ensure exactly one grid file is present or explicitly specify which grid to use."
            )
        return spg_files[0]

    def _load(self, grid_filename: Optional[str]) -> None:
        """Deserializes the SPG file and populates transformation attributes.

        :param grid_filename: Optional grid filename or path
        :type grid_filename: Optional[str]
        """
        grid_path = self._resolve_grid_file(grid_filename)
        self.grid_path: str = grid_path
        self.grid_filename: str = os.path.basename(grid_path)

        with open(grid_path, "rb") as f:
            spg_data: dict[str, Any] = pickle.load(f)

        # 1. Metadata and Parameters
        self.params: dict[str, Any] = spg_data.get("params", {})
        self.metadata: dict[str, Any] = spg_data.get("metadata", {})

        # Interpolation strategy
        interp_cfg: dict[str, Any] = self.params["interpolation"]
        self.interp_horizontal: int = interp_cfg["horizontal"]
        self.interp_vertical: int = interp_cfg["vertical"]

        # Helmert parameters (st70 <-> os intermediate system)
        helmert_cfg: dict[str, Any] = self.params["helmert"]
        self.helmert_st70_os: dict[str, Any] = helmert_cfg["st70_os"]
        self.helmert_os_st70: dict[str, Any] = helmert_cfg["os_st70"]

        # Precompute Helmert parameters for st70 -> os (standard Stereo70 to StereoGRS80 direction)
        self.helmert_tn: float = float(self.helmert_st70_os["tN"])
        self.helmert_te: float = float(self.helmert_st70_os["tE"])
        self.helmert_ppm: float = float(self.helmert_st70_os["dm"])

        # Rz in SPG is in arcseconds: convert to sexagesimal string '0 0 Rz' or radians directly via utils
        rz_sec = float(self.helmert_st70_os["Rz"])
        self.helmert_rot_rad: float = utils.sexa_to_rad(f"0 0 {abs(rz_sec)}")
        if rz_sec < 0:
            self.helmert_rot_rad = -self.helmert_rot_rad

        # 2. Geodetic shifts (2D horizontal)
        shifts_data = spg_data["grids"]["geodetic_shifts"]
        s_meta = shifts_data["metadata"]

        # Determine actual array boundaries based on numpy shape to prevent index overflow
        s_grid = shifts_data["grid"]
        c_count, r_count = s_grid[0].shape[1], s_grid[0].shape[0]
        maxe = s_meta["mine"] + (c_count - 1) * s_meta["stepe"]
        maxn = s_meta["minn"] + (r_count - 1) * s_meta["stepn"]
        self.shift_bounds: Tuple[float, float, float, float, float, float] = (
            s_meta["mine"],
            maxe,
            s_meta["minn"],
            maxn,
            s_meta["stepe"],
            s_meta["stepn"],
        )

        # Convert numpy array to flat pure python tuple for processing
        self.shift_e_flat: Tuple[float, ...] = tuple(float(x) for x in s_grid[0].flatten())
        self.shift_n_flat: Tuple[float, ...] = tuple(float(x) for x in s_grid[1].flatten())

        # 3. Geoid heights (1D vertical)
        heights_data = spg_data["grids"]["geoid_heights"]
        h_meta = heights_data["metadata"]
        h_grid = heights_data["grid"][0]

        # Determine actual array boundaries for height grid
        hc_count, hr_count = h_grid.shape[1], h_grid.shape[0]
        maxla = h_meta["minla"] + (hc_count - 1) * h_meta["stepla"]
        maxphi = h_meta["minphi"] + (hr_count - 1) * h_meta["stepphi"]
        self.height_bounds: Tuple[float, float, float, float, float, float] = (
            h_meta["minla"],
            maxla,
            h_meta["minphi"],
            maxphi,
            h_meta["stepla"],
            h_meta["stepphi"],
        )
        self.height_flat: Tuple[float, ...] = tuple(float(x) for x in h_grid.flatten())
