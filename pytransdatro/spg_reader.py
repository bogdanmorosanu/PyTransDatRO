import os
import glob
import pickle
import numpy as np
from pytransdatro.exceptions import MissingGridError, AmbiguousGridError
from pytransdatro import utils

class SpgReader:
    """Singleton class to load the .spg binary file exactly once and keep 
    its flat tuple data and transformation metadata in memory for pure-python processing."""
    
    _instance = None

    def __new__(cls, grid_filename=None):
        if cls._instance is None:
            cls._instance = super(SpgReader, cls).__new__(cls)
            cls._instance._load(grid_filename)
        return cls._instance

    @classmethod
    def reset_instance(cls):
        """Reset the singleton instance (useful for unit testing)."""
        cls._instance = None

    def _resolve_grid_file(self, grid_filename):
        grids_dir = os.path.join(os.path.dirname(__file__), 'grids')
        
        if grid_filename is not None:
            if os.path.isabs(grid_filename):
                grid_path = grid_filename
            else:
                grid_path = os.path.join(grids_dir, grid_filename)
            if not os.path.isfile(grid_path):
                raise MissingGridError(f"Grid file '{grid_path}' not found.")
            return grid_path
            
        spg_files = glob.glob(os.path.join(grids_dir, '*.spg'))
        if len(spg_files) == 0:
            raise MissingGridError(f"No .spg grid file found in grids directory: '{grids_dir}'.")
        elif len(spg_files) > 1:
            file_names = [os.path.basename(f) for f in spg_files]
            raise AmbiguousGridError(
                f"Multiple .spg grid files found in '{grids_dir}': {file_names}. "
                f"Please ensure exactly one grid file is present or explicitly specify which grid to use."
            )
        return spg_files[0]

    def _load(self, grid_filename):
        grid_path = self._resolve_grid_file(grid_filename)
        self.grid_path = grid_path
        self.grid_filename = os.path.basename(grid_path)
        
        with open(grid_path, 'rb') as f:
            spg_data = pickle.load(f)
        
        # 1. Metadata and Parameters
        self.params = spg_data.get('params', {})
        self.metadata = spg_data.get('metadata', {})
        
        # Interpolation strategy
        interp_cfg = self.params['interpolation']
        self.interp_horizontal = interp_cfg['horizontal']
        self.interp_vertical = interp_cfg['vertical']
        
        # Helmert parameters (st70 <-> os intermediate system)
        helmert_cfg = self.params['helmert']
        self.helmert_st70_os = helmert_cfg['st70_os']
        self.helmert_os_st70 = helmert_cfg['os_st70']
        
        # Precompute Helmert parameters for st70 -> os (standard Stereo70 to StereoGRS80 direction)
        self.helmert_tn = float(self.helmert_st70_os['tN'])
        self.helmert_te = float(self.helmert_st70_os['tE'])
        self.helmert_ppm = float(self.helmert_st70_os['dm'])
        
        # Rz in SPG is in arcseconds: convert to sexagesimal string '0 0 Rz' or radians directly via utils
        rz_sec = float(self.helmert_st70_os['Rz'])
        # utils.sexa_to_rad accepts '0 0 <seconds>'
        self.helmert_rot_rad = utils.sexa_to_rad(f"0 0 {abs(rz_sec)}")
        if rz_sec < 0:
            self.helmert_rot_rad = -self.helmert_rot_rad

        # 2. Geodetic shifts (2D horizontal)
        shifts_data = spg_data['grids']['geodetic_shifts']
        s_meta = shifts_data['metadata']
        
        # Determine actual array boundaries based on numpy shape to prevent index overflow
        s_grid = shifts_data['grid']
        c_count, r_count = s_grid[0].shape[1], s_grid[0].shape[0]
        maxe = s_meta['mine'] + (c_count - 1) * s_meta['stepe']
        maxn = s_meta['minn'] + (r_count - 1) * s_meta['stepn']
        self.shift_bounds = (s_meta['mine'], maxe, s_meta['minn'], maxn, s_meta['stepe'], s_meta['stepn'])
        
        # Convert numpy array to flat pure python tuple for processing
        self.shift_e_flat = tuple(float(x) for x in s_grid[0].flatten())
        self.shift_n_flat = tuple(float(x) for x in s_grid[1].flatten())
        
        # 3. Geoid heights (1D vertical)
        heights_data = spg_data['grids']['geoid_heights']
        h_meta = heights_data['metadata']
        h_grid = heights_data['grid'][0]
        
        # Determine actual array boundaries for height grid
        hc_count, hr_count = h_grid.shape[1], h_grid.shape[0]
        maxla = h_meta['minla'] + (hc_count - 1) * h_meta['stepla']
        maxphi = h_meta['minphi'] + (hr_count - 1) * h_meta['stepphi']
        self.height_bounds = (h_meta['minla'], maxla, h_meta['minphi'], maxphi, h_meta['stepla'], h_meta['stepphi'])
        self.height_flat = tuple(float(x) for x in h_grid.flatten())
