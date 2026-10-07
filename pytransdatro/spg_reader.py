import os
import pickle
import numpy as np

class SpgReader:
    """Singleton class to load the .spg binary file exactly once and keep 
    its flat tuple data in memory for pure-python interpolation processing."""
    
    _instance = None

    def __new__(cls, grid_filename='rom_grid3d_25.09.spg'):
        if cls._instance is None:
            cls._instance = super(SpgReader, cls).__new__(cls)
            cls._instance._load(grid_filename)
        return cls._instance

    def _load(self, grid_filename):
        grid_path = os.path.join(os.path.dirname(__file__), 'grids', grid_filename)
        with open(grid_path, 'rb') as f:
            spg_data = pickle.load(f)
        
        # Geodetic shifts (2D)
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
        
        # Geoid heights (1D)
        heights_data = spg_data['grids']['geoid_heights']
        h_meta = heights_data['metadata']
        h_grid = heights_data['grid'][0]
        
        # Determine actual array boundaries for height grid
        hc_count, hr_count = h_grid.shape[1], h_grid.shape[0]
        maxla = h_meta['minla'] + (hc_count - 1) * h_meta['stepla']
        maxphi = h_meta['minphi'] + (hr_count - 1) * h_meta['stepphi']
        self.height_bounds = (h_meta['minla'], maxla, h_meta['minphi'], maxphi, h_meta['stepla'], h_meta['stepphi'])
        self.height_flat = tuple(float(x) for x in h_grid.flatten())
