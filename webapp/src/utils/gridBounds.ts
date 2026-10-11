/**
 * Official Romania TransDatRO SPG Grid Envelope Constants
 * Derived from the union of:
 * 1. Planimetric 2D Grid (etrs89_krasovschi42 in Stereo 70)
 * 2. Vertical 1D Quasigeoid Grid (romhybqgeoid_cs1 in ETRS89)
 * Extracted via tools/extract_spg_metadata.py from pytransdatro/grids/rom_grid3d_25.09.spg
 */

export const STEREO70_GRID_BOUNDS = {
  minN: 213600.0,
  maxN: 785640.0,
  minE: 101260.0,
  maxE: 890790.0,
};

export const ETRS89_GRID_BOUNDS = {
  minLat: 43.30,
  maxLat: 48.58,
  minLon: 19.70,
  maxLon: 30.30,
};

/**
 * Coarse pre-flight envelope check for Stereo 70 coordinates (meters).
 */
export function isWithinStereo70Grid(northing: number, easting: number): boolean {
  if (isNaN(northing) || isNaN(easting)) return false;
  return (
    northing >= STEREO70_GRID_BOUNDS.minN &&
    northing <= STEREO70_GRID_BOUNDS.maxN &&
    easting >= STEREO70_GRID_BOUNDS.minE &&
    easting <= STEREO70_GRID_BOUNDS.maxE
  );
}

/**
 * Coarse pre-flight envelope check for ETRS89 coordinates (decimal degrees).
 */
export function isWithinETRS89Grid(lat: number, lon: number): boolean {
  if (isNaN(lat) || isNaN(lon)) return false;
  return (
    lat >= ETRS89_GRID_BOUNDS.minLat &&
    lat <= ETRS89_GRID_BOUNDS.maxLat &&
    lon >= ETRS89_GRID_BOUNDS.minLon &&
    lon <= ETRS89_GRID_BOUNDS.maxLon
  );
}
