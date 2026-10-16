import { describe, it, expect } from 'vitest';
import {
  STEREO70_GRID_BOUNDS,
  ETRS89_GRID_BOUNDS,
  isWithinStereo70Grid,
  isWithinETRS89Grid,
} from './gridBounds';

describe('gridBounds', () => {
  describe('constants', () => {
    it('defines valid envelope ranges for Stereo 70', () => {
      expect(STEREO70_GRID_BOUNDS.minN).toBe(213600.0);
      expect(STEREO70_GRID_BOUNDS.maxN).toBe(785640.0);
      expect(STEREO70_GRID_BOUNDS.minE).toBe(101260.0);
      expect(STEREO70_GRID_BOUNDS.maxE).toBe(890790.0);
      expect(STEREO70_GRID_BOUNDS.minN).toBeLessThan(STEREO70_GRID_BOUNDS.maxN);
      expect(STEREO70_GRID_BOUNDS.minE).toBeLessThan(STEREO70_GRID_BOUNDS.maxE);
    });

    it('defines valid geographic envelope ranges for ETRS89', () => {
      expect(ETRS89_GRID_BOUNDS.minLat).toBe(43.30);
      expect(ETRS89_GRID_BOUNDS.maxLat).toBe(48.58);
      expect(ETRS89_GRID_BOUNDS.minLon).toBe(19.70);
      expect(ETRS89_GRID_BOUNDS.maxLon).toBe(30.30);
      expect(ETRS89_GRID_BOUNDS.minLat).toBeLessThan(ETRS89_GRID_BOUNDS.maxLat);
      expect(ETRS89_GRID_BOUNDS.minLon).toBeLessThan(ETRS89_GRID_BOUNDS.maxLon);
    });
  });

  describe('isWithinStereo70Grid', () => {
    it('returns true for known Romanian city coordinates in Stereo 70', () => {
      // Bucharest: ~324000 N, ~588000 E
      expect(isWithinStereo70Grid(324000, 588000)).toBe(true);
      // Cluj-Napoca: ~585000 N, ~392000 E
      expect(isWithinStereo70Grid(585000, 392000)).toBe(true);
      // Timișoara: ~480000 N, ~206000 E
      expect(isWithinStereo70Grid(480000, 206000)).toBe(true);
      // Constanța: ~300000 N, ~780000 E
      expect(isWithinStereo70Grid(300000, 780000)).toBe(true);
      // Iași: ~640000 N, ~690000 E
      expect(isWithinStereo70Grid(640000, 690000)).toBe(true);
    });

    it('returns true on exact boundary threshold limits', () => {
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.minN, STEREO70_GRID_BOUNDS.minE)).toBe(true);
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.maxN, STEREO70_GRID_BOUNDS.maxE)).toBe(true);
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.minN, STEREO70_GRID_BOUNDS.maxE)).toBe(true);
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.maxN, STEREO70_GRID_BOUNDS.minE)).toBe(true);
    });

    it('returns false for points outside Stereo 70 grid envelope', () => {
      // Just beyond borders
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.minN - 0.1, 500000)).toBe(false);
      expect(isWithinStereo70Grid(STEREO70_GRID_BOUNDS.maxN + 0.1, 500000)).toBe(false);
      expect(isWithinStereo70Grid(500000, STEREO70_GRID_BOUNDS.minE - 0.1)).toBe(false);
      expect(isWithinStereo70Grid(500000, STEREO70_GRID_BOUNDS.maxE + 0.1)).toBe(false);

      // Negative values or far away coordinates
      expect(isWithinStereo70Grid(-1000, 500000)).toBe(false);
      expect(isWithinStereo70Grid(0, 0)).toBe(false);
      expect(isWithinStereo70Grid(1000000, 1000000)).toBe(false);
    });

    it('handles NaN and non-numeric inputs gracefully', () => {
      expect(isWithinStereo70Grid(NaN, 500000)).toBe(false);
      expect(isWithinStereo70Grid(500000, NaN)).toBe(false);
      expect(isWithinStereo70Grid(NaN, NaN)).toBe(false);
    });
  });

  describe('isWithinETRS89Grid', () => {
    it('returns true for known Romanian locations in ETRS89', () => {
      // Bucharest: ~44.43° N, ~26.10° E
      expect(isWithinETRS89Grid(44.4323, 26.1063)).toBe(true);
      // Cluj-Napoca: ~46.77° N, ~23.60° E
      expect(isWithinETRS89Grid(46.7712, 23.6236)).toBe(true);
      // Timișoara: ~45.75° N, ~21.21° E
      expect(isWithinETRS89Grid(45.7489, 21.2087)).toBe(true);
      // Constanța: ~44.18° N, ~28.65° E
      expect(isWithinETRS89Grid(44.1792, 28.6498)).toBe(true);
      // Iași: ~47.16° N, ~27.60° E
      expect(isWithinETRS89Grid(47.1585, 27.6014)).toBe(true);
    });

    it('returns true on exact boundary threshold limits', () => {
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.minLat, ETRS89_GRID_BOUNDS.minLon)).toBe(true);
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.maxLat, ETRS89_GRID_BOUNDS.maxLon)).toBe(true);
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.minLat, ETRS89_GRID_BOUNDS.maxLon)).toBe(true);
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.maxLat, ETRS89_GRID_BOUNDS.minLon)).toBe(true);
    });

    it('returns false for locations outside Romania', () => {
      // Paris, France
      expect(isWithinETRS89Grid(48.8566, 2.3522)).toBe(false);
      // London, UK
      expect(isWithinETRS89Grid(51.5074, -0.1278)).toBe(false);
      // Equator / Null Island
      expect(isWithinETRS89Grid(0, 0)).toBe(false);

      // Just beyond envelope
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.minLat - 0.01, 25.0)).toBe(false);
      expect(isWithinETRS89Grid(ETRS89_GRID_BOUNDS.maxLat + 0.01, 25.0)).toBe(false);
      expect(isWithinETRS89Grid(45.0, ETRS89_GRID_BOUNDS.minLon - 0.01)).toBe(false);
      expect(isWithinETRS89Grid(45.0, ETRS89_GRID_BOUNDS.maxLon + 0.01)).toBe(false);
    });

    it('handles NaN gracefully', () => {
      expect(isWithinETRS89Grid(NaN, 25.0)).toBe(false);
      expect(isWithinETRS89Grid(45.0, NaN)).toBe(false);
      expect(isWithinETRS89Grid(NaN, NaN)).toBe(false);
    });
  });
});
