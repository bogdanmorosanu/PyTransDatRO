import { describe, it, expect } from 'vitest';
import {
  formatAngularCoordinate,
  formatPlanarCoordinate,
  formatCoordinatePair,
  computeQuasigeoidZeta,
  formatQuasigeoidZeta,
} from './coordinateFormat';

describe('coordinateFormat', () => {
  describe('formatAngularCoordinate', () => {
    it('formats DMS, radians, and decimal degrees with canonical precisions', () => {
      const lat = 45.5;
      expect(formatAngularCoordinate(lat, 'dms')).toBe('45°30\'00.00000"');
      expect(formatAngularCoordinate(lat, 'degrees')).toBe('45.50000000');
      expect(formatAngularCoordinate(lat, 'radians')).toBe('0.7941248097');
    });

    it('returns empty string on NaN input', () => {
      expect(formatAngularCoordinate(NaN, 'degrees')).toBe('');
    });
  });

  describe('formatPlanarCoordinate', () => {
    it('formats planar coordinates with 3 decimals by default', () => {
      expect(formatPlanarCoordinate(500000.1234)).toBe('500000.123');
      expect(formatPlanarCoordinate(500000)).toBe('500000.000');
    });

    it('returns empty string on NaN input', () => {
      expect(formatPlanarCoordinate(NaN)).toBe('');
    });
  });

  describe('formatCoordinatePair', () => {
    it('formats geographic coordinates respecting NE and EN order', () => {
      const lat = 46.0;
      const lon = 25.0;

      const ne = formatCoordinatePair(lat, lon, 'degrees', 'NE', true);
      expect(ne.col1).toBe('46.00000000');
      expect(ne.col2).toBe('25.00000000');

      const en = formatCoordinatePair(lat, lon, 'degrees', 'EN', true);
      expect(en.col1).toBe('25.00000000');
      expect(en.col2).toBe('46.00000000');
    });

    it('formats planar coordinates respecting NE and EN order', () => {
      const n = 500000.0;
      const e = 600000.0;

      const ne = formatCoordinatePair(n, e, 'degrees', 'NE', false);
      expect(ne.col1).toBe('500000.000');
      expect(ne.col2).toBe('600000.000');

      const en = formatCoordinatePair(n, e, 'degrees', 'EN', false);
      expect(en.col1).toBe('600000.000');
      expect(en.col2).toBe('500000.000');
    });
  });

  describe('computeQuasigeoidZeta', () => {
    it('computes zeta = h - H accurately', () => {
      const h = 135.25;
      const H = 100.0;
      expect(computeQuasigeoidZeta(h, H)).toBeCloseTo(35.25, 4);
    });

    it('returns null if either height is missing or NaN', () => {
      expect(computeQuasigeoidZeta(undefined, 100)).toBeNull();
      expect(computeQuasigeoidZeta(100, undefined)).toBeNull();
      expect(computeQuasigeoidZeta(null, 100)).toBeNull();
      expect(computeQuasigeoidZeta(NaN, 100)).toBeNull();
    });
  });

  describe('formatQuasigeoidZeta', () => {
    it('formats positive zeta with + sign and 3 decimals', () => {
      expect(formatQuasigeoidZeta(0.3456)).toBe('+0.346');
      expect(formatQuasigeoidZeta(0)).toBe('+0.000');
    });

    it('formats negative zeta with - sign and 3 decimals', () => {
      expect(formatQuasigeoidZeta(-0.1234)).toBe('-0.123');
    });

    it('returns empty string for null, undefined, or NaN', () => {
      expect(formatQuasigeoidZeta(null)).toBe('');
      expect(formatQuasigeoidZeta(undefined)).toBe('');
      expect(formatQuasigeoidZeta(NaN)).toBe('');
    });
  });
});
