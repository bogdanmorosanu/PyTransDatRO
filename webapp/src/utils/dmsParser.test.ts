import { describe, it, expect } from 'vitest';
import {
  degToRad,
  radToDeg,
  parseAngleToDeg,
  parseAngleToRad,
  formatDegToDMS,
  parseAngularTokenToDeg,
} from './dmsParser';

describe('dmsParser', () => {
  describe('degToRad and radToDeg', () => {
    it('converts degrees to radians and back accurately', () => {
      const deg = 45.0;
      const rad = degToRad(deg);
      expect(rad).toBeCloseTo(Math.PI / 4, 10);
      expect(radToDeg(rad)).toBeCloseTo(deg, 10);
    });
  });

  describe('parseAngleToDeg', () => {
    it('parses pure decimal degrees', () => {
      expect(parseAngleToDeg('45.5')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('-45.5')).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg('+45.5')).toBeCloseTo(45.5, 8);
    });

    it('parses DMS with standard symbols', () => {
      expect(parseAngleToDeg('45° 30\' 00"')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45°30\'00"')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45d 30m 00s')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45 30 00')).toBeCloseTo(45.5, 8);
    });

    it('parses cardinal directions correctly', () => {
      expect(parseAngleToDeg('45.5 N')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45.5 S')).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg('S 45.5')).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg('45.5 W')).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg('45.5 E')).toBeCloseTo(45.5, 8);
    });

    it('throws on invalid or empty inputs', () => {
      expect(() => parseAngleToDeg('')).toThrow();
      expect(() => parseAngleToDeg('abc')).toThrow();
      expect(() => parseAngleToDeg('45° 70\' 00"')).toThrow();
    });
  });

  describe('parseAngleToRad', () => {
    it('parses DMS directly to radians', () => {
      expect(parseAngleToRad('180°')).toBeCloseTo(Math.PI, 10);
      expect(parseAngleToRad('90°')).toBeCloseTo(Math.PI / 2, 10);
    });
  });

  describe('formatDegToDMS', () => {
    it('formats degrees to standard DMS strings', () => {
      expect(formatDegToDMS(45.5)).toBe('45°30\'00.00000"');
      expect(formatDegToDMS(-45.5)).toBe('-45°30\'00.00000"');
    });

    it('respects requested precision', () => {
      expect(formatDegToDMS(45.5, 2)).toBe('45°30\'00.00"');
      expect(formatDegToDMS(45.5, 0)).toBe('45°30\'00"');
      expect(formatDegToDMS(45.999719444, 4)).toBe('45°59\'58.9900"');
    });

    it('handles seconds rollover cleanly', () => {
      const angleNearRollover = 45 + 59 / 60 + 59.9999 / 3600;
      expect(formatDegToDMS(angleNearRollover, 2)).toBe('46°00\'00.00"');
    });

    it('throws on NaN input', () => {
      expect(() => formatDegToDMS(NaN)).toThrow('Cannot format NaN to DMS');
    });
  });

  describe('parseAngularTokenToDeg', () => {
    it('parses DMS and decimal tokens in dms / degrees mode', () => {
      expect(parseAngularTokenToDeg('45.5', 'degrees')).toBeCloseTo(45.5, 8);
      expect(parseAngularTokenToDeg('45° 30\' 00"', 'dms')).toBeCloseTo(45.5, 8);
    });

    it('parses radians tokens with and without rad suffix', () => {
      const piOverFour = Math.PI / 4;
      expect(parseAngularTokenToDeg(String(piOverFour), 'radians')).toBeCloseTo(45.0, 8);
      expect(parseAngularTokenToDeg(`${piOverFour} rad`, 'radians')).toBeCloseTo(45.0, 8);
    });

    it('falls back to DMS parser if unit is radians but input contains DMS symbols', () => {
      expect(parseAngularTokenToDeg('45° 30\' 00"', 'radians')).toBeCloseTo(45.5, 8);
      expect(parseAngularTokenToDeg("45°30'00\"", 'radians')).toBeCloseTo(45.5, 8);
    });

    it('returns NaN on empty, null, or invalid input', () => {
      expect(isNaN(parseAngularTokenToDeg(''))).toBe(true);
      expect(isNaN(parseAngularTokenToDeg(null))).toBe(true);
      expect(isNaN(parseAngularTokenToDeg(undefined))).toBe(true);
      expect(isNaN(parseAngularTokenToDeg('not-a-number', 'degrees'))).toBe(true);
      expect(isNaN(parseAngularTokenToDeg('not-a-number', 'radians'))).toBe(true);
    });
  });
});
