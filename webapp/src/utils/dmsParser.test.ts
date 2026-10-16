import { describe, it, expect } from 'vitest';
import {
  degToRad,
  radToDeg,
  parseAngleToDeg,
  parseAngleToRad,
  formatDegToDMS,
  DEG_TO_RAD,
  RAD_TO_DEG,
} from './dmsParser';

describe('dmsParser', () => {
  describe('constants & conversions', () => {
    it('correctly converts degrees to radians and vice versa', () => {
      expect(degToRad(180)).toBeCloseTo(Math.PI, 10);
      expect(degToRad(90)).toBeCloseTo(Math.PI / 2, 10);
      expect(degToRad(0)).toBe(0);

      expect(radToDeg(Math.PI)).toBeCloseTo(180, 10);
      expect(radToDeg(Math.PI / 2)).toBeCloseTo(90, 10);

      const angle = 45.54321;
      expect(radToDeg(degToRad(angle))).toBeCloseTo(angle, 10);
    });

    it('exposes accurate DEG_TO_RAD and RAD_TO_DEG factors', () => {
      expect(DEG_TO_RAD * 180).toBeCloseTo(Math.PI, 12);
      expect(RAD_TO_DEG * Math.PI).toBeCloseTo(180, 12);
    });
  });

  describe('parseAngleToDeg', () => {
    it('parses standard decimal degree strings', () => {
      expect(parseAngleToDeg('45.999718')).toBeCloseTo(45.999718, 6);
      expect(parseAngleToDeg('-45.999718')).toBeCloseTo(-45.999718, 6);
      expect(parseAngleToDeg('+45.0')).toBeCloseTo(45.0, 6);
      expect(parseAngleToDeg('45')).toBe(45);
      expect(parseAngleToDeg('  26.1025   ')).toBeCloseTo(26.1025, 6);
    });

    it('parses DMS with degrees, minutes, and seconds symbols', () => {
      // 45 + 59/60 + 58.99/3600 = 45.999719444...
      const expected = 45 + 59 / 60 + 58.99 / 3600;
      expect(parseAngleToDeg("45° 59' 58.99\"")).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg("45°59'58.99\"")).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg("45° 30' 00\"")).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg("45° 30'")).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45°')).toBe(45);
    });

    it('parses DMS with d/m/s letters in both lowercase and uppercase', () => {
      const expected = 45 + 59 / 60 + 58.99 / 3600;
      expect(parseAngleToDeg('45d 59m 58.99s')).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg('45D 59M 58.99S')).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg('45d59m58.99s')).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg('45D30M')).toBeCloseTo(45.5, 8);
    });

    it('parses whitespace-separated DMS strings', () => {
      const expected = 45 + 59 / 60 + 58.99 / 3600;
      expect(parseAngleToDeg('45 59 58.99')).toBeCloseTo(expected, 8);
      expect(parseAngleToDeg('45 30 00')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('45 30')).toBeCloseTo(45.5, 8);
    });

    it('handles cardinal directions (N, S, E, W)', () => {
      expect(parseAngleToDeg("45° 30' 00\" N")).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg("45° 30' 00\" S")).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg("26° 00' 00\" E")).toBeCloseTo(26.0, 8);
      expect(parseAngleToDeg("26° 00' 00\" W")).toBeCloseTo(-26.0, 8);
      expect(parseAngleToDeg('N 45.5')).toBeCloseTo(45.5, 8);
      expect(parseAngleToDeg('S 45.5')).toBeCloseTo(-45.5, 8);
      expect(parseAngleToDeg('26.0 E')).toBeCloseTo(26.0, 8);
      expect(parseAngleToDeg('26.0 W')).toBeCloseTo(-26.0, 8);
    });

    it('throws errors on invalid inputs', () => {
      expect(() => parseAngleToDeg('')).toThrow('Empty coordinate input');
      expect(() => parseAngleToDeg('   ')).toThrow('Empty coordinate input');
      expect(() => parseAngleToDeg('+')).toThrow('Invalid coordinate input');
      expect(() => parseAngleToDeg('-')).toThrow('Invalid coordinate input');
      expect(() => parseAngleToDeg('N')).toThrow('Invalid coordinate input');
      expect(() => parseAngleToDeg('foo')).toThrow();
      expect(() => parseAngleToDeg("45° 60' 00\"")).toThrow('Minutes must be between 0 and 59');
      expect(() => parseAngleToDeg("45° 30' 60\"")).toThrow('Seconds must be between 0 and 59.999');
      expect(() => parseAngleToDeg('45 30 15 10')).toThrow('Invalid DMS format');
    });
  });

  describe('parseAngleToRad', () => {
    it('parses DMS and decimal strings directly to radians', () => {
      expect(parseAngleToRad("180° 00' 00\"")).toBeCloseTo(Math.PI, 8);
      expect(parseAngleToRad("90° 00' 00\"")).toBeCloseTo(Math.PI / 2, 8);
      expect(parseAngleToRad('180')).toBeCloseTo(Math.PI, 8);
    });
  });

  describe('formatDegToDMS', () => {
    it('formats decimal degrees to standard DMS representation', () => {
      expect(formatDegToDMS(45.5)).toBe('45°30\'00.00000"');
      expect(formatDegToDMS(-45.5)).toBe('-45°30\'00.00000"');
    });

    it('respects requested precision', () => {
      expect(formatDegToDMS(45.5, 2)).toBe('45°30\'00.00"');
      expect(formatDegToDMS(45.5, 0)).toBe('45°30\'00"');
      expect(formatDegToDMS(45.999719444, 4)).toBe('45°59\'58.9900"');
    });

    it('handles seconds rollover cleanly', () => {
      // 45 deg, 59 min, 59.9999 sec with 2 decimals rounds to 60.00 sec -> 46 deg 00 min 00.00 sec
      const angleNearRollover = 45 + 59 / 60 + 59.9999 / 3600;
      expect(formatDegToDMS(angleNearRollover, 2)).toBe('46°00\'00.00"');
    });

    it('throws on NaN input', () => {
      expect(() => formatDegToDMS(NaN)).toThrow('Cannot format NaN to DMS');
    });
  });
});
