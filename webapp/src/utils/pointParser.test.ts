import { describe, it, expect } from 'vitest';
import {
  isValidDecimal,
  parseStrictDecimal,
  detectDelimiter,
  parseInputLines,
} from './pointParser';

describe('pointParser', () => {
  describe('isValidDecimal', () => {
    it('returns true for valid decimal strings', () => {
      expect(isValidDecimal('123')).toBe(true);
      expect(isValidDecimal('123.456')).toBe(true);
      expect(isValidDecimal('-123.456')).toBe(true);
      expect(isValidDecimal('+123.456')).toBe(true);
      expect(isValidDecimal('.5')).toBe(true);
      expect(isValidDecimal('0.0')).toBe(true);
      expect(isValidDecimal('  42.5  ')).toBe(true);
    });

    it('returns false for invalid decimal strings', () => {
      expect(isValidDecimal('')).toBe(false);
      expect(isValidDecimal('   ')).toBe(false);
      expect(isValidDecimal('abc')).toBe(false);
      expect(isValidDecimal('12.34.56')).toBe(false);
      expect(isValidDecimal('--5')).toBe(false);
      expect(isValidDecimal('+-5')).toBe(false);
      expect(isValidDecimal('NaN')).toBe(false);
      expect(isValidDecimal('Infinity')).toBe(false);
    });
  });

  describe('parseStrictDecimal', () => {
    it('converts valid strings to numbers', () => {
      expect(parseStrictDecimal('123.45')).toBe(123.45);
      expect(parseStrictDecimal('-500.25')).toBe(-500.25);
      expect(parseStrictDecimal('  0  ')).toBe(0);
    });

    it('throws error on non-decimal inputs', () => {
      expect(() => parseStrictDecimal('invalid')).toThrow('Invalid numeric coordinate');
      expect(() => parseStrictDecimal('')).toThrow('Invalid numeric coordinate');
    });
  });

  describe('detectDelimiter', () => {
    it('detects semicolon delimiter', () => {
      const lines = ['# Comment', '101;324000;588000', '102;324100;588100'];
      expect(detectDelimiter(lines)).toBe(';');
    });

    it('detects comma delimiter', () => {
      const lines = ['', '# Header', '101,324000,588000'];
      expect(detectDelimiter(lines)).toBe(',');
    });

    it('detects tab delimiter', () => {
      const lines = ['101\t324000\t588000', '102\t324100\t588100'];
      expect(detectDelimiter(lines)).toBe('\t');
    });

    it('defaults to space delimiter when no punctuation delimiter is present', () => {
      const lines = ['101 324000 588000', '102 324100 588100'];
      expect(detectDelimiter(lines)).toBe(' ');
    });

    it('returns space if all lines are empty or comments', () => {
      expect(detectDelimiter(['', '   ', '# comment only'])).toBe(' ');
    });
  });

  describe('parseInputLines', () => {
    describe('Stereo 70 (planar) coordinates', () => {
      it('parses standard 2D coordinate lines (NE order)', () => {
        const text = '324000.500 588000.250\n325000 589000';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(2);
        expect(rows[0].pointId).toBeNull();
        expect(rows[0].coords).toEqual([324000.5, 588000.25]);
        expect(rows[0].error).toBeUndefined();

        expect(rows[1].pointId).toBeNull();
        expect(rows[1].coords).toEqual([325000, 589000]);
      });

      it('normalizes EN order to internal NE ([North, East]) format', () => {
        // Input: East, North -> output must be North, East
        const text = '588000.250 324000.500';
        const rows = parseInputLines(text, 'auto', 'EN', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].coords).toEqual([324000.5, 588000.25]);
        expect(rows[0].originalOrder).toBe('EN');
      });

      it('parses 3D coordinates (North, East, Elevation)', () => {
        const text = '324000 588000 120.75';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBeNull();
        expect(rows[0].coords).toEqual([324000, 588000, 120.75]);
      });

      it('parses 4 tokens: Point ID with 3D coordinates', () => {
        const text = 'P1 324000.5 588000.5 105.2';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBe('P1');
        expect(rows[0].coords).toEqual([324000.5, 588000.5, 105.2]);
      });

      it('distinguishes string Point ID + 2D coordinates from 3D coordinates', () => {
        const text = 'StationA 324000 588000';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBe('StationA');
        expect(rows[0].coords).toEqual([324000, 588000]);
      });

      it('distinguishes numeric Point ID when planar coordinates are > 50000', () => {
        // 101 is ID, 324000 and 588000 are Stereo 70 coordinates
        const text = '101 324000 588000';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBe('101');
        expect(rows[0].coords).toEqual([324000, 588000]);
      });

      it('handles various delimiters explicitly', () => {
        const csvText = 'P1,324000.0,588000.0\nP2,324100.0,588100.0';
        const csvRows = parseInputLines(csvText, ',', 'NE', false);
        expect(csvRows[0].pointId).toBe('P1');
        expect(csvRows[0].coords).toEqual([324000, 588000]);

        const semiText = 'P1;324000.0;588000.0';
        const semiRows = parseInputLines(semiText, ';', 'NE', false);
        expect(semiRows[0].coords).toEqual([324000, 588000]);

        const tabText = 'P1\t324000.0\t588000.0';
        const tabRows = parseInputLines(tabText, '\t', 'NE', false);
        expect(tabRows[0].coords).toEqual([324000, 588000]);
      });
    });

    describe('Geographic (ETRS89) coordinates', () => {
      it('parses decimal degrees (Lat, Lon)', () => {
        const text = '44.432312 26.106312';
        const rows = parseInputLines(text, 'auto', 'NE', true, 'degrees');

        expect(rows).toHaveLength(1);
        expect(rows[0].coords).toEqual([44.432312, 26.106312]);
      });

      it('parses DMS coordinate format with angularUnit = "dms"', () => {
        const text = '44°25\'56.28" 26°06\'22.72"';
        const rows = parseInputLines(text, ' ', 'NE', true, 'dms');

        expect(rows).toHaveLength(1);
        expect(rows[0].coords[0]).toBeCloseTo(44.4323, 3);
        expect(rows[0].coords[1]).toBeCloseTo(26.1063, 3);
      });

      it('parses radians format with angularUnit = "radians"', () => {
        // ~45° is 0.785398 rad
        const text = '0.785398 0.45564';
        const rows = parseInputLines(text, 'auto', 'NE', true, 'radians');

        expect(rows).toHaveLength(1);
        expect(rows[0].coords[0]).toBeCloseTo(45.0, 3);
      });

      it('recognizes Point ID with geographic coordinates', () => {
        const text = 'GPS_BASE 44.4323 26.1063';
        const rows = parseInputLines(text, 'auto', 'NE', true, 'degrees');

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBe('GPS_BASE');
        expect(rows[0].coords).toEqual([44.4323, 26.1063]);
      });

      it('parses 4-token geographic line (ID, Lat, Lon, Elevation)', () => {
        const text = 'STATION_9 44.432312 26.106312 110.5';
        const rows = parseInputLines(text, 'auto', 'NE', true, 'degrees');

        expect(rows).toHaveLength(1);
        expect(rows[0].pointId).toBe('STATION_9');
        expect(rows[0].coords).toEqual([44.432312, 26.106312, 110.5]);
      });

      it('normalizes EN geographic input order to [Lat, Lon]', () => {
        // Input: Lon (26.1063), Lat (44.4323) -> normalized: [44.4323, 26.1063]
        const text = '26.1063 44.4323';
        const rows = parseInputLines(text, 'auto', 'EN', true, 'degrees');

        expect(rows).toHaveLength(1);
        expect(rows[0].coords).toEqual([44.4323, 26.1063]);
      });

      it('handles mixed Windows CRLF and Unix LF line endings', () => {
        const text = "324000 588000\r\n325000 589000\n326000 590000";
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(3);
        expect(rows[0].coords).toEqual([324000, 588000]);
        expect(rows[1].coords).toEqual([325000, 589000]);
        expect(rows[2].coords).toEqual([326000, 590000]);
      });
    });

    describe('blank lines, comments, and errors', () => {
      it('flags empty lines without throwing', () => {
        const text = '324000 588000\n\n   \n325000 589000';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(4);
        expect(rows[0].coords).toEqual([324000, 588000]);
        expect(rows[1].isEmpty).toBe(true);
        expect(rows[2].isEmpty).toBe(true);
        expect(rows[3].coords).toEqual([325000, 589000]);
      });

      it('flags comment lines starting with #', () => {
        const text = '# This is a header\n324000 588000\n# End comment';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(3);
        expect(rows[0].isComment).toBe(true);
        expect(rows[1].coords).toEqual([324000, 588000]);
        expect(rows[2].isComment).toBe(true);
      });

      it('records line-level errors for lines with too few columns', () => {
        const text = '324000';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].error).toBe('Too few columns (minimum 2 coordinates required)');
        expect(rows[0].coords).toHaveLength(0);
      });

      it('records line-level errors for invalid coordinate tokens', () => {
        const text = 'P1 324000 not_a_number';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows).toHaveLength(1);
        expect(rows[0].error).toBeDefined();
        expect(rows[0].coords).toHaveLength(0);
      });

      it('preserves line indices correctly across mixed lines', () => {
        const text = '# Comment (0)\n\n324000 588000 (2)\nBadLine (3)';
        const rows = parseInputLines(text, 'auto', 'NE', false);

        expect(rows[0].lineIndex).toBe(0);
        expect(rows[1].lineIndex).toBe(1);
        expect(rows[2].lineIndex).toBe(2);
        expect(rows[3].lineIndex).toBe(3);
      });
    });
  });
});
