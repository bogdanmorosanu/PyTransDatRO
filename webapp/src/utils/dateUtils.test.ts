import { describe, it, expect } from 'vitest';
import { calculatePresetDateRange, formatDisplayDate } from './dateUtils';

describe('dateUtils', () => {
  const mockBaseDate = new Date('2026-10-15T12:00:00Z');

  describe('calculatePresetDateRange', () => {
    it('calculates 7d preset (7 days before base date)', () => {
      const range = calculatePresetDateRange('7d', '2011-01-09', undefined, mockBaseDate);
      expect(range.startDate).toBe('2026-10-08');
      expect(range.endDate).toBe('2026-10-15');
    });

    it('calculates 30d preset (30 days before base date)', () => {
      const range = calculatePresetDateRange('30d', '2011-01-09', undefined, mockBaseDate);
      expect(range.startDate).toBe('2026-09-15');
      expect(range.endDate).toBe('2026-10-15');
    });

    it('calculates ytd preset (January 1st of current year)', () => {
      const range = calculatePresetDateRange('ytd', '2011-01-09', undefined, mockBaseDate);
      expect(range.startDate).toBe('2026-01-01');
      expect(range.endDate).toBe('2026-10-15');
    });

    it('calculates 12m preset (1 year before base date)', () => {
      const range = calculatePresetDateRange('12m', '2011-01-09', undefined, mockBaseDate);
      expect(range.startDate).toBe('2025-10-15');
      expect(range.endDate).toBe('2026-10-15');
    });

    it('calculates all preset (defaults to minDate)', () => {
      const range = calculatePresetDateRange('all', '2023-11-01', undefined, mockBaseDate);
      expect(range.startDate).toBe('2023-11-01');
      expect(range.endDate).toBe('2026-10-15');
    });

    it('honors custom maxDate if provided', () => {
      const range = calculatePresetDateRange('7d', '2011-01-09', '2026-09-30', mockBaseDate);
      expect(range.endDate).toBe('2026-09-30');
    });
  });

  describe('formatDisplayDate', () => {
    it('formats YYYY-MM-DD into DD/MM/YY', () => {
      expect(formatDisplayDate('2026-10-15')).toBe('15/10/26');
      expect(formatDisplayDate('2011-01-09')).toBe('09/01/11');
    });

    it('formats YYYY-MM into 01/MM/YY', () => {
      expect(formatDisplayDate('2026-10')).toBe('01/10/26');
    });

    it('returns empty string for null, undefined, or empty', () => {
      expect(formatDisplayDate(null)).toBe('');
      expect(formatDisplayDate(undefined)).toBe('');
      expect(formatDisplayDate('')).toBe('');
    });
  });
});
