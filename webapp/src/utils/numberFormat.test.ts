import { describe, it, expect, beforeEach } from 'vitest';
import { setLanguage } from '../i18n';
import { formatTelemetryNumber, formatTelemetryCompact, formatTelemetryDecimal, formatTelemetryPercent, formatDisplayDate } from './numberFormat';

describe('numberFormat utility', () => {
  describe('English locale (en)', () => {
    beforeEach(() => {
      setLanguage('en');
    });

    it('formats thousands with comma and decimals with dot', () => {
      expect(formatTelemetryNumber(1234567)).toBe('1,234,567');
      expect(formatTelemetryNumber(1234.56, { minimumFractionDigits: 2 })).toBe('1,234.56');
    });

    it('formats compact millions with dot', () => {
      expect(formatTelemetryCompact(1200000)).toBe('1.2M');
      expect(formatTelemetryCompact(5500000)).toBe('5.5M');
    });

    it('formats compact thousands', () => {
      expect(formatTelemetryCompact(2500)).toBe('2.5k');
      expect(formatTelemetryCompact(2000)).toBe('2k');
    });

    it('handles small numbers and 0', () => {
      expect(formatTelemetryCompact(42)).toBe('42');
      expect(formatTelemetryCompact(0)).toBe('0');
    });

    it('formats decimal numbers with dot', () => {
      expect(formatTelemetryDecimal(12.34, 1)).toBe('12.3');
      expect(formatTelemetryDecimal(0, 1)).toBe('0.0');
    });

    it('formats percentages with dot', () => {
      expect(formatTelemetryPercent(24.5, 1)).toBe('24.5%');
      expect(formatTelemetryPercent(100, 1)).toBe('100%');
    });
  });

  describe('Romanian locale (ro)', () => {
    beforeEach(() => {
      setLanguage('ro');
    });

    it('formats thousands with dot and decimals with comma', () => {
      const formatted = formatTelemetryNumber(1234.56, { minimumFractionDigits: 2 });
      // In Romanian locale, decimal separator is comma (,)
      expect(formatted).toContain(',');
    });

    it('formats compact millions with comma', () => {
      const compactM = formatTelemetryCompact(1200000);
      expect(compactM).toBe('1,2M');
    });

    it('formats compact thousands with comma if fractional', () => {
      const compactK = formatTelemetryCompact(2500);
      expect(compactK).toBe('2,5k');
    });

    it('formats decimal numbers with comma', () => {
      expect(formatTelemetryDecimal(12.34, 1)).toBe('12,3');
      expect(formatTelemetryDecimal(0, 1)).toBe('0,0');
    });

    it('formats percentages with comma', () => {
      expect(formatTelemetryPercent(24.5, 1)).toBe('24,5%');
      expect(formatTelemetryPercent(100, 1)).toBe('100%');
    });

    it('handles numeric string inputs correctly in Romanian', () => {
      expect(formatTelemetryDecimal('12.34', 1)).toBe('12,3');
      expect(formatTelemetryPercent('24.5', 1)).toBe('24,5%');
    });
  });

  describe('formatDisplayDate utility', () => {
    it('formats standard ISO YYYY-MM-DD dates to UK standard DD/MM/YY', () => {
      expect(formatDisplayDate('2026-10-24')).toBe('24/10/26');
      expect(formatDisplayDate('2023-11-01')).toBe('01/11/23');
      expect(formatDisplayDate('2024-01-05')).toBe('05/01/24');
      expect(formatDisplayDate('2024-12-31')).toBe('31/12/24');
    });

    it('formats month-only ISO YYYY-MM dates to 01/MM/YY', () => {
      expect(formatDisplayDate('2023-11')).toBe('01/11/23');
      expect(formatDisplayDate('2025-06')).toBe('01/06/25');
    });

    it('handles whitespace, null, undefined, or empty strings gracefully', () => {
      expect(formatDisplayDate('  2026-05-15  ')).toBe('15/05/26');
      expect(formatDisplayDate(null)).toBe('');
      expect(formatDisplayDate(undefined)).toBe('');
      expect(formatDisplayDate('')).toBe('');
    });
  });
});
