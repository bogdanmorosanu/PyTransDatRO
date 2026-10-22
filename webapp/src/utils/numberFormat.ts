import { getLanguage } from '../i18n';

/**
 * Formats a numeric value using locale-aware rules for the telemetry/usage dashboard.
 * Romanian ('ro') formats decimals with comma (',') and thousands with dot ('.').
 * English ('en') formats decimals with dot ('.') and thousands with comma (',').
 *
 * NOTE: This function is strictly for Telemetry/Usage UI displays.
 * Geodetic coordinates and CAD/GIS export files must NEVER use this formatter.
 */
export function formatTelemetryNumber(num: number | string, options?: Intl.NumberFormatOptions): string {
  const n = typeof num === 'string' ? parseFloat(num) : num;
  if (isNaN(n)) return '0';
  const locale = getLanguage() === 'ro' ? 'ro-RO' : 'en-US';
  return n.toLocaleString(locale, options);
}

/**
 * Formats large values with localized compact notation (e.g. '1,2M' in RO, '1.2M' in EN).
 */
export function formatTelemetryCompact(num: number | string, fractionDigits: number = 1): string {
  const n = typeof num === 'string' ? parseFloat(num) : num;
  if (isNaN(n)) return '0';
  const isRo = getLanguage() === 'ro';
  const locale = isRo ? 'ro-RO' : 'en-US';

  if (n >= 1_000_000) {
    let val = (n / 1_000_000).toLocaleString(locale, {
      minimumFractionDigits: fractionDigits,
      maximumFractionDigits: fractionDigits,
    });
    if (isRo && val.includes('.')) {
      const parts = (n / 1_000_000).toFixed(fractionDigits).split('.');
      val = `${parts[0]},${parts[1]}`;
    }
    return `${val}M`;
  }
  if (n >= 1_000) {
    let val = (n / 1_000).toLocaleString(locale, {
      minimumFractionDigits: 0,
      maximumFractionDigits: fractionDigits,
    });
    if (isRo && val.includes('.')) {
      const fixed = (n / 1_000).toFixed(fractionDigits);
      const parts = fixed.split('.');
      const frac = parts[1]?.replace(/0+$/, '');
      val = frac && frac.length > 0 ? `${parts[0]},${frac}` : parts[0];
    }
    return `${val}k`;
  }
  return n.toLocaleString(locale);
}

/**
 * Formats a decimal number with fixed fraction digits according to active locale.
 * Romanian ('ro') formats decimals with comma (',') (e.g. '12,3').
 * English ('en') formats decimals with dot ('.') (e.g. '12.3').
 */
export function formatTelemetryDecimal(num: number | string, fractionDigits: number = 1): string {
  const n = typeof num === 'string' ? parseFloat(num) : num;
  if (isNaN(n)) return '0';
  const isRo = getLanguage() === 'ro';
  const locale = isRo ? 'ro-RO' : 'en-US';
  let formatted = n.toLocaleString(locale, {
    minimumFractionDigits: fractionDigits,
    maximumFractionDigits: fractionDigits,
  });
  if (isRo && formatted.includes('.')) {
    const parts = n.toFixed(fractionDigits).split('.');
    const intFormatted = Number(parts[0]).toLocaleString('ro-RO');
    formatted = `${intFormatted},${parts[1]}`;
  }
  return formatted;
}

/**
 * Formats a percentage value according to active locale.
 * Romanian ('ro') formats e.g. '24,5%'.
 * English ('en') formats e.g. '24.5%'.
 */
export function formatTelemetryPercent(num: number | string, fractionDigits: number = 1): string {
  const n = typeof num === 'string' ? parseFloat(num) : num;
  if (isNaN(n)) return '0%';
  const isRo = getLanguage() === 'ro';
  const locale = isRo ? 'ro-RO' : 'en-US';
  let formatted = n.toLocaleString(locale, {
    minimumFractionDigits: 0,
    maximumFractionDigits: fractionDigits,
  });
  if (isRo && formatted.includes('.')) {
    const fixed = n.toFixed(fractionDigits);
    const parts = fixed.split('.');
    const frac = parts[1]?.replace(/0+$/, '');
    const intFormatted = Number(parts[0]).toLocaleString('ro-RO');
    formatted = frac && frac.length > 0 ? `${intFormatted},${frac}` : intFormatted;
  }
  return `${formatted}%`;
}
