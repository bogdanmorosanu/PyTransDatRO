/**
 * Canonical Date & Time Range Utilities for Telemetry & Dashboard Analytics
 * Eliminates duplicate preset calculations across computation, legacy, and spatial views.
 */

import { formatDisplayDate } from './numberFormat';

export { formatDisplayDate };

export type TelemetryPreset = 'all' | '12m' | 'ytd' | '30d' | '7d';

export interface PresetDateRange {
  startDate: string;
  endDate: string;
}

/**
 * Calculates start and end ISO dates (YYYY-MM-DD) for a given relative telemetry preset.
 *
 * @param preset - The selected preset ('7d', '30d', 'ytd', '12m', or 'all')
 * @param minDate - Earliest available date in database (used for 'all')
 * @param maxDate - Custom upper bound (defaults to today's date)
 * @param baseDate - Reference date object (defaults to current wall clock)
 */
export function calculatePresetDateRange(
  preset: TelemetryPreset,
  minDate: string = '2011-01-09',
  maxDate?: string,
  baseDate: Date = new Date()
): PresetDateRange {
  const endDateStr = maxDate || baseDate.toISOString().slice(0, 10);

  if (preset === '7d') {
    const d = new Date(baseDate);
    d.setDate(d.getDate() - 7);
    return { startDate: d.toISOString().slice(0, 10), endDate: endDateStr };
  }

  if (preset === '30d') {
    const d = new Date(baseDate);
    d.setDate(d.getDate() - 30);
    return { startDate: d.toISOString().slice(0, 10), endDate: endDateStr };
  }

  if (preset === 'ytd') {
    return { startDate: `${baseDate.getFullYear()}-01-01`, endDate: endDateStr };
  }

  if (preset === '12m') {
    const d = new Date(baseDate);
    d.setFullYear(d.getFullYear() - 1);
    return { startDate: d.toISOString().slice(0, 10), endDate: endDateStr };
  }

  // 'all' or fallback
  return { startDate: minDate, endDate: endDateStr };
}
