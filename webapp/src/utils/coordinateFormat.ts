/**
 * Canonical Coordinate and Geodetic Formatting Utilities
 * Eliminates duplicate formatting routines across Single Point, Batch, File, and Map views.
 */

import { degToRad, formatDegToDMS, type AngularUnit } from './dmsParser';

export type { AngularUnit };

/**
 * Formats a geographic angular coordinate (latitude or longitude) based on the active unit.
 * - 'dms': Formats as DD°MM'SS.sssss" (default 5 decimal seconds)
 * - 'radians': Formats in decimal radians (default 10 decimals)
 * - 'degrees': Formats in decimal degrees (default 8 decimals)
 */
export function formatAngularCoordinate(
  deg: number,
  unit: AngularUnit,
  precision?: number
): string {
  if (isNaN(deg)) return '';
  if (unit === 'dms') {
    return formatDegToDMS(deg, precision ?? 5);
  }
  if (unit === 'radians') {
    return degToRad(deg).toFixed(precision ?? 10);
  }
  return deg.toFixed(precision ?? 8);
}

/**
 * Formats a planar projected coordinate (Stereo 70 Easting or Northing) in meters (default 3 decimals).
 */
export function formatPlanarCoordinate(meters: number, precision: number = 3): string {
  if (isNaN(meters)) return '';
  return meters.toFixed(precision);
}

/**
 * Formats a pair of 2D coordinates respecting angular unit and coordinate ordering ('NE' or 'EN').
 * If isGeographic is true, c1 is Latitude and c2 is Longitude.
 * If isGeographic is false, c1 is Northing and c2 is Easting.
 */
export function formatCoordinatePair(
  c1: number,
  c2: number,
  unit: AngularUnit,
  order: 'NE' | 'EN' = 'NE',
  isGeographic: boolean = true
): {
  col1: string;
  col2: string;
  label1: string;
  label2: string;
  latStr?: string;
  lonStr?: string;
  nStr?: string;
  eStr?: string;
} {
  if (isGeographic) {
    const latStr = formatAngularCoordinate(c1, unit);
    const lonStr = formatAngularCoordinate(c2, unit);
    return {
      col1: order === 'EN' ? lonStr : latStr,
      col2: order === 'EN' ? latStr : lonStr,
      label1: order === 'EN' ? 'Lon' : 'Lat',
      label2: order === 'EN' ? 'Lat' : 'Lon',
      latStr,
      lonStr,
    };
  } else {
    const nStr = formatPlanarCoordinate(c1);
    const eStr = formatPlanarCoordinate(c2);
    return {
      col1: order === 'EN' ? eStr : nStr,
      col2: order === 'EN' ? nStr : eStr,
      label1: order === 'EN' ? 'E' : 'N',
      label2: order === 'EN' ? 'N' : 'E',
      nStr,
      eStr,
    };
  }
}

/**
 * Calculates quasigeoid undulation / height anomaly (zeta = h - H).
 * - h: Ellipsoidal height (ETRS89 / GRS80) in meters
 * - H: Normal / orthometric height (Black Sea 1975) in meters
 * Returns null if either height is undefined, null, or NaN.
 */
export function computeQuasigeoidZeta(
  ellipsoidalHeight?: number | null,
  normalHeight?: number | null
): number | null {
  if (
    ellipsoidalHeight === undefined ||
    ellipsoidalHeight === null ||
    isNaN(ellipsoidalHeight) ||
    normalHeight === undefined ||
    normalHeight === null ||
    isNaN(normalHeight)
  ) {
    return null;
  }
  return ellipsoidalHeight - normalHeight;
}

/**
 * Formats quasigeoid undulation (zeta) with explicit positive sign and 3 decimals (e.g. "+0.345" or "-0.120").
 * Returns an empty string if zeta is null, undefined, or NaN.
 */
export function formatQuasigeoidZeta(zeta: number | null | undefined): string {
  if (zeta === null || zeta === undefined || isNaN(zeta)) return '';
  const sign = zeta >= 0 ? '+' : '';
  return `${sign}${zeta.toFixed(3)}`;
}
