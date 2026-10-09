/**
 * Geodetic Angular Conversion Utilities
 * Handles decimal degrees, radians, and degrees-minutes-seconds (DMS) string parsing.
 */

export const DEG_TO_RAD = Math.PI / 180.0;
export const RAD_TO_DEG = 180.0 / Math.PI;

export function degToRad(deg: number): number {
  return deg * DEG_TO_RAD;
}

export function radToDeg(rad: number): number {
  return rad * RAD_TO_DEG;
}

/**
 * Parses arbitrary angular representations into decimal degrees.
 * Supports:
 * - Decimal numbers: "45.999718"
 * - DMS with symbols: "45° 59' 58.99\" N", "45d 59m 58.99s", "45 59 58.99"
 * - Cardinal direction indicators (N, S, E, W)
 */
export function parseAngleToDeg(input: string): number {
  const raw = input.trim();
  if (!raw) throw new Error('Empty coordinate input');

  // Check cardinal directions
  let sign = 1;
  let clean = raw.toUpperCase();
  if (clean.includes('S') || clean.includes('W') || clean.startsWith('-')) {
    sign = -1;
  }
  clean = clean.replace(/[NSWE\-+]/g, '').trim();

  // If simple double
  const simpleDouble = Number(clean);
  if (!isNaN(simpleDouble) && !clean.includes('°') && !clean.includes("'") && !clean.includes('"') && !clean.includes('D') && !clean.includes(' ')) {
    return sign * simpleDouble;
  }

  // Split DMS components by delimiters
  const tokens = clean
    .split(/[°'"`dms\s]+/)
    .map((s) => s.trim())
    .filter((s) => s.length > 0);

  if (tokens.length === 0) {
    throw new Error(`Invalid DMS format: "${input}"`);
  }

  const d = Math.abs(parseFloat(tokens[0]) || 0);
  const m = tokens.length > 1 ? parseFloat(tokens[1]) || 0 : 0;
  const s = tokens.length > 2 ? parseFloat(tokens[2]) || 0 : 0;

  if (m < 0 || m >= 60) {
    throw new Error(`Minutes must be between 0 and 59: ${m}`);
  }
  if (s < 0 || s >= 60) {
    throw new Error(`Seconds must be between 0 and 59.999: ${s}`);
  }

  const totalDeg = d + m / 60.0 + s / 3600.0;
  return sign * totalDeg;
}

export function parseAngleToRad(input: string): number {
  return degToRad(parseAngleToDeg(input));
}

/**
 * Formats decimal degrees into standard geodetic DMS string: DD° MM' SS.sssss"
 */
export function formatDegToDMS(deg: number, precision: number = 5): string {
  const sign = deg < 0 ? '-' : '';
  let absDeg = Math.abs(deg);

  let d = Math.floor(absDeg);
  let rem = (absDeg - d) * 60.0;
  let m = Math.floor(rem);
  let s = (rem - m) * 60.0;

  // Round seconds according to requested precision
  const factor = Math.pow(10, precision);
  s = Math.round(s * factor) / factor;

  // Handle rounding rollover
  if (s >= 60.0) {
    s = 0;
    m += 1;
    if (m >= 60) {
      m = 0;
      d += 1;
    }
  }

  const sStr = s.toFixed(precision).padStart(precision + 3, '0');
  const mStr = String(m).padStart(2, '0');
  const dStr = String(d);

  return `${sign}${dStr}° ${mStr}' ${sStr}"`;
}
