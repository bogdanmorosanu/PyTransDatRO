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

  // Check negative sign and cardinal directions
  let sign = 1;
  let clean = raw.trim();

  if (/^[-]\s*/.test(clean) || /^[SW]\s*\d/i.test(clean) || /^[SW]\b/i.test(clean)) {
    // Leading negative or S/W (e.g., "-45.5", "S 45.5", "S45°...")
    sign = -1;
    clean = clean.replace(/^[-+]\s*|^[NSWE]\s*/i, '');
  } else if (/\s+[SW]$/i.test(clean) || /["'`]\s*[SW]$/i.test(clean)) {
    // Trailing cardinal direction after whitespace or quote: e.g. "45° 30' S", "45.5 S", "45° 30' 00\" S"
    sign = -1;
    clean = clean.replace(/\s*[NSWE]$/i, '');
  } else if (/[SW]$/i.test(clean) && !/[d°m']/i.test(clean)) {
    // Trailing S/W directly attached to a decimal without dms symbols: e.g. "45.5S"
    sign = -1;
    clean = clean.replace(/[NSWE]$/i, '');
  } else {
    // Strip trailing N or E if present
    clean = clean.replace(/\s*[NE]$/i, '');
    if (/[NE]$/i.test(clean) && !/[d°m']/i.test(clean)) {
      clean = clean.replace(/[NE]$/i, '');
    }
    // Strip leading + or N/E
    clean = clean.replace(/^[+]\s*|^[NE]\s*/i, '');
  }

  clean = clean.trim();
  if (!clean) {
    throw new Error(`Invalid coordinate input: "${input}"`);
  }

  // If simple double without any symbols or separators
  const simpleDouble = Number(clean);
  if (!isNaN(simpleDouble) && !/[°'"`dms\s]/i.test(clean)) {
    return sign * simpleDouble;
  }

  // Split DMS components by delimiters (case-insensitive for d, m, s)
  const tokens = clean
    .split(/[°'"`dms\s]+/i)
    .map((s) => s.trim())
    .filter((s) => s.length > 0);

  if (tokens.length === 0 || tokens.length > 3) {
    throw new Error(`Invalid DMS format: "${input}"`);
  }

  for (const tok of tokens) {
    if (!/^[+-]?(?:\d+\.?\d*|\.\d+)$/.test(tok) || isNaN(Number(tok))) {
      throw new Error('Invalid numeric coordinate');
    }
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
 * Formats decimal degrees into standard geodetic DMS string: DD°MM'SS.sssss" (no spaces)
 */
export function formatDegToDMS(deg: number, precision: number = 5): string {
  if (isNaN(deg)) throw new Error('Cannot format NaN to DMS');
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

  const sStr = precision > 0
    ? s.toFixed(precision).padStart(precision + 3, '0')
    : String(Math.round(s)).padStart(2, '0');
  const mStr = String(m).padStart(2, '0');
  const dStr = String(d);

  return `${sign}${dStr}°${mStr}'${sStr}"`;
}
