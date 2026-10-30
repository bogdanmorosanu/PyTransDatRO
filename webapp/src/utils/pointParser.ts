import { parseAngleToDeg, radToDeg } from './dmsParser';

export type CoordinateOrder = 'NE' | 'EN';

export interface ParsedPointRow {
  lineIndex: number;
  rawLine: string;
  pointId: string | null;
  coords: number[]; // Normalized to [North/Lat, East/Long, (Elevation)]
  originalOrder: CoordinateOrder;
  error?: string;
  errorType?: 'invalid_input' | 'invalid_angle';
  isComment?: boolean;
  isEmpty?: boolean;
}

export function isValidDecimal(s: string): boolean {
  const trimmed = s.trim();
  return /^[+-]?(?:\d+\.?\d*|\.\d+)$/.test(trimmed) && !isNaN(Number(trimmed));
}

/**
 * Checks if a token contains explicit DMS or sexagesimal symbols.
 * Distinguishes true angular markers (degree/minute/second symbols or digits followed by d/m/s)
 * from arbitrary alphabetical words/place names like 'Bucuresti' or 'Timisoara'.
 */
export function hasDmsSymbols(s: string): boolean {
  return /[°'"\u2018\u2019\u201C\u201D`]|(?:\d\s*[dms])/i.test(s);
}

/**
 * Checks whether a token can realistically represent an angle (either decimal, DMS with symbols, or signed compass heading).
 */
export function isPotentialAngle(s: string): boolean {
  const trimmed = s.trim();
  if (!trimmed) return false;
  if (isValidDecimal(trimmed)) return true;
  if (hasDmsSymbols(trimmed)) {
    try {
      const val = parseAngleToDeg(trimmed);
      return !isNaN(val);
    } catch {
      return false;
    }
  }
  return false;
}

export interface ExtractedPointTokens {
  pointId: string | null;
  coordTokens: string[];
}

/**
 * Centralized extractor that analyzes delimited tokens to determine whether
 * the first column represents a Point ID or the start of coordinate values.
 * Shared across Text (batch) parsing and File upload mapping.
 */
export function extractPointIdAndCoords(
  tokens: string[],
  isSourceGeo: boolean,
  order: CoordinateOrder = 'NE'
): ExtractedPointTokens {
  if (tokens.length === 0) {
    return { pointId: null, coordTokens: [] };
  }
  if (tokens.length === 1) {
    return { pointId: null, coordTokens: tokens };
  }
  if (tokens.length === 2) {
    return { pointId: null, coordTokens: tokens };
  }
  if (tokens.length >= 4) {
    return { pointId: tokens[0], coordTokens: tokens.slice(1, 4) };
  }

  // Exactly 3 tokens: [t0, t1, t2]
  const [t0, t1, t2] = tokens;

  if (isSourceGeo) {
    const hasAngleSymbol0 = hasDmsSymbols(t0);
    const pot0 = isPotentialAngle(t0);
    const pot1 = isPotentialAngle(t1);
    const pot2 = isPotentialAngle(t2);

    let t0Deg = NaN;
    try {
      t0Deg = parseAngleToDeg(t0);
    } catch {}

    // Case 1: t0 is non-numeric string label (e.g. "P1", "Station_A")
    // or numeric ID > 90 (which cannot be Latitude in NE order), while t1 and t2 can be angles:
    if (!hasAngleSymbol0 && pot1 && pot2 && (!pot0 || (order === 'NE' && !isNaN(t0Deg) && Math.abs(t0Deg) > 90))) {
      return { pointId: t0, coordTokens: [t1, t2] };
    }

    // Case 2: All 3 tokens are numbers/angles -> treat as 3D (Lat, Lon, Height)
    return { pointId: null, coordTokens: tokens };
  } else {
    // Stereo 70 (metric coordinates)
    const t0Num = isValidDecimal(t0) ? Number(t0) : null;
    const t1Num = isValidDecimal(t1) ? Number(t1) : null;
    const t2Num = isValidDecimal(t2) ? Number(t2) : null;

    if (t1Num !== null && t1Num > 50000) {
      if (t2Num !== null && t2Num > 50000) {
        // Both t1 and t2 are planimetric coordinates -> t0 is Point ID
        return { pointId: t0, coordTokens: [t1, t2] };
      } else {
        // t2 is elevation or other -> 3D coordinates without ID
        return { pointId: null, coordTokens: tokens };
      }
    } else if (t0Num !== null && t0Num > 50000) {
      // t0 is already a planimetric coordinate -> 3D coordinates without ID
      return { pointId: null, coordTokens: tokens };
    } else if (t0Num !== null && t1Num !== null && t2Num !== null) {
      // All 3 are numbers -> 3D coordinates
      return { pointId: null, coordTokens: tokens };
    } else if (!isValidDecimal(t0) && isValidDecimal(t1) && isValidDecimal(t2)) {
      // t0 is string label (e.g. "Bucuresti", "P1"), t1 and t2 are coordinates
      return { pointId: t0, coordTokens: [t1, t2] };
    } else {
      return { pointId: null, coordTokens: tokens };
    }
  }
}

export function parseStrictDecimal(s: string): number {
  const trimmed = s.trim();
  if (!isValidDecimal(trimmed)) {
    throw new Error('Invalid numeric coordinate');
  }
  return Number(trimmed);
}

/**
 * Autodetects the primary column delimiter from non-empty lines.
 */
export function detectDelimiter(lines: string[]): string {
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;
    if (trimmed.includes(';')) return ';';
    if (trimmed.includes(',')) return ',';
    if (trimmed.includes('\t')) return '\t';
  }
  return ' '; // Default to whitespace if no other delimiter found
}

/**
 * Parses multiline input text with optional Point ID and delimiter detection.
 */
export function parseInputLines(
  rawText: string,
  delimiterConfig: string = 'auto',
  order: CoordinateOrder = 'NE',
  isSourceGeo: boolean = false,
  angularUnit: 'degrees' | 'dms' | 'radians' = 'degrees'
): ParsedPointRow[] {
  const lines = rawText.split(/\r\n|\r|\n/);
  const rows: ParsedPointRow[] = [];

  const delimiter = delimiterConfig === 'auto' ? detectDelimiter(lines) : delimiterConfig;
  const splitRegex = delimiter === ' ' ? /\s+/ : new RegExp(delimiter === ';' ? ';+' : delimiter === ',' ? ',+' : '\\t+');

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const trimmed = line.trim();
    if (!trimmed) {
      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId: null,
        coords: [],
        originalOrder: order,
        isEmpty: true,
      });
      continue;
    }
    if (trimmed.startsWith('#')) {
      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId: null,
        coords: [],
        originalOrder: order,
        isComment: true,
      });
      continue;
    }

    const tokens = trimmed.split(splitRegex).map((t) => t.trim()).filter((t) => t.length > 0);

    if (tokens.length < 2) {
      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId: null,
        coords: [],
        originalOrder: order,
        error: 'Too few columns (minimum 2 coordinates required)',
        errorType: 'invalid_input',
      });
      continue;
    }

    const { pointId, coordTokens } = extractPointIdAndCoords(tokens, isSourceGeo, order);

    if (coordTokens.length < 2 || coordTokens.length > 3) {
      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId,
        coords: [],
        originalOrder: order,
        error: `Invalid number of coordinate columns: ${coordTokens.length}`,
        errorType: 'invalid_input',
      });
      continue;
    }

    try {
      let c1: number;
      let c2: number;
      let c3: number | undefined;

      if (isSourceGeo) {
        const parseGeoToken = (tok: string): number => {
          if (angularUnit === 'radians') {
            if (tok.includes('°') || tok.includes("'") || tok.includes('"')) {
              return parseAngleToDeg(tok);
            }
            const clean = tok.replace(/rad/gi, '').trim();
            const num = Number(clean);
            if (!isNaN(num)) return radToDeg(num);
          }
          return parseAngleToDeg(tok);
        };
        c1 = parseGeoToken(coordTokens[0]);
        c2 = parseGeoToken(coordTokens[1]);

        if (isNaN(c1) || isNaN(c2)) {
          rows.push({
            lineIndex: i,
            rawLine: line,
            pointId,
            coords: [],
            originalOrder: order,
            error: 'Invalid coordinate syntax',
            errorType: 'invalid_input',
          });
          continue;
        }

        const latVal = order === 'EN' ? c2 : c1;
        const lonVal = order === 'EN' ? c1 : c2;
        if (Math.abs(latVal) > 90 || Math.abs(lonVal) > 180) {
          rows.push({
            lineIndex: i,
            rawLine: line,
            pointId,
            coords: [],
            originalOrder: order,
            error: 'Invalid angle',
            errorType: 'invalid_angle',
          });
          continue;
        }
      } else {
        c1 = parseStrictDecimal(coordTokens[0]);
        c2 = parseStrictDecimal(coordTokens[1]);
      }

      if (coordTokens.length === 3) {
        c3 = parseStrictDecimal(coordTokens[2]);
      }

      // Normalize to [North/Lat, East/Long, (Elev)]
      let normalizedCoords: number[];
      if (order === 'EN') {
        // Input is East, North -> normalize to North, East
        normalizedCoords = c3 !== undefined ? [c2, c1, c3] : [c2, c1];
      } else {
        // Input is North, East
        normalizedCoords = c3 !== undefined ? [c1, c2, c3] : [c1, c2];
      }

      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId,
        coords: normalizedCoords,
        originalOrder: order,
      });
    } catch (err: any) {
      rows.push({
        lineIndex: i,
        rawLine: line,
        pointId,
        coords: [],
        originalOrder: order,
        error: err.message || 'Invalid coordinate format',
        errorType: 'invalid_input',
      });
    }
  }

  return rows;
}
