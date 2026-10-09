import { parseAngleToDeg } from './dmsParser';

export type CoordinateOrder = 'NE' | 'EN';

export interface ParsedPointRow {
  rawLine: string;
  pointId: string | null;
  coords: number[]; // Normalized to [North/Lat, East/Long, (Elevation)]
  originalOrder: CoordinateOrder;
  error?: string;
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
  isSourceGeo: boolean = false
): ParsedPointRow[] {
  const lines = rawText.split(/\r\n|\r|\n/);
  const rows: ParsedPointRow[] = [];

  const delimiter = delimiterConfig === 'auto' ? detectDelimiter(lines) : delimiterConfig;
  const splitRegex = delimiter === ' ' ? /\s+/ : new RegExp(delimiter === ';' ? ';+' : delimiter === ',' ? ',+' : '\\t+');

  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (trimmed.startsWith('#')) continue; // Skip comments

    const tokens = trimmed.split(splitRegex).map((t) => t.trim()).filter((t) => t.length > 0);

    if (tokens.length < 2) {
      rows.push({
        rawLine: line,
        pointId: null,
        coords: [],
        originalOrder: order,
        error: 'Too few columns (minimum 2 coordinates required)',
      });
      continue;
    }

    let pointId: string | null = null;
    let coordTokens: string[] = [];

    // Check if first column is Point ID:
    // If we have 3 or 4 tokens, and token 0 is not easily parsed as a coordinate
    // or if the line has 4 tokens (PointId, X, Y, Z)
    if (tokens.length === 4) {
      pointId = tokens[0];
      coordTokens = tokens.slice(1);
    } else if (tokens.length === 3) {
      // Could be (X, Y, Z) or (PointId, X, Y)
      const firstIsNumeric = !isNaN(Number(tokens[0])) && !isNaN(parseFloat(tokens[0]));
      const secondIsNumeric = !isNaN(Number(tokens[1])) && !isNaN(parseFloat(tokens[1]));
      const thirdIsNumeric = !isNaN(Number(tokens[2])) && !isNaN(parseFloat(tokens[2]));

      if (firstIsNumeric && secondIsNumeric && thirdIsNumeric) {
        // All numeric -> 3D coordinates (X, Y, Z)
        coordTokens = tokens;
      } else {
        // First is likely a label -> (PointId, X, Y)
        pointId = tokens[0];
        coordTokens = tokens.slice(1);
      }
    } else {
      // 2 tokens -> 2D coordinates (X, Y)
      coordTokens = tokens;
    }

    if (coordTokens.length < 2 || coordTokens.length > 3) {
      rows.push({
        rawLine: line,
        pointId,
        coords: [],
        originalOrder: order,
        error: `Invalid number of coordinate columns: ${coordTokens.length}`,
      });
      continue;
    }

    try {
      let c1: number;
      let c2: number;
      let c3: number | undefined;

      if (isSourceGeo) {
        c1 = parseAngleToDeg(coordTokens[0]);
        c2 = parseAngleToDeg(coordTokens[1]);
      } else {
        c1 = parseFloat(coordTokens[0]);
        c2 = parseFloat(coordTokens[1]);
        if (isNaN(c1) || isNaN(c2)) throw new Error('Non-numeric coordinate value');
      }

      if (coordTokens.length === 3) {
        c3 = parseFloat(coordTokens[2]);
        if (isNaN(c3)) throw new Error('Non-numeric elevation value');
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
        rawLine: line,
        pointId,
        coords: normalizedCoords,
        originalOrder: order,
      });
    } catch (err: any) {
      rows.push({
        rawLine: line,
        pointId,
        coords: [],
        originalOrder: order,
        error: err.message || 'Invalid coordinate format',
      });
    }
  }

  return rows;
}
