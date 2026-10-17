/**
 * Reference points from the official ANCPI benchmark files in tests/data/.
 *
 * Source rows (id = first CSV column):
 *  - st70_to_etrs89_input.csv     rows 1-3: Stereo 70 N, E (m) and quasigeoid height (m)
 *  - st70_to_etrs89_expected.csv  rows 1-3: ETRS89 lat/lon (DMS, Windows-1252) and ellipsoidal h
 *  - etrs89_to_st70_input.csv / etrs89_to_st70_expected.csv rows 1-3: the reverse direction; its
 *    expected Stereo 70 values are rounded to 0.01 m, so the full-precision values from the
 *    forward input file are used here for both directions.
 *
 * DMS to decimal degrees: deg + min/60 + sec/3600, 12 decimals.
 * Stereo 70 column order in the CSVs is `id, N, E, H` (northing first).
 */
export interface ReferencePoint {
  /** Row id in the benchmark CSVs. */
  id: number;
  /** Stereo 70 northing (m). */
  north: number;
  /** Stereo 70 easting (m). */
  east: number;
  /** Quasigeoid height (m), the third Stereo 70 coordinate. */
  height: number;
  /** ETRS89 latitude, decimal degrees. */
  lat: number;
  /** ETRS89 longitude, decimal degrees. */
  lon: number;
  /** ETRS89 ellipsoidal height (m). */
  ellipsoidalHeight: number;
}

export const REFERENCE_POINTS: readonly ReferencePoint[] = [
  {
    id: 1,
    north: 230134.564,
    east: 434283.0401,
    height: 100,
    lat: 43.568075511111, // 43°34'05.07184"
    lon: 24.185193983333, // 24°11'06.69834"
    ellipsoidalHeight: 141.091,
  },
  {
    id: 2,
    north: 230134.564,
    east: 439783.0401,
    height: 100,
    lat: 43.568551188889, // 43°34'06.78428"
    lon: 24.253258238889, // 24°15'11.72966"
    ellipsoidalHeight: 141.02,
  },
  {
    id: 3,
    north: 230134.564,
    east: 445283.0401,
    height: 100,
    lat: 43.568985311111, // 43°34'08.34712"
    lon: 24.321323708333, // 24°19'16.76535"
    ellipsoidalHeight: 140.971,
  },
];

/** Decimals shown by the UI: ETRS89 degrees (8), metres (3). See BatchConverter / PointConverter. */
export const DISPLAY_DECIMALS = { degrees: 8, metres: 3 } as const;

/** True if `actual` equals `expected` at `decimals` places, allowing +-1 in the last digit. */
export function matchesAtDisplayPrecision(
  actual: number | string,
  expected: number,
  decimals: number,
): boolean {
  const value = typeof actual === 'string' ? Number(actual.trim()) : actual;
  if (!Number.isFinite(value)) return false;
  const unit = 10 ** -decimals;
  // Small epsilon absorbs binary floating-point error at the +-1 unit boundary.
  return Math.abs(value - expected) <= unit * (1 + 1e-9);
}
