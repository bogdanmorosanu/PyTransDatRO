/**
 * Core SVG and Canvas Charting Geometry Math
 * Eliminates duplicate path smoothing, grid tick generation, and mouse interaction math.
 */

export interface SvgPoint {
  x: number;
  y: number;
}

/**
 * Builds a smooth SVG cubic Bézier path (`C` commands) through a sequence of points.
 */
export function buildSmoothSvgPath(coords: SvgPoint[]): string {
  if (coords.length === 0) return '';
  if (coords.length === 1) return `M ${coords[0].x} ${coords[0].y}`;

  let d = `M ${coords[0].x} ${coords[0].y}`;
  for (let i = 1; i < coords.length; i++) {
    const prev = coords[i - 1];
    const curr = coords[i];
    const cx1 = prev.x + (curr.x - prev.x) / 2;
    const cy1 = prev.y;
    const cx2 = prev.x + (curr.x - prev.x) / 2;
    const cy2 = curr.y;
    d += ` C ${cx1} ${cy1}, ${cx2} ${cy2}, ${curr.x} ${curr.y}`;
  }
  return d;
}

/**
 * Calculates the zero-based index of the nearest data point to the cursor X position.
 */
export function calculateNearestXIndex(
  mouseSvgX: number,
  paddingLeft: number,
  stepX: number,
  itemCount: number
): number {
  if (itemCount <= 1) return 0;
  const index = Math.round((mouseSvgX - paddingLeft) / stepX);
  return Math.max(0, Math.min(itemCount - 1, index));
}

/**
 * Calculates the exact X pixel coordinate for an index on an SVG timeline.
 */
export function calculateXPosition(
  index: number,
  paddingLeft: number,
  stepX: number,
  itemCount: number,
  chartWidth: number
): number {
  return paddingLeft + (itemCount === 1 ? chartWidth / 2 : index * stepX);
}
