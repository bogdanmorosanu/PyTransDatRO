// Types only, with no imports, so the e2e project (tsc -p e2e) can share them without pulling in
// Leaflet, CSS or Vite env typings.

export interface BoundsBox {
  south: number;
  west: number;
  north: number;
  east: number;
}

/** Plain-data snapshot of the map state, exposed to e2e tests only (see `window.__e2e`). */
export interface MapViewSnapshot {
  center: { lat: number; lng: number };
  zoom: number;
  bounds: BoundsBox;
  maxBounds: BoundsBox | null;
  markerCount: number;
  markersByLayer: { point: number; text: number; files: number };
}
