import type { MapViewSnapshot } from './components/mapViewSnapshot';

declare global {
  interface Window {
    /** Read-only e2e handle; only assigned in builds made with VITE_E2E=true. */
    __e2e?: {
      mapView(): MapViewSnapshot;
      fireMapClick?(lat: number, lng: number): void;
    };
  }
}
