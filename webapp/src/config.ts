/**
 * Centralized application configuration for the PyTransDatRO Web Application.
 *
 * Exposes runtime settings loaded from Vite environment variables (VITE_*)
 * with resilient fallbacks for local and production deployment modes.
 */

declare const __APP_VERSION__: string | undefined;

export interface AppConfig {
  /** Application display title. Defaults to 'PyTransDatRO'. */
  appName: string;
  /** Application semantic version (injected at build time from package.json). */
  appVersion: string;
  /** Base URL for API requests. Defaults to empty string (same-origin / Vite dev proxy). */
  apiBaseUrl: string;
  /** Default geographic map center: Romania centroid [Latitude, Longitude]. */
  defaultMapCenter: [number, number];
  /** Default initial map zoom level. */
  defaultMapZoom: number;
  /** Point volume warning threshold triggering sample/render confirmation dialog. */
  batchWarningThreshold: number;
  /** Number of coordinates to sample when the user selects 'sample' in the volume warning dialog. */
  batchSampleSize: number;
  /** Maximum allowable coordinates in an in-memory batch conversion. */
  maxBatchPoints: number;
  /** Maximum allowable coordinates in a streaming file conversion. */
  maxFilePoints: number;
  /** Maximum allowable file upload size in megabytes. */
  maxFileSizeMb: number;
}

export const APP_CONFIG: AppConfig = {
  appName: import.meta.env.VITE_APP_NAME || 'PyTransDatRO',
  appVersion: typeof __APP_VERSION__ !== 'undefined' ? __APP_VERSION__ : (import.meta.env.VITE_APP_VERSION || '1.0.0-beta'),
  apiBaseUrl: import.meta.env.VITE_API_BASE_URL || '',
  defaultMapCenter: [45.9432, 24.9668],
  defaultMapZoom: 7,
  batchWarningThreshold: Number(import.meta.env.VITE_BATCH_WARNING_THRESHOLD) || 2000,
  batchSampleSize: Number(import.meta.env.VITE_BATCH_SAMPLE_SIZE) || 1000,
  maxBatchPoints: 20000,
  maxFilePoints: 100000,
  maxFileSizeMb: 50,
};

/** True only in builds made with VITE_E2E=true (Playwright harness); exposes `window.__e2e`. */
export const E2E_ENABLED: boolean = import.meta.env.VITE_E2E === 'true';
