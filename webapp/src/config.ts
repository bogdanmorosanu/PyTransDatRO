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
  /** Source repository URL for help, documentation, and issue reporting. */
  repoUrl: string;
  /** Default geographic map center: Romania centroid [Latitude, Longitude]. */
  defaultMapCenter: [number, number];
  /** Default initial map zoom level. */
  defaultMapZoom: number;
  /** Zoom level applied when centering the map on an individual point. */
  pointZoomLevel: number;
  /** Minimum map zoom constraint. */
  minMapZoom: number;
  /** Maximum map zoom constraint. */
  maxMapZoom: number;
  /** Background health check polling interval in milliseconds. */
  healthPollIntervalMs: number;
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
  appVersion: typeof __APP_VERSION__ !== 'undefined' ? __APP_VERSION__ : (import.meta.env.VITE_APP_VERSION || '1.0.0'),
  apiBaseUrl: (import.meta.env.VITE_API_BASE_URL || '').replace(/\/+$/, ''),
  repoUrl: import.meta.env.VITE_REPO_URL || 'https://github.com/bogdanmorosanu/PyTransDatRO',
  defaultMapCenter: [
    Number(import.meta.env.VITE_DEFAULT_MAP_LAT) || 45.9432,
    Number(import.meta.env.VITE_DEFAULT_MAP_LON) || 24.9668,
  ],
  defaultMapZoom: Number(import.meta.env.VITE_DEFAULT_MAP_ZOOM) || 7,
  pointZoomLevel: Number(import.meta.env.VITE_POINT_ZOOM_LEVEL) || 14,
  minMapZoom: 3,
  maxMapZoom: 18,
  healthPollIntervalMs: Number(import.meta.env.VITE_HEALTH_POLL_MS) || 30000,
  batchWarningThreshold: Number(import.meta.env.VITE_BATCH_WARNING_THRESHOLD) || 2000,
  batchSampleSize: Number(import.meta.env.VITE_BATCH_SAMPLE_SIZE) || 1000,
  maxBatchPoints: Number(import.meta.env.VITE_MAX_BATCH_POINTS) || 20000,
  maxFilePoints: Number(import.meta.env.VITE_MAX_FILE_POINTS) || 100000,
  maxFileSizeMb: Number(import.meta.env.VITE_MAX_FILE_SIZE_MB) || 50,
};

export interface RuntimeLimits {
  max_batch_points?: number;
  max_file_points?: number;
  max_file_upload_bytes?: number;
}

/**
 * Synchronizes client-side limits dynamically from backend GET /health response.
 * Allows server .env configuration to take immediate effect in client UI without rebuild.
 */
export function syncRuntimeLimits(limits: RuntimeLimits | null | undefined): void {
  if (!limits) return;
  if (typeof limits.max_batch_points === 'number' && limits.max_batch_points > 0) {
    APP_CONFIG.maxBatchPoints = limits.max_batch_points;
  }
  if (typeof limits.max_file_points === 'number' && limits.max_file_points > 0) {
    APP_CONFIG.maxFilePoints = limits.max_file_points;
  }
  if (typeof limits.max_file_upload_bytes === 'number' && limits.max_file_upload_bytes > 0) {
    APP_CONFIG.maxFileSizeMb = Math.round(limits.max_file_upload_bytes / (1024 * 1024));
  }
}

/** True only in builds made with VITE_E2E=true (Playwright harness); exposes `window.__e2e`. */
export const E2E_ENABLED: boolean = import.meta.env.VITE_E2E === 'true';
