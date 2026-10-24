export type TransformationOp = 'Stereo70ToETRS89' | 'ETRS89ToStereo70';
export type AngleUnit = 'radians' | 'degrees' | 'dms';

export interface PointTransformRequest {
  op: TransformationOp;
  coos: number[];
  unit?: AngleUnit;
}

export interface PointTransformResponse {
  coos: number[];
  warning: string | null;
}

export interface BatchTransformRequest {
  op: TransformationOp;
  points: number[][];
  unit?: AngleUnit;
}

export interface BatchTransformResponse {
  results: PointTransformResponse[];
  count: number;
}

export interface SpgFileInfo {
  name: string;
  size_bytes: number;
  size_formatted: string;
  sha256: string;
}

export interface GridMetadataResponse {
  _generator?: string;
  _generated_at?: string;
  file_info: SpgFileInfo;
  metadata: Record<string, any>;
  params: Record<string, any>;
  grids: Record<string, any>;
  [key: string]: any;
}

export interface GridInfoResponse {
  grid_file: string;
  bounds: {
    min_northing: number;
    min_easting: number;
    max_northing: number;
    max_easting: number;
  };
  crs: {
    projected: string;
    geographic: string;
  };
  quasigeoid_model: string;
  metadata?: Record<string, any>;
  params?: Record<string, any>;
  grids?: Record<string, any>;
  file_info?: SpgFileInfo;
  [key: string]: any;
}

export interface TelemetryStatsResponse {
  summary: {
    total_calls: number;
    total_points: number;
  };
  daily_trend: Array<{
    date: string;
    source: string;
    calls: number;
    points: number;
  }>;
}

export interface SourceItem {
  id: number;
  name: string;
  description?: string | null;
}

export interface TimelineSeriesMetrics {
  web_map?: number;
  rest_api?: number;
  desktop_cli?: number;
  unknown?: number;
  st70_to_etrs?: number;
  etrs_to_st70?: number;
  dim_2d?: number;
  dim_3d?: number;
  total: number;
  sources?: Record<string, number>;
  [key: string]: number | Record<string, number> | undefined;
}

export interface TimelineSeriesItem {
  date: string;
  calls: TimelineSeriesMetrics;
  points: TimelineSeriesMetrics;
}

export interface TimelineSummary {
  total_calls: number;
  total_points: number;
  avg_points_per_call: number;
  from_date: string | null;
  to_date: string | null;
  effective_from_date?: string | null;
  effective_to_date?: string | null;
  days_count: number | null;
  period_count?: number | null;
}

export interface DirectionComparisonSummary {
  direction: string;
  total_points: number;
  status: string;
  max_abs_error: Record<string, number>;
  mean_abs_error: Record<string, number>;
  flags: Record<string, Record<string, number>>;
}

export interface GridComparisonResponse {
  total_points: number;
  overall_status: string;
  etrs89_to_st70: DirectionComparisonSummary;
  st70_to_etrs89: DirectionComparisonSummary;
}

export interface TimelineResponse {
  summary: TimelineSummary;
  series: TimelineSeriesItem[];
  available_sources?: SourceItem[];
}

export interface TimelineQueryParams {
  from_date?: string;
  to_date?: string;
  granularity?: 'daily' | 'weekly' | 'monthly';
  source?: number;
  operation?: string;
  is_3d?: number;
}

export interface SpatialCellItem {
  lat: number;
  lon: number;
  hits: number;
}

export interface QuadrantItem {
  value: number;
  percentage: number;
}

export interface QuadrantSummary {
  center: { lat: number; lon: number };
  NW: QuadrantItem;
  NE: QuadrantItem;
  SW: QuadrantItem;
  SE: QuadrantItem;
}

export interface SpatialBounds {
  min_lat: number;
  max_lat: number;
  min_lon: number;
  max_lon: number;
}

export interface SpatialMatrixSummary {
  active_cells_count: number;
  total_metric_value: number;
  max_cell_hits: number;
  metric: 'points' | 'calls';
  from_date?: string | null;
  to_date?: string | null;
  effective_from_date?: string | null;
  effective_to_date?: string | null;
  quadrants: QuadrantSummary;
}

export interface SpatialMatrixResponse {
  resolution: number;
  bounds: SpatialBounds;
  summary: SpatialMatrixSummary;
  cells: SpatialCellItem[];
}


export interface LegacyStatsResponse {
  total_points: number;
  app_st70: number;
  app_st30: number;
  service_st70: number;
  service_st30: number;
  from_date: string | null;
  to_date: string | null;
  recorded_days: number;
}

export interface LegacyTimelineMetrics {
  app_st70: number;
  app_st30: number;
  service_st70: number;
  service_st30: number;
  app_total: number;
  service_total: number;
  st70_total: number;
  st30_total: number;
  total: number;
}

export interface LegacyTimelineSeriesItem {
  date: string;
  points: LegacyTimelineMetrics;
}

export interface LegacyTimelineSummary {
  total_points: number;
  from_date: string | null;
  to_date: string | null;
  days_count: number | null;
  period_count: number | null;
}

export interface LegacyTimelineResponse {
  summary: LegacyTimelineSummary;
  series: LegacyTimelineSeriesItem[];
  bounds: {
    min_date?: string;
    max_date?: string;
  };
}

export interface LegacyTimelineQueryParams {
  from_date?: string;
  to_date?: string;
  granularity?: 'daily' | 'weekly' | 'monthly';
  source?: string;
  operation?: string;
}
export interface SpatialMatrixQueryParams {
  from_date?: string;
  to_date?: string;
  metric?: 'points' | 'calls';
  sources?: string;
  operation?: number;
  is_3d?: number;
}

import { APP_CONFIG } from '../config';

const COMMON_HEADERS = {
  'Content-Type': 'application/json',
  'X-Client-Source': 'web_app',
};

const apiUrl = (path: string): string => `${APP_CONFIG.apiBaseUrl}${path}`;

export const apiClient = {
  async transformPoint(req: PointTransformRequest): Promise<PointTransformResponse> {
    const res = await fetch(apiUrl('/api/v1/transform/point'), {
      method: 'POST',
      headers: COMMON_HEADERS,
      body: JSON.stringify(req),
    });
    if (!res.ok) {
      throw new Error(`API Error: ${res.status} ${res.statusText}`);
    }
    return res.json();
  },

  async transformBatch(req: BatchTransformRequest): Promise<BatchTransformResponse> {
    const res = await fetch(apiUrl('/api/v1/transform/batch'), {
      method: 'POST',
      headers: COMMON_HEADERS,
      body: JSON.stringify(req),
    });
    if (!res.ok) {
      throw new Error(`API Error: ${res.status} ${res.statusText}`);
    }
    return res.json();
  },

  async transformFile(
    file: File,
    op: TransformationOp,
    unit: AngleUnit = 'degrees',
    delimiter: string = 'auto',
    skipHeader: boolean = false
  ): Promise<Response> {
    const formData = new FormData();
    formData.append('file', file);
    formData.append('op', op);
    formData.append('unit', unit);
    formData.append('delimiter', delimiter);
    formData.append('skip_header', String(skipHeader));

    const res = await fetch(apiUrl('/api/v1/transform/file'), {
      method: 'POST',
      headers: {
        'X-Client-Source': 'web_app',
      },
      body: formData,
    });
    if (!res.ok) {
      throw new Error(`Upload Error: ${res.status} ${res.statusText}`);
    }
    return res;
  },

  async getGridInfo(): Promise<GridInfoResponse> {
    const res = await fetch(apiUrl('/api/v1/grid/info'), {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Grid Info Error: ${res.status}`);
    }
    return res.json();
  },

  async getGridMetadata(): Promise<GridMetadataResponse> {
    const res = await fetch(apiUrl('/api/v1/grid/metadata'), {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Grid Metadata Error: ${res.status}`);
    }
    return res.json();
  },

  async getGridComparison(): Promise<GridComparisonResponse> {
    const res = await fetch(apiUrl('/api/v1/grid/comparison'), {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Grid Comparison Error: ${res.status}`);
    }
    return res.json();
  },

  async getTelemetryStats(): Promise<TelemetryStatsResponse> {
    const res = await fetch(apiUrl('/api/v1/telemetry/stats'), {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Telemetry Error: ${res.status}`);
    }
    return res.json();
  },

  async getTelemetryTimeline(params?: TimelineQueryParams): Promise<TimelineResponse> {
    const query = new URLSearchParams();
    if (params?.from_date) query.set('from_date', params.from_date);
    if (params?.to_date) query.set('to_date', params.to_date);
    if (params?.granularity) query.set('granularity', params.granularity);
    if (params?.source !== undefined) query.set('source', String(params.source));
    if (params?.operation) query.set('operation', params.operation);
    if (params?.is_3d !== undefined) query.set('is_3d', String(params.is_3d));

    const qs = query.toString();
    const url = apiUrl(`/api/v1/telemetry/timeline${qs ? `?${qs}` : ''}`);
    const res = await fetch(url, {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Timeline Error: ${res.status}`);
    }
    return res.json();
  },

  async getSpatialMatrix(params?: SpatialMatrixQueryParams): Promise<SpatialMatrixResponse> {
    const query = new URLSearchParams();
    if (params?.from_date) query.set('from_date', params.from_date);
    if (params?.to_date) query.set('to_date', params.to_date);
    if (params?.metric) query.set('metric', params.metric);
    if (params?.sources) query.set('sources', params.sources);
    if (params?.operation !== undefined) query.set('operation', String(params.operation));
    if (params?.is_3d !== undefined) query.set('is_3d', String(params.is_3d));

    const qs = query.toString();
    const url = apiUrl(`/api/v1/telemetry/spatial-matrix${qs ? `?${qs}` : ''}`);
    const res = await fetch(url, {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Spatial Matrix Error: ${res.status}`);
    }
    return res.json();
  },


  async getLegacyStats(): Promise<LegacyStatsResponse> {
    const res = await fetch(apiUrl('/api/v1/telemetry/legacy/stats'), {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Legacy Stats Error: ${res.status}`);
    }
    return res.json();
  },

  async getLegacyTimeline(params?: LegacyTimelineQueryParams): Promise<LegacyTimelineResponse> {
    const query = new URLSearchParams();
    if (params?.from_date) query.set('from_date', params.from_date);
    if (params?.to_date) query.set('to_date', params.to_date);
    if (params?.granularity) query.set('granularity', params.granularity);
    if (params?.source) query.set('source', params.source);
    if (params?.operation) query.set('operation', params.operation);

    const qs = query.toString();
    const url = apiUrl(`/api/v1/telemetry/legacy/timeline${qs ? `?${qs}` : ''}`);
    const res = await fetch(url, {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Legacy Timeline Error: ${res.status}`);
    }
    return res.json();
  },
  async checkHealth(): Promise<{ status: string; grid_loaded: boolean }> {
    const res = await fetch(apiUrl('/health'));
    return res.json();
  },

  async submitFeedback(payload: {
    report_type: 'discrepancy' | 'grid_update' | 'bug' | 'feature' | 'api_down';
    description: string;
    client_timestamp?: string;
    app_version?: string;
    active_grid?: string;
    honeypot?: string;
  }): Promise<{ status: 'success' | 'logged'; delivered_email: boolean; message: string }> {
    const res = await fetch(apiUrl('/api/v1/feedback'), {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'X-Client-Source': 'web_app',
      },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      throw new Error(`Feedback submission failed: ${res.status}`);
    }
    return res.json();
  },
};
