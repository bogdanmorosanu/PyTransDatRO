export type TransformationOp = 'Stereo70ToETRS89' | 'ETRS89ToStereo70';
export type AngleUnit = 'radians' | 'degrees';

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

const COMMON_HEADERS = {
  'Content-Type': 'application/json',
  'X-Client-Source': 'web_app',
};

export const apiClient = {
  async transformPoint(req: PointTransformRequest): Promise<PointTransformResponse> {
    const res = await fetch('/api/v1/transform/point', {
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
    const res = await fetch('/api/v1/transform/batch', {
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
    delimiter: string = 'auto'
  ): Promise<Response> {
    const formData = new FormData();
    formData.append('file', file);
    formData.append('op', op);
    formData.append('unit', unit);
    formData.append('delimiter', delimiter);

    const res = await fetch('/api/v1/transform/file', {
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
    const res = await fetch('/api/v1/grid/info', {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Grid Info Error: ${res.status}`);
    }
    return res.json();
  },

  async getTelemetryStats(): Promise<TelemetryStatsResponse> {
    const res = await fetch('/api/v1/telemetry/stats', {
      headers: { 'X-Client-Source': 'web_app' },
    });
    if (!res.ok) {
      throw new Error(`Telemetry Error: ${res.status}`);
    }
    return res.json();
  },

  async checkHealth(): Promise<{ status: string; grid_loaded: boolean }> {
    const res = await fetch('/health');
    return res.json();
  },
};
