import { formatTelemetryNumber, formatTelemetryCompact } from '../utils/numberFormat';
import { t } from '../i18n';
import type { LegacyTimelineSeriesItem } from '../utils/apiClient';

export interface LegacySeriesConfig {
  key: string;
  labelKey: string;
  defaultColor: string;
}

export const LEGACY_SERIES_DEFS: Record<string, LegacySeriesConfig> = {
  app_st70: {
    key: 'app_st70',
    labelKey: 'stats.legacySeriesAppSt70',
    defaultColor: '#34d399', // Emerald
  },
  app_st30: {
    key: 'app_st30',
    labelKey: 'stats.legacySeriesAppSt30',
    defaultColor: '#059669', // Deep Emerald
  },
  service_st70: {
    key: 'service_st70',
    labelKey: 'stats.legacySeriesServiceSt70',
    defaultColor: '#818cf8', // Bright Indigo
  },
  service_st30: {
    key: 'service_st30',
    labelKey: 'stats.legacySeriesServiceSt30',
    defaultColor: '#4f46e5', // Deep Indigo
  },
};

export const LEGACY_TOTAL_CONFIG = {
  key: 'total',
  labelKey: 'stats.seriesTotal',
  defaultColor: '#0ea5e9', // Vibrant Sky Cyan dashed line
};

const LEGACY_COLORS_STORAGE_KEY = 'pytransdat_legacy_chart_colors';

export interface LegacyChartOptions {
  activeKeys: Set<string>; // subset of ['app_st70', 'app_st30', 'service_st70', 'service_st30']
  showTotal: boolean;
}

export class LegacyTimelineChart {
  private container: HTMLElement;
  private data: LegacyTimelineSeriesItem[] = [];
  private customColors: Record<string, string> = {};
  private options: LegacyChartOptions = {
    activeKeys: new Set(['app_st70', 'service_st70']),
    showTotal: true,
  };

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'timeline-chart-wrapper';
    this.loadCustomColors();

    window.addEventListener('themechange', () => {
      this.render();
    });
  }

  private loadCustomColors(): void {
    try {
      const raw = localStorage.getItem(LEGACY_COLORS_STORAGE_KEY);
      if (raw) {
        const parsed = JSON.parse(raw);
        if (parsed && typeof parsed === 'object') {
          this.customColors = parsed as Record<string, string>;
        }
      }
    } catch (err) {
      console.warn('Failed to load custom legacy colors from localStorage:', err);
    }
  }

  public getLineColor(seriesKey: string, defaultColor: string): string {
    return this.customColors[seriesKey] || defaultColor;
  }

  public setLineColor(seriesKey: string, color: string): void {
    this.customColors[seriesKey] = color;
    try {
      localStorage.setItem(LEGACY_COLORS_STORAGE_KEY, JSON.stringify(this.customColors));
    } catch (err) {
      console.warn('Failed to save legacy colors to localStorage:', err);
    }
    this.render();
  }

  public getTotalColor(): string {
    return this.getLineColor(LEGACY_TOTAL_CONFIG.key, LEGACY_TOTAL_CONFIG.defaultColor);
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public update(data: LegacyTimelineSeriesItem[], options?: Partial<LegacyChartOptions>): void {
    this.data = data;
    if (options) {
      if (options.activeKeys) this.options.activeKeys = options.activeKeys;
      if (options.showTotal !== undefined) this.options.showTotal = options.showTotal;
    }
    this.render();
  }

  public render(): void {
    if (!this.data || this.data.length === 0) {
      this.container.innerHTML = `
        <div class="timeline-empty-state">
          <span>${t('stats.noLegacyDataFound')}</span>
        </div>
      `;
      return;
    }

    const width = 880;
    const height = 340;
    const padding = { top: 30, right: 80, bottom: 45, left: 70 };
    const chartWidth = width - padding.left - padding.right;
    const chartHeight = height - padding.top - padding.bottom;

    const n = this.data.length;
    const stepX = chartWidth / Math.max(n - 1, 1);

    const activeKeysList = Array.from(this.options.activeKeys);

    // Calculate maximum points across selected series and total
    let maxPoints = 0;
    for (const item of this.data) {
      let bucketTotal = 0;
      for (const k of activeKeysList) {
        const val = (item.points as any)[k] ?? 0;
        bucketTotal += val;
        if (val > maxPoints) maxPoints = val;
      }
      if (this.options.showTotal && bucketTotal > maxPoints) {
        maxPoints = bucketTotal;
      }
    }

    maxPoints = Math.max(maxPoints, 100);
    maxPoints = Math.ceil(maxPoints * 1.15);

    // Y Grid lines
    const yTicks = 4;
    let gridLinesSvg = '';
    let axisLabelsSvg = '';

    for (let i = 0; i <= yTicks; i++) {
      const ratio = i / yTicks;
      const y = padding.top + chartHeight - ratio * chartHeight;
      const pointTickVal = Math.round(ratio * maxPoints);

      gridLinesSvg += `<line x1="${padding.left}" y1="${y}" x2="${width - padding.right}" y2="${y}" class="chart-grid-line" />`;
      axisLabelsSvg += `<text x="${padding.left - 8}" y="${y + 4}" text-anchor="end" class="chart-axis-text">${this.formatNumber(pointTickVal)}</text>`;
    }

    // X axis ticks
    let xAxisLabelsSvg = '';
    const labelInterval = Math.max(1, Math.floor(n / 6));

    for (let i = 0; i < n; i += labelInterval) {
      const item = this.data[i];
      const x = padding.left + (n === 1 ? chartWidth / 2 : i * stepX);
      xAxisLabelsSvg += `<text x="${x}" y="${height - 12}" text-anchor="middle" class="chart-axis-text">${item.date}</text>`;
    }
    if ((n - 1) % labelInterval !== 0) {
      const lastItem = this.data[n - 1];
      const lastX = padding.left + chartWidth;
      xAxisLabelsSvg += `<text x="${lastX}" y="${height - 12}" text-anchor="middle" class="chart-axis-text">${lastItem.date}</text>`;
    }

    // Series paths
    let pathsSvg = '';

    const buildSmoothPath = (coords: { x: number; y: number }[]): string => {
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
    };

    for (const key of activeKeysList) {
      const def = LEGACY_SERIES_DEFS[key];
      if (!def) continue;
      const color = this.getLineColor(key, def.defaultColor);
      const coords: { x: number; y: number }[] = [];

      this.data.forEach((item, index) => {
        const x = padding.left + (n === 1 ? chartWidth / 2 : index * stepX);
        const val = (item.points as any)[key] ?? 0;
        const y = padding.top + chartHeight - (val / maxPoints) * chartHeight;
        coords.push({ x, y });
      });

      const pathD = buildSmoothPath(coords);
      pathsSvg += `
        <path 
          d="${pathD}" 
          fill="none" 
          stroke="${color}" 
          stroke-width="2.2" 
          stroke-linecap="round" 
          stroke-linejoin="round" 
          class="chart-series-line line-points"
          data-series="${key}"
        />
      `;
    }

    // Total line
    if (this.options.showTotal && activeKeysList.length > 0) {
      const totalColor = this.getTotalColor();
      const totalCoords: { x: number; y: number }[] = [];

      this.data.forEach((item, index) => {
        const x = padding.left + (n === 1 ? chartWidth / 2 : index * stepX);
        let sum = 0;
        for (const k of activeKeysList) {
          sum += (item.points as any)[k] ?? 0;
        }
        const y = padding.top + chartHeight - (sum / maxPoints) * chartHeight;
        totalCoords.push({ x, y });
      });

      const totalD = buildSmoothPath(totalCoords);
      pathsSvg += `
        <path 
          d="${totalD}" 
          fill="none" 
          stroke="${totalColor}" 
          stroke-width="2.6" 
          stroke-dasharray="7 5" 
          stroke-linecap="round" 
          stroke-linejoin="round" 
          class="chart-total-line total-points" 
          data-series="total-points" 
        />
      `;
    }

    this.container.innerHTML = `
      <div class="timeline-chart-inner">
        <svg viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" class="timeline-svg" id="legacy-timeline-svg">
          <!-- Grid Lines -->
          <g class="chart-grid">${gridLinesSvg}</g>

          <!-- Base Axis Lines -->
          <line x1="${padding.left}" y1="${padding.top + chartHeight}" x2="${width - padding.right}" y2="${padding.top + chartHeight}" class="chart-axis-base" />

          <!-- Series & Total Lines -->
          <g class="chart-lines">${pathsSvg}</g>

          <!-- Axes Text -->
          <g class="chart-labels">${axisLabelsSvg}</g>
          <g class="chart-labels">${xAxisLabelsSvg}</g>

          <!-- Crosshair Guide (hidden by default) -->
          <line id="legacy-crosshair-x" x1="0" y1="${padding.top}" x2="0" y2="${padding.top + chartHeight}" class="chart-crosshair" style="display: none;" />
        </svg>

        <!-- Floating Tooltip -->
        <div id="legacy-chart-tooltip" class="timeline-tooltip" style="display: none;"></div>
      </div>
    `;

    this.attachInteractivity(padding, chartWidth, chartHeight, width, stepX, activeKeysList);
  }

  private attachInteractivity(
    padding: { top: number; right: number; bottom: number; left: number },
    chartWidth: number,
    _chartHeight: number,
    totalWidth: number,
    stepX: number,
    activeKeysList: string[]
  ): void {
    const svg = this.container.querySelector('#legacy-timeline-svg') as SVGSVGElement | null;
    const tooltip = this.container.querySelector('#legacy-chart-tooltip') as HTMLElement | null;
    const crosshair = this.container.querySelector('#legacy-crosshair-x') as SVGLineElement | null;

    if (!svg || !tooltip || !crosshair) return;

    svg.addEventListener('mousemove', (e: MouseEvent) => {
      const rect = svg.getBoundingClientRect();
      const mouseSvgX = ((e.clientX - rect.left) / rect.width) * totalWidth;

      if (mouseSvgX < padding.left || mouseSvgX > totalWidth - padding.right) {
        crosshair.style.display = 'none';
        tooltip.style.display = 'none';
        return;
      }

      const n = this.data.length;
      let nearestIndex = n === 1 ? 0 : Math.round((mouseSvgX - padding.left) / stepX);
      nearestIndex = Math.max(0, Math.min(n - 1, nearestIndex));

      const item = this.data[nearestIndex];
      const targetX = padding.left + (n === 1 ? chartWidth / 2 : nearestIndex * stepX);

      crosshair.setAttribute('x1', String(targetX));
      crosshair.setAttribute('x2', String(targetX));
      crosshair.style.display = 'block';

      // Tooltip rows
      let sourceRowsHtml = '';
      let totalSelected = 0;

      for (const k of activeKeysList) {
        const def = LEGACY_SERIES_DEFS[k];
        if (!def) continue;
        const color = this.getLineColor(k, def.defaultColor);
        const val = (item.points as any)[k] ?? 0;
        totalSelected += val;

        sourceRowsHtml += `
          <div class="tooltip-source-block">
            <div class="tooltip-source-title">
              <span class="tooltip-dot" style="background-color: ${color};"></span>
              <strong>${t(def.labelKey)}</strong>
            </div>
            <div class="tooltip-source-details">
              <span><span class="legend-indicator" style="background-color: ${color};"></span> ${t('stats.metricTypePoints')}: <strong>${formatTelemetryNumber(val)}</strong></span>
            </div>
          </div>
        `;
      }

      let totalRowHtml = '';
      if (this.options.showTotal && activeKeysList.length > 0) {
        const totalColor = this.getTotalColor();
        totalRowHtml = `
          <div class="tooltip-row total-row">
            <div class="tooltip-source-title">
              <span class="tooltip-dot dashed" style="background-color: ${totalColor};"></span>
              <strong>${t('stats.seriesTotal')}</strong>
            </div>
            <div class="tooltip-source-details">
              <span><span class="legend-indicator dashed" style="background-color: ${totalColor};"></span> ${t('stats.metricTypePoints')}: <strong>${formatTelemetryNumber(totalSelected)}</strong></span>
            </div>
          </div>
        `;
      }

      tooltip.innerHTML = `
        <div class="tooltip-header">${item.date}</div>
        <div class="tooltip-body">
          ${sourceRowsHtml}
          ${totalRowHtml}
        </div>
      `;

      tooltip.style.display = 'block';
      const tooltipRect = tooltip.getBoundingClientRect();
      const clientX = e.clientX - rect.left;

      let leftPx = clientX + 15;
      if (leftPx + tooltipRect.width > rect.width) {
        leftPx = clientX - tooltipRect.width - 15;
      }
      tooltip.style.left = `${Math.max(10, leftPx)}px`;
      tooltip.style.top = `${Math.max(10, e.clientY - rect.top - 20)}px`;
    });

    svg.addEventListener('mouseleave', () => {
      crosshair.style.display = 'none';
      tooltip.style.display = 'none';
    });
  }

  private formatNumber(num: number): string {
    return formatTelemetryCompact(num, 1);
  }
}
