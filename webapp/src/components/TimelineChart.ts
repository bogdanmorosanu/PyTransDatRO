import { formatTelemetryNumber, formatTelemetryCompact } from '../utils/numberFormat';
import { t } from '../i18n';
import type { TimelineSeriesItem, SourceItem } from '../utils/apiClient';

export interface SourceDisplayConfig {
  key: string;
  label: string;
  callsColor: string;
  pointsColor: string;
}

export const KNOWN_SOURCE_COLORS: Record<string, { calls: string; points: string; labelKey: string }> = {
  rest_api: {
    calls: '#818cf8',    // Bright Indigo (calls)
    points: '#4f46e5',   // Deep Indigo (points)
    labelKey: 'stats.sourceRestApi',
  },
  web_app: {
    calls: '#34d399',    // Emerald (calls)
    points: '#059669',   // Deep Emerald (points)
    labelKey: 'stats.sourceWebMap',
  },
  web_map: {
    calls: '#34d399',    // Emerald (calls)
    points: '#059669',   // Deep Emerald (points)
    labelKey: 'stats.sourceWebMap',
  },
  desktop_cli: {
    calls: '#fb7185',    // Coral Rose (calls — high-contrast on Dark & Light)
    points: '#e11d48',   // Deep Crimson Rose (points — high-contrast on Dark & Light)
    labelKey: 'stats.sourceDesktopCli',
  },
  unknown: {
    calls: '#94a3b8',    // Medium Slate
    points: '#64748b',   // Deep Slate
    labelKey: 'stats.sourceUnknown',
  },
};

export const TIMELINE_COLOR_PRESETS = [
  '#818cf8', // Bright Indigo
  '#4f46e5', // Deep Indigo
  '#34d399', // Emerald
  '#059669', // Deep Emerald
  '#fb7185', // Coral Rose
  '#e11d48', // Crimson Rose
  '#c084fc', // Bright Purple
  '#9333ea', // Deep Purple
  '#38bdf8', // Sky Cyan
  '#f97316', // Vibrant Orange
];

const TIMELINE_COLORS_STORAGE_KEY = 'pytransdat_timeline_colors';

// Distinct fallback palette pairs for dynamic sources beyond the known ones (calls lighter, points darker)
const FALLBACK_PALETTE_PAIRS = [
  { calls: '#f472b6', points: '#db2777' }, // Pink
  { calls: '#c084fc', points: '#9333ea' }, // Purple
  { calls: '#2dd4bf', points: '#0d9488' }, // Teal
  { calls: '#fb923c', points: '#ea580c' }, // Bright Orange
  { calls: '#38bdf8', points: '#0284c7' }, // Cyan
];

export const TOTAL_SERIES_CONFIG = {
  key: 'total',
  labelKey: 'stats.seriesTotal',
  callsColor: 'var(--text-main, #f8fafc)',  // Theme-adaptive dashed line for Total Invocations
  pointsColor: '#0ea5e9', // Vibrant Sky Cyan dashed line for Total Coordinates
};

export type TimelineMetricFilter = 'both' | 'points' | 'calls';

export interface TimelineChartOptions {
  activeSourceKeys: Set<string>;
  availableSources: SourceItem[];
  showTotal: boolean;
  metricFilter: TimelineMetricFilter;
}

export class TimelineChart {
  private container: HTMLElement;
  private data: TimelineSeriesItem[] = [];
  private customColors: Record<string, string> = {};
  private options: TimelineChartOptions = {
    activeSourceKeys: new Set(['rest_api', 'web_map', 'desktop_cli']),
    availableSources: [],
    showTotal: true,
    metricFilter: 'both',
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
      const raw = localStorage.getItem(TIMELINE_COLORS_STORAGE_KEY);
      if (raw) {
        const parsed = JSON.parse(raw);
        if (parsed && typeof parsed === 'object') {
          this.customColors = parsed as Record<string, string>;
        }
      }
    } catch (err) {
      console.warn('Failed to load custom timeline colors from localStorage:', err);
    }
  }

  public getLineColor(seriesKey: string, defaultColor: string): string {
    return this.customColors[seriesKey] || defaultColor;
  }

  public setLineColor(seriesKey: string, color: string): void {
    this.customColors[seriesKey] = color;
    try {
      localStorage.setItem(TIMELINE_COLORS_STORAGE_KEY, JSON.stringify(this.customColors));
    } catch (err) {
      console.warn('Failed to save custom timeline colors to localStorage:', err);
    }
    this.render();
  }

  public getTotalSeriesColors(): { callsColor: string; pointsColor: string } {
    return {
      callsColor: this.getLineColor('total:calls', TOTAL_SERIES_CONFIG.callsColor),
      pointsColor: this.getLineColor('total:points', TOTAL_SERIES_CONFIG.pointsColor),
    };
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public update(data: TimelineSeriesItem[], options?: Partial<TimelineChartOptions>): void {
    this.data = data;
    if (options) {
      if (options.activeSourceKeys) this.options.activeSourceKeys = options.activeSourceKeys;
      if (options.availableSources) this.options.availableSources = options.availableSources;
      if (options.showTotal !== undefined) this.options.showTotal = options.showTotal;
      if (options.metricFilter) this.options.metricFilter = options.metricFilter;
    }
    this.render();
  }

  public getSourceConfig(sourceName: string, index: number = 0): SourceDisplayConfig {
    if (KNOWN_SOURCE_COLORS[sourceName]) {
      const known = KNOWN_SOURCE_COLORS[sourceName];
      return {
        key: sourceName,
        label: t(known.labelKey),
        callsColor: this.getLineColor(`${sourceName}:calls`, known.calls),
        pointsColor: this.getLineColor(`${sourceName}:points`, known.points),
      };
    }

    const pair = FALLBACK_PALETTE_PAIRS[index % FALLBACK_PALETTE_PAIRS.length];
    // Capitalize source name as label
    const formatted = sourceName
      .split('_')
      .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
      .join(' ');

    return {
      key: sourceName,
      label: formatted,
      callsColor: this.getLineColor(`${sourceName}:calls`, pair.calls),
      pointsColor: this.getLineColor(`${sourceName}:points`, pair.points),
    };
  }

  private render(): void {
    if (!this.data || this.data.length === 0) {
      this.container.innerHTML = `
        <div class="timeline-empty-state">
          <span>${t('stats.noDataFound')}</span>
        </div>
      `;
      return;
    }

    const metricMode = this.options.metricFilter;
    const showCalls = metricMode === 'both' || metricMode === 'calls';
    const showPoints = metricMode === 'both' || metricMode === 'points';

    // Chart dimensions
    const width = 880;
    const height = 350;
    const padding = { 
      top: 30, 
      right: showPoints ? 80 : 30, 
      bottom: 45, 
      left: showCalls ? 70 : 30 
    };
    const chartWidth = width - padding.left - padding.right;
    const chartHeight = height - padding.top - padding.bottom;

    const n = this.data.length;
    const stepX = chartWidth / Math.max(n - 1, 1);

    // Build active source configs
    const activeConfigs: SourceDisplayConfig[] = [];
    let idx = 0;
    for (const srcKey of this.options.activeSourceKeys) {
      activeConfigs.push(this.getSourceConfig(srcKey, idx++));
    }

    // Compute maximums for axes
    let maxCalls = 0;
    let maxPoints = 0;

    for (const item of this.data) {
      let bucketSelectedCalls = 0;
      let bucketSelectedPoints = 0;

      for (const s of activeConfigs) {
        const cVal = item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
        const pVal = item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;
        bucketSelectedCalls += cVal;
        bucketSelectedPoints += pVal;
        if (cVal > maxCalls) maxCalls = cVal;
        if (pVal > maxPoints) maxPoints = pVal;
      }

      // Total line strictly sums selected sources
      if (bucketSelectedCalls > maxCalls) maxCalls = bucketSelectedCalls;
      if (bucketSelectedPoints > maxPoints) maxPoints = bucketSelectedPoints;
    }

    maxCalls = Math.max(maxCalls, 10);
    maxPoints = Math.max(maxPoints, 100);
    maxCalls = Math.ceil(maxCalls * 1.15);
    maxPoints = Math.ceil(maxPoints * 1.15);

    // Y Grid lines (4 intervals)
    const yTicks = 4;
    let gridLinesSvg = '';
    let leftAxisLabelsSvg = '';
    let rightAxisLabelsSvg = '';

    for (let i = 0; i <= yTicks; i++) {
      const ratio = i / yTicks;
      const y = padding.top + chartHeight - ratio * chartHeight;
      const callTickVal = Math.round(ratio * maxCalls);
      const pointTickVal = Math.round(ratio * maxPoints);

      gridLinesSvg += `<line x1="${padding.left}" y1="${y}" x2="${width - padding.right}" y2="${y}" class="chart-grid-line" />`;
      if (showCalls) {
        leftAxisLabelsSvg += `<text x="${padding.left - 8}" y="${y + 4}" text-anchor="end" class="chart-axis-text">${this.formatNumber(callTickVal)}</text>`;
      }
      if (showPoints) {
        rightAxisLabelsSvg += `<text x="${width - padding.right + 8}" y="${y + 4}" text-anchor="start" class="chart-axis-text">${this.formatCompact(pointTickVal)}</text>`;
      }
    }

    // Add Axis Title Labels (Neutral grey matching X-axis)
    if (showCalls) {
      leftAxisLabelsSvg += `<text x="${padding.left - 8}" y="${padding.top - 10}" text-anchor="end" class="chart-axis-title calls-title">${t('stats.metricTypeCalls')}</text>`;
    }
    if (showPoints) {
      rightAxisLabelsSvg += `<text x="${width - padding.right + 8}" y="${padding.top - 10}" text-anchor="start" class="chart-axis-title points-title">${t('stats.metricTypePoints')}</text>`;
    }

    // Build series paths
    let pathsSvg = '';

    // 1. Individual source lines based on toggles
    for (const s of activeConfigs) {
      const callCoords: { x: number; y: number }[] = [];
      const pointCoords: { x: number; y: number }[] = [];

      this.data.forEach((item, index) => {
        const x = padding.left + (n === 1 ? chartWidth / 2 : index * stepX);
        const cVal = item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
        const pVal = item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;

        const yC = padding.top + chartHeight - (cVal / maxCalls) * chartHeight;
        const yP = padding.top + chartHeight - (pVal / maxPoints) * chartHeight;

        callCoords.push({ x, y: yC });
        pointCoords.push({ x, y: yP });
      });

      if (showCalls) {
        const callPathD = this.buildSmoothPath(callCoords);
        pathsSvg += `
          <path 
            d="${callPathD}" 
            fill="none" 
            stroke="${s.callsColor}" 
            stroke-width="2.2" 
            stroke-linecap="round" 
            stroke-linejoin="round"
            class="chart-series-line line-calls"
            data-series="${s.key}-calls"
          />
        `;
      }

      if (showPoints) {
        const pointPathD = this.buildSmoothPath(pointCoords);
        pathsSvg += `
          <path 
            d="${pointPathD}" 
            fill="none" 
            stroke="${s.pointsColor}" 
            stroke-width="2.2" 
            stroke-linecap="round" 
            stroke-linejoin="round"
            class="chart-series-line line-points"
            data-series="${s.key}-points"
          />
        `;
      }
    }

    // 2. Persistent Total Lines (calculated strictly from selected sources)
    const totalColors = this.getTotalSeriesColors();
    if (this.options.showTotal && activeConfigs.length > 0) {
      const totalCallCoords: { x: number; y: number }[] = [];
      const totalPointCoords: { x: number; y: number }[] = [];

      this.data.forEach((item, index) => {
        const x = padding.left + (n === 1 ? chartWidth / 2 : index * stepX);
        let cVal = 0;
        let pVal = 0;
        for (const s of activeConfigs) {
          cVal += item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
          pVal += item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;
        }

        const yC = padding.top + chartHeight - (cVal / maxCalls) * chartHeight;
        const yP = padding.top + chartHeight - (pVal / maxPoints) * chartHeight;

        totalCallCoords.push({ x, y: yC });
        totalPointCoords.push({ x, y: yP });
      });

      if (showCalls) {
        const totalCallPathD = this.buildSmoothPath(totalCallCoords);
        pathsSvg += `
          <path 
            d="${totalCallPathD}" 
            fill="none" 
            stroke="${totalColors.callsColor}" 
            stroke-width="2.6" 
            stroke-dasharray="7 5"
            stroke-linecap="round" 
            stroke-linejoin="round"
            class="chart-total-line total-calls"
            data-series="total-calls"
          />
        `;
      }

      if (showPoints) {
        const totalPointPathD = this.buildSmoothPath(totalPointCoords);
        pathsSvg += `
          <path 
            d="${totalPointPathD}" 
            fill="none" 
            stroke="${totalColors.pointsColor}" 
            stroke-width="2.6" 
            stroke-dasharray="7 5"
            stroke-linecap="round" 
            stroke-linejoin="round"
            class="chart-total-line total-points"
            data-series="total-points"
          />
        `;
      }
    }

    // X Axis Dates
    let xAxisLabelsSvg = '';
    const maxLabels = Math.min(8, n);
    const labelInterval = Math.max(1, Math.floor(n / maxLabels));

    for (let i = 0; i < n; i += labelInterval) {
      const item = this.data[i];
      const x = padding.left + (n === 1 ? chartWidth / 2 : i * stepX);
      const dateLabel = this.formatDateLabel(item.date);
      xAxisLabelsSvg += `<text x="${x}" y="${height - 12}" text-anchor="middle" class="chart-axis-text">${dateLabel}</text>`;
    }
    if ((n - 1) % labelInterval !== 0) {
      const lastItem = this.data[n - 1];
      const lastX = padding.left + chartWidth;
      const dateLabel = this.formatDateLabel(lastItem.date);
      xAxisLabelsSvg += `<text x="${lastX}" y="${height - 12}" text-anchor="middle" class="chart-axis-text">${dateLabel}</text>`;
    }

    this.container.innerHTML = `
      <div class="timeline-chart-inner">
        <svg viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" class="timeline-svg" id="timeline-svg">
          <!-- Grid Lines -->
          <g class="chart-grid">${gridLinesSvg}</g>

          <!-- Base Axis Lines -->
          <line x1="${padding.left}" y1="${padding.top + chartHeight}" x2="${width - padding.right}" y2="${padding.top + chartHeight}" class="chart-axis-base" />

          <!-- Series & Total Lines -->
          <g class="chart-lines">${pathsSvg}</g>

          <!-- Axes Text -->
          <g class="chart-labels">${leftAxisLabelsSvg}</g>
          <g class="chart-labels">${rightAxisLabelsSvg}</g>
          <g class="chart-labels">${xAxisLabelsSvg}</g>

          <!-- Crosshair Guide (hidden by default) -->
          <line id="crosshair-x" x1="0" y1="${padding.top}" x2="0" y2="${padding.top + chartHeight}" class="chart-crosshair" style="display: none;" />
        </svg>

        <!-- Floating Tooltip -->
        <div id="chart-tooltip" class="timeline-tooltip" style="display: none;"></div>
      </div>
    `;

    this.attachInteractivity(padding, chartWidth, chartHeight, width, stepX, activeConfigs);
  }

  private buildSmoothPath(coords: { x: number; y: number }[]): string {
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

  private attachInteractivity(
    padding: { top: number; right: number; bottom: number; left: number },
    chartWidth: number,
    _chartHeight: number,
    totalWidth: number,
    stepX: number,
    activeConfigs: SourceDisplayConfig[]
  ): void {
    const svg = this.container.querySelector('#timeline-svg') as SVGSVGElement | null;
    const tooltip = this.container.querySelector('#chart-tooltip') as HTMLElement | null;
    const crosshair = this.container.querySelector('#crosshair-x') as SVGLineElement | null;

    if (!svg || !tooltip || !crosshair) return;

    const metricMode = this.options.metricFilter;
    const showCalls = metricMode === 'both' || metricMode === 'calls';
    const showPoints = metricMode === 'both' || metricMode === 'points';
    const totalColors = this.getTotalSeriesColors();

    const hideCrosshairAndTooltip = () => {
      crosshair.style.display = 'none';
      tooltip.style.display = 'none';
    };

    const updateCrosshairAndTooltip = (clientX: number, clientY: number) => {
      const rect = svg.getBoundingClientRect();
      const mouseSvgX = ((clientX - rect.left) / rect.width) * totalWidth;

      if (mouseSvgX < padding.left || mouseSvgX > totalWidth - padding.right) {
        hideCrosshairAndTooltip();
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

      // Tooltip rows: Sources breakdown
      let sourceRowsHtml = '';
      for (const s of activeConfigs) {
        const cVal = item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
        const pVal = item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;
        
        let detailsHtml = '';
        if (showCalls) {
          detailsHtml += `<span><span class="legend-indicator" style="background-color: ${s.callsColor};"></span> ${t('stats.metricTypeCalls')}: <strong>${formatTelemetryNumber(cVal)}</strong></span>`;
        }
        if (showPoints) {
          detailsHtml += `<span><span class="legend-indicator" style="background-color: ${s.pointsColor};"></span> ${t('stats.metricTypePoints')}: <strong>${formatTelemetryNumber(pVal)}</strong></span>`;
        }

        sourceRowsHtml += `
          <div class="tooltip-source-block">
            <div class="tooltip-source-title">
              <span class="tooltip-dot" style="background-color: ${showCalls ? s.callsColor : s.pointsColor};"></span>
              <strong>${s.label}</strong>
            </div>
            <div class="tooltip-source-details">
              ${detailsHtml}
            </div>
          </div>
        `;
      }

      // Total row for selected sources
      let totalSelectedCalls = 0;
      let totalSelectedPoints = 0;
      for (const s of activeConfigs) {
        totalSelectedCalls += item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
        totalSelectedPoints += item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;
      }

      let totalDetailsHtml = '';
      if (showCalls) {
        totalDetailsHtml += `<span><span class="legend-indicator dashed" style="background-color: ${totalColors.callsColor};"></span> ${t('stats.metricTypeCalls')}: <strong>${formatTelemetryNumber(totalSelectedCalls)}</strong></span>`;
      }
      if (showPoints) {
        totalDetailsHtml += `<span><span class="legend-indicator dashed" style="background-color: ${totalColors.pointsColor};"></span> ${t('stats.metricTypePoints')}: <strong>${formatTelemetryNumber(totalSelectedPoints)}</strong></span>`;
      }

      const totalHtml = `
        <div class="tooltip-row total-row">
          <div class="tooltip-source-title">
            <span class="tooltip-dot dashed" style="background-color: ${showCalls ? totalColors.callsColor : totalColors.pointsColor};"></span>
            <strong>${t('stats.seriesTotal')}</strong>
          </div>
          <div class="tooltip-source-details">
            ${totalDetailsHtml}
          </div>
        </div>
      `;

      tooltip.innerHTML = `
        <div class="tooltip-header">${item.date}</div>
        ${sourceRowsHtml}
        ${totalHtml}
      `;

      tooltip.style.display = 'block';
      const tooltipRect = tooltip.getBoundingClientRect();
      const relX = clientX - rect.left;

      let leftPx = relX + 15;
      if (leftPx + tooltipRect.width > rect.width) {
        leftPx = relX - tooltipRect.width - 15;
      }
      leftPx = Math.max(10, Math.min(leftPx, rect.width - tooltipRect.width - 10));
      tooltip.style.left = `${leftPx}px`;

      // Ergonomic vertical positioning: position ABOVE finger/pointer by default
      // so finger doesn't obscure tooltip content; flip below if near top of chart
      let topPx = clientY - rect.top - tooltipRect.height - 18;
      if (topPx < 10) {
        topPx = clientY - rect.top + 24;
      }
      topPx = Math.max(5, Math.min(topPx, rect.height - tooltipRect.height - 5));
      tooltip.style.top = `${topPx}px`;
    };

    svg.addEventListener('mousemove', (e: MouseEvent) => {
      updateCrosshairAndTooltip(e.clientX, e.clientY);
    });

    svg.addEventListener('mouseleave', () => {
      hideCrosshairAndTooltip();
    });

    svg.addEventListener(
      'touchstart',
      (e: TouchEvent) => {
        if (e.touches.length > 0) {
          updateCrosshairAndTooltip(e.touches[0].clientX, e.touches[0].clientY);
        }
      },
      { passive: true },
    );

    svg.addEventListener(
      'touchmove',
      (e: TouchEvent) => {
        if (e.touches.length > 0) {
          updateCrosshairAndTooltip(e.touches[0].clientX, e.touches[0].clientY);
        }
      },
      { passive: true },
    );

    svg.addEventListener('touchend', () => {
      hideCrosshairAndTooltip();
    });

    svg.addEventListener('touchcancel', () => {
      hideCrosshairAndTooltip();
    });
  }

  private formatDateLabel(dateStr: string): string {
    if (dateStr.length === 7) return dateStr;
    const parts = dateStr.split('-');
    if (parts.length === 3) return `${parts[1]}/${parts[2]}`;
    return dateStr;
  }

  private formatNumber(num: number): string {
    return formatTelemetryCompact(num, 1);
  }

  private formatCompact(num: number): string {
    return formatTelemetryCompact(num, 0);
  }
}
