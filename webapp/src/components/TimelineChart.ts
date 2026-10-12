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
    calls: '#6366f1',    // Tech Indigo
    points: '#a5b4fc',   // Light Indigo / Lilac
    labelKey: 'stats.sourceRestApi',
  },
  web_map: {
    calls: '#10b981',    // Emerald / Green
    points: '#6ee7b7',   // Mint / Light Emerald
    labelKey: 'stats.sourceWebMap',
  },
  desktop_cli: {
    calls: '#f59e0b',    // Amber / Orange
    points: '#fde68a',   // Warm Yellow
    labelKey: 'stats.sourceDesktopCli',
  },
  unknown: {
    calls: '#94a3b8',    // Silver Slate
    points: '#cbd5e1',   // Light Slate
    labelKey: 'stats.sourceUnknown',
  },
};

// Distinct fallback palette pairs for dynamic sources beyond the known ones
const FALLBACK_PALETTE_PAIRS = [
  { calls: '#ec4899', points: '#fbcfe8' }, // Pink
  { calls: '#8b5cf6', points: '#ddd6fe' }, // Purple
  { calls: '#14b8a6', points: '#99f6e4' }, // Teal
  { calls: '#f97316', points: '#fed7aa' }, // Bright Orange
  { calls: '#06b6d4', points: '#a5f3fc' }, // Cyan
];

export const TOTAL_SERIES_CONFIG = {
  key: 'total',
  labelKey: 'stats.seriesTotal',
  callsColor: '#f8fafc',  // Crisp White dashed line for Total Invocations
  pointsColor: '#38bdf8', // Vibrant Cyan dashed line for Total Coordinates
};

export interface TimelineChartOptions {
  activeSourceKeys: Set<string>;
  availableSources: SourceItem[];
  showTotal: boolean;
}

export class TimelineChart {
  private container: HTMLElement;
  private data: TimelineSeriesItem[] = [];
  private options: TimelineChartOptions = {
    activeSourceKeys: new Set(['rest_api', 'web_map', 'desktop_cli']),
    availableSources: [],
    showTotal: true,
  };

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'timeline-chart-wrapper';
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
    }
    this.render();
  }

  public getSourceConfig(sourceName: string, index: number = 0): SourceDisplayConfig {
    if (KNOWN_SOURCE_COLORS[sourceName]) {
      const known = KNOWN_SOURCE_COLORS[sourceName];
      return {
        key: sourceName,
        label: t(known.labelKey),
        callsColor: known.calls,
        pointsColor: known.points,
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
      callsColor: pair.calls,
      pointsColor: pair.points,
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

    // Chart dimensions
    const width = 880;
    const height = 350;
    const padding = { top: 25, right: 75, bottom: 45, left: 65 };
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

    // Compute maximums for dual Y-axes
    // Left axis = Invocations / Calls
    // Right axis = Coordinates processed / Points
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

      // Total line for selected sources
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
      leftAxisLabelsSvg += `<text x="${padding.left - 8}" y="${y + 4}" text-anchor="end" class="chart-axis-text axis-calls">${this.formatNumber(callTickVal)}</text>`;
      rightAxisLabelsSvg += `<text x="${width - padding.right + 8}" y="${y + 4}" text-anchor="start" class="chart-axis-text axis-points">${this.formatCompact(pointTickVal)}</text>`;
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

      // Calls line (Solid line in Calls Color)
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

      // Coordinates line (Solid line in Points Color)
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

    // 2. Persistent Total Lines (calculated from selected sources)
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

    const totalCallPathD = this.buildSmoothPath(totalCallCoords);
    const totalPointPathD = this.buildSmoothPath(totalPointCoords);

    // Total Invocations line (Dashed white/light)
    pathsSvg += `
      <path 
        d="${totalCallPathD}" 
        fill="none" 
        stroke="${TOTAL_SERIES_CONFIG.callsColor}" 
        stroke-width="2.2" 
        stroke-dasharray="4 3"
        stroke-linecap="round" 
        stroke-linejoin="round"
        class="chart-total-line total-calls"
        data-series="total-calls"
      />
    `;

    // Total Coordinates processed line (Dashed cyan)
    pathsSvg += `
      <path 
        d="${totalPointPathD}" 
        fill="none" 
        stroke="${TOTAL_SERIES_CONFIG.pointsColor}" 
        stroke-width="2.2" 
        stroke-dasharray="4 3"
        stroke-linecap="round" 
        stroke-linejoin="round"
        class="chart-total-line total-points"
        data-series="total-points"
      />
    `;

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

      // Tooltip rows: Sources breakdown
      let sourceRowsHtml = '';
      for (const s of activeConfigs) {
        const cVal = item.calls.sources?.[s.key] ?? (item.calls[s.key] as number) ?? 0;
        const pVal = item.points.sources?.[s.key] ?? (item.points[s.key] as number) ?? 0;
        sourceRowsHtml += `
          <div class="tooltip-source-block">
            <div class="tooltip-source-title">
              <span class="tooltip-dot" style="background-color: ${s.callsColor};"></span>
              <strong>${s.label}</strong>
            </div>
            <div class="tooltip-source-details">
              <span><span class="legend-indicator" style="background-color: ${s.callsColor};"></span> ${t('stats.metricTypeCalls')}: <strong>${cVal.toLocaleString()}</strong></span>
              <span><span class="legend-indicator" style="background-color: ${s.pointsColor};"></span> ${t('stats.metricTypePoints')}: <strong>${pVal.toLocaleString()}</strong></span>
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

      const totalHtml = `
        <div class="tooltip-row total-row">
          <div class="tooltip-source-title">
            <span class="tooltip-dot dashed" style="background-color: ${TOTAL_SERIES_CONFIG.callsColor};"></span>
            <strong>${t('stats.seriesTotal')}</strong>
          </div>
          <div class="tooltip-source-details">
            <span><span class="legend-indicator dashed" style="background-color: ${TOTAL_SERIES_CONFIG.callsColor};"></span> ${t('stats.metricTypeCalls')}: <strong>${totalSelectedCalls.toLocaleString()}</strong></span>
            <span><span class="legend-indicator dashed" style="background-color: ${TOTAL_SERIES_CONFIG.pointsColor};"></span> ${t('stats.metricTypePoints')}: <strong>${totalSelectedPoints.toLocaleString()}</strong></span>
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

  private formatDateLabel(dateStr: string): string {
    if (dateStr.length === 7) return dateStr;
    const parts = dateStr.split('-');
    if (parts.length === 3) return `${parts[1]}/${parts[2]}`;
    return dateStr;
  }

  private formatNumber(num: number): string {
    if (num >= 1_000_000) return `${(num / 1_000_000).toFixed(1)}M`;
    if (num >= 1_000) return `${(num / 1_000).toFixed(1)}k`;
    return String(num);
  }

  private formatCompact(num: number): string {
    if (num >= 1_000_000) return `${(num / 1_000_000).toFixed(1)}M`;
    if (num >= 1_000) return `${(num / 1_000).toFixed(0)}k`;
    return String(num);
  }
}
