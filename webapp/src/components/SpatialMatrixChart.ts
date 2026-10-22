import { formatTelemetryNumber } from '../utils/numberFormat';
import { t, getLanguage } from '../i18n';
import type { SpatialMatrixResponse } from '../utils/apiClient';

export interface SpatialMatrixChartOptions {
  metric: 'points' | 'calls';
  scaleMode: 'log' | 'linear';
}

// Fixed Romanian bounding box for macroscopic punchcard grid
export const MATRIX_BOUNDS = {
  minLat: 43.5,
  maxLat: 48.5,
  minLon: 20.0,
  maxLon: 30.0,
  step: 0.1,
};

// Center point for mathematical quadrants
export const QUADRANT_CENTER = {
  lat: 45.9,
  lon: 24.9,
};

export class SpatialMatrixChart {
  private container: HTMLElement;
  private canvas: HTMLCanvasElement;
  private ctx: CanvasRenderingContext2D | null = null;
  private tooltipEl: HTMLElement;
  private data: SpatialMatrixResponse | null = null;
  private options: SpatialMatrixChartOptions = {
    metric: 'points',
    scaleMode: 'log',
  };

  // Lookup map: "lat:lon" -> hits
  private cellMap: Map<string, number> = new Map();
  private maxHits: number = 0;

  // Grid dimensions
  private cols: number = 0;
  private rows: number = 0;
  private cellSize: number = 8;
  private cellGap: number = 2;
  private paddingLeft: number = 42;
  private paddingRight: number = 16;
  private paddingTop: number = 16;
  private paddingBottom: number = 28;

  // Hover state
  private hoveredCell: { lat: number; lon: number; hits: number } | null = null;

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'spatial-matrix-chart-wrapper';

    this.canvas = document.createElement('canvas');
    this.canvas.className = 'spatial-matrix-canvas';
    this.ctx = this.canvas.getContext('2d');

    this.tooltipEl = document.createElement('div');
    this.tooltipEl.className = 'spatial-cell-tooltip';
    this.tooltipEl.style.display = 'none';

    this.container.appendChild(this.canvas);
    this.container.appendChild(this.tooltipEl);

    this.cols = Math.round((MATRIX_BOUNDS.maxLon - MATRIX_BOUNDS.minLon) / MATRIX_BOUNDS.step) + 1;
    this.rows = Math.round((MATRIX_BOUNDS.maxLat - MATRIX_BOUNDS.minLat) / MATRIX_BOUNDS.step) + 1;

    this.initEvents();
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public update(data: SpatialMatrixResponse, options?: Partial<SpatialMatrixChartOptions>): void {
    this.data = data;
    if (options) {
      if (options.metric) this.options.metric = options.metric;
      if (options.scaleMode) this.options.scaleMode = options.scaleMode;
    }

    this.cellMap.clear();
    this.maxHits = data.summary.max_cell_hits || 0;

    for (const cell of data.cells) {
      const key = `${cell.lat.toFixed(1)}:${cell.lon.toFixed(1)}`;
      this.cellMap.set(key, cell.hits);
    }

    this.render();
  }

  private initEvents(): void {
    this.canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => this.handleMouseLeave());
    window.addEventListener('resize', () => {
      if (this.data) this.render();
    });

    if (typeof MutationObserver !== 'undefined') {
      const observer = new MutationObserver(() => {
        if (this.data) this.render();
      });
      observer.observe(document.documentElement, {
        attributes: true,
        attributeFilter: ['data-theme', 'data-style'],
      });
    }
  }

  private isLightTheme(): boolean {
    return document.documentElement.getAttribute('data-theme') === 'light';
  }

  private render(): void {
    if (!this.ctx) return;

    const isLight = this.isLightTheme();
    const ghostFill = isLight ? 'rgba(0, 0, 0, 0.05)' : 'rgba(255, 255, 255, 0.04)';
    const ghostHoverStroke = isLight ? 'rgba(0, 0, 0, 0.35)' : 'rgba(255, 255, 255, 0.3)';
    const populatedHoverStroke = isLight ? '#18181b' : '#ffffff';

    // Responsive sizing based on container width
    const containerWidth = Math.max(this.container.clientWidth || 800, 320);
    // 101 cols -> calculate ideal cell size with gaps
    const availableWidth = containerWidth - this.paddingLeft - this.paddingRight;
    this.cellSize = Math.max(Math.floor((availableWidth - (this.cols * this.cellGap)) / this.cols), 6);
    this.cellGap = 2;

    const totalWidth = (this.cols * (this.cellSize + this.cellGap)) + this.paddingLeft + this.paddingRight;
    const totalHeight = (this.rows * (this.cellSize + this.cellGap)) + this.paddingTop + this.paddingBottom;

    // Handle HiDPI displays
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = totalWidth * dpr;
    this.canvas.height = totalHeight * dpr;
    this.canvas.style.width = `${totalWidth}px`;
    this.canvas.style.height = `${totalHeight}px`;

    this.ctx.resetTransform();
    this.ctx.scale(dpr, dpr);

    // Clear canvas
    this.ctx.clearRect(0, 0, totalWidth, totalHeight);

    // Draw background
    this.ctx.fillStyle = 'transparent';
    this.ctx.fillRect(0, 0, totalWidth, totalHeight);

    // Iterate rows (Lat: maxLat down to minLat so North is on Top)
    for (let r = 0; r < this.rows; r++) {
      const lat = Number((MATRIX_BOUNDS.maxLat - (r * MATRIX_BOUNDS.step)).toFixed(1));
      const y = this.paddingTop + (r * (this.cellSize + this.cellGap));

      for (let c = 0; c < this.cols; c++) {
        const lon = Number((MATRIX_BOUNDS.minLon + (c * MATRIX_BOUNDS.step)).toFixed(1));
        const x = this.paddingLeft + (c * (this.cellSize + this.cellGap));

        const key = `${lat.toFixed(1)}:${lon.toFixed(1)}`;
        const hits = this.cellMap.get(key) || 0;

        // Is hovered
        const isHovered = this.hoveredCell &&
          Math.abs(this.hoveredCell.lat - lat) < 0.05 &&
          Math.abs(this.hoveredCell.lon - lon) < 0.05;

        if (hits > 0) {
          // Populated cell
          this.ctx.fillStyle = this.getColorForHits(hits);
          this.ctx.fillRect(x, y, this.cellSize, this.cellSize);

          if (isHovered) {
            this.ctx.strokeStyle = populatedHoverStroke;
            this.ctx.lineWidth = 1.5;
            this.ctx.strokeRect(x - 0.5, y - 0.5, this.cellSize + 1, this.cellSize + 1);
          }
        } else {
          // Ghost cell (unvisited)
          this.ctx.fillStyle = ghostFill;
          this.ctx.fillRect(x, y, this.cellSize, this.cellSize);

          if (isHovered) {
            this.ctx.strokeStyle = ghostHoverStroke;
            this.ctx.lineWidth = 1;
            this.ctx.strokeRect(x - 0.5, y - 0.5, this.cellSize + 1, this.cellSize + 1);
          }
        }
      }
    }

    // Draw Subtle Mathematical Quadrant Divider Lines (Lat: 45.9°N, Lon: 24.9°E)
    this.renderQuadrantDividers();

    // Draw Subtle Axis Coordinate Labels (Every 1.0 degree)
    this.renderAxisLabels();
  }

  private renderQuadrantDividers(): void {
    if (!this.ctx) return;

    // Calculate Y position for Latitude 45.9°N
    const latOffset = (MATRIX_BOUNDS.maxLat - QUADRANT_CENTER.lat) / MATRIX_BOUNDS.step;
    const splitY = this.paddingTop + (latOffset * (this.cellSize + this.cellGap));

    // Calculate X position for Longitude 24.9°E
    const lonOffset = (QUADRANT_CENTER.lon - MATRIX_BOUNDS.minLon) / MATRIX_BOUNDS.step;
    const splitX = this.paddingLeft + (lonOffset * (this.cellSize + this.cellGap));

    const gridRight = this.paddingLeft + (this.cols * (this.cellSize + this.cellGap)) - this.cellGap;
    const gridBottom = this.paddingTop + (this.rows * (this.cellSize + this.cellGap)) - this.cellGap;

    this.ctx.save();
    this.ctx.strokeStyle = this.isLightTheme()
      ? 'rgba(2, 132, 199, 0.5)'
      : 'rgba(56, 189, 248, 0.45)';
    this.ctx.lineWidth = 1.5;
    this.ctx.setLineDash([4, 4]);

    // Horizontal line (45.9°N)
    this.ctx.beginPath();
    this.ctx.moveTo(this.paddingLeft, splitY);
    this.ctx.lineTo(gridRight, splitY);
    this.ctx.stroke();

    // Vertical line (24.9°E)
    this.ctx.beginPath();
    this.ctx.moveTo(splitX, this.paddingTop);
    this.ctx.lineTo(splitX, gridBottom);
    this.ctx.stroke();

    this.ctx.restore();
  }

  private renderAxisLabels(): void {
    if (!this.ctx) return;

    this.ctx.fillStyle = this.isLightTheme()
      ? 'rgba(87, 96, 106, 0.85)'
      : 'rgba(139, 148, 158, 0.8)';
    this.ctx.font = '10px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, monospace';

    // Longitude labels at bottom
    this.ctx.textAlign = 'center';
    this.ctx.textBaseline = 'top';
    for (let c = 0; c < this.cols; c++) {
      const lon = Number((MATRIX_BOUNDS.minLon + (c * MATRIX_BOUNDS.step)).toFixed(1));
      if (Math.abs(lon - Math.round(lon)) < 0.01) {
        const x = this.paddingLeft + (c * (this.cellSize + this.cellGap)) + (this.cellSize / 2);
        const y = this.paddingTop + (this.rows * (this.cellSize + this.cellGap)) + 4;
        this.ctx.fillText(`${lon.toFixed(0)}°E`, x, y);
      }
    }

    // Latitude labels at LEFT (standard cartographic placement, no clipping)
    this.ctx.textAlign = 'right';
    this.ctx.textBaseline = 'middle';
    for (let r = 0; r < this.rows; r++) {
      const lat = Number((MATRIX_BOUNDS.maxLat - (r * MATRIX_BOUNDS.step)).toFixed(1));
      if (Math.abs(lat - Math.round(lat)) < 0.01) {
        const x = this.paddingLeft - 6;
        const y = this.paddingTop + (r * (this.cellSize + this.cellGap)) + (this.cellSize / 2);
        this.ctx.fillText(`${lat.toFixed(0)}°N`, x, y);
      }
    }
  }

  private getColorForHits(hits: number): string {
    const isLight = this.isLightTheme();
    if (hits <= 0) return isLight ? 'rgba(0, 0, 0, 0.05)' : 'rgba(255, 255, 255, 0.04)';

    // Compute normalized ratio (0.0 to 1.0)
    let ratio: number;
    if (this.options.scaleMode === 'log') {
      const logMax = Math.log10(Math.max(this.maxHits, 10));
      const logVal = Math.log10(Math.max(hits, 1));
      ratio = Math.min(Math.max(logVal / logMax, 0.05), 1.0);
    } else {
      ratio = Math.min(Math.max(hits / Math.max(this.maxHits, 1), 0.05), 1.0);
    }

    if (isLight) {
      if (ratio < 0.25) return '#bae6fd';
      if (ratio < 0.50) return '#38bdf8';
      if (ratio < 0.75) return '#0284c7';
      return '#0369a1';
    }

    // 5-Tier Discrete Color Palette (Obsidian Punchcard / Cyan Glow Ramp)
    if (ratio < 0.25) {
      return '#1e3a5f';
    } else if (ratio < 0.50) {
      return '#0284c7';
    } else if (ratio < 0.75) {
      return '#38bdf8';
    } else {
      return '#a5f3fc';
    }
  }

  private handleMouseMove(e: MouseEvent): void {
    const rect = this.canvas.getBoundingClientRect();
    const x = e.clientX - rect.left - this.paddingLeft;
    const y = e.clientY - rect.top - this.paddingTop;

    const colPitch = this.cellSize + this.cellGap;
    const rowPitch = this.cellSize + this.cellGap;

    const c = Math.floor(x / colPitch);
    const r = Math.floor(y / rowPitch);

    if (c >= 0 && c < this.cols && r >= 0 && r < this.rows) {
      const lat = Number((MATRIX_BOUNDS.maxLat - (r * MATRIX_BOUNDS.step)).toFixed(1));
      const lon = Number((MATRIX_BOUNDS.minLon + (c * MATRIX_BOUNDS.step)).toFixed(1));
      const key = `${lat.toFixed(1)}:${lon.toFixed(1)}`;
      const hits = this.cellMap.get(key) || 0;

      this.hoveredCell = { lat, lon, hits };
      this.showTooltip(e.clientX, e.clientY, lat, lon, hits);
      this.render();
    } else {
      this.handleMouseLeave();
    }
  }

  private handleMouseLeave(): void {
    if (this.hoveredCell) {
      this.hoveredCell = null;
      this.tooltipEl.style.display = 'none';
      this.render();
    }
  }

  private showTooltip(clientX: number, clientY: number, lat: number, lon: number, hits: number): void {
    // Quadrant determination
    const quad = lat >= QUADRANT_CENTER.lat
      ? (lon < QUADRANT_CENTER.lon ? 'NW' : 'NE')
      : (lon < QUADRANT_CENTER.lon ? 'SW' : 'SE');

    const isRo = getLanguage() === 'ro';
    const fmt = (v: number) => isRo ? v.toFixed(1).replace('.', ',') : v.toFixed(1);

    const latMin = fmt(lat);
    const latMax = fmt(lat + 0.1);
    const lonMin = fmt(lon);
    const lonMax = fmt(lon + 0.1);

    const metricLabel = this.options.metric === 'points'
      ? t('stats.metricTypePoints')
      : t('stats.metricTypeCalls');

    this.tooltipEl.innerHTML = `
      <div class="tooltip-header">
        <span class="tooltip-quadrant-badge">${quad}</span>
        <span class="tooltip-coords-center">${fmt(lat)}°N, ${fmt(lon)}°E</span>
      </div>
      <div class="tooltip-body">
        <div class="tooltip-row">
          <span class="muted">${t('stats.cellTooltipBounds')}:</span>
          <span class="mono">${latMin}°–${latMax}°N, ${lonMin}°–${lonMax}°E</span>
        </div>
        <div class="tooltip-row">
          <span class="muted">${metricLabel}:</span>
          <span class="mono font-bold primary-light">${hits > 0 ? formatTelemetryNumber(hits) : '0'}</span>
        </div>
      </div>
    `;

    this.tooltipEl.style.display = 'block';

    const rect = this.container.getBoundingClientRect();
    const tooltipX = clientX - rect.left + 14;
    const tooltipY = clientY - rect.top + 14;

    this.tooltipEl.style.left = `${tooltipX}px`;
    this.tooltipEl.style.top = `${tooltipY}px`;
  }
}
