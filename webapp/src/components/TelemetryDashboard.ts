import noUiSlider, { type API as NoUiSliderAPI } from 'nouislider';
import 'nouislider/dist/nouislider.css';
import { t, onLanguageChange } from '../i18n';
import {
  apiClient,
  type GridInfoResponse,
  type TelemetryStatsResponse,
  type TimelineResponse,
  type SourceItem,
  type SpatialMatrixResponse,
} from '../utils/apiClient';
import {
  TimelineChart,
  TOTAL_SERIES_CONFIG,
} from './TimelineChart';
import { SpatialMatrixChart } from './SpatialMatrixChart';

export type TelemetryPreset = 'all' | '12m' | 'ytd' | '30d' | '7d';
export type Granularity = 'daily' | 'weekly' | 'monthly';

export interface TelemetryDashboardCallbacks {
  onBackToStudio: () => void;
}

export class TelemetryDashboard {
  private element: HTMLElement;
  private callbacks: TelemetryDashboardCallbacks;

  private gridInfo: GridInfoResponse | null = null;
  private telemetry: TelemetryStatsResponse | null = null;
  private timelineData: TimelineResponse | null = null;
  private spatialData: SpatialMatrixResponse | null = null;

  // View state
  private activeSubTab: 'timeline' | 'spatial' | 'grid' = 'timeline';
  private activePreset: TelemetryPreset = 'all';
  private activeGranularity: Granularity = 'daily';

  // Spatial view state
  private spatialMetric: 'points' | 'calls' = 'points';
  private spatialScaleMode: 'log' | 'linear' = 'log';
  private spatialPreset: TelemetryPreset = 'all';
  private spatialActiveSourceKeys: Set<string> = new Set();
  private isSpatialLoading: boolean = false;
  private spatialHasCustomRange: boolean = false;
  private spatialCustomStartDate: string = '2023-11-01';
  private spatialCustomEndDate: string = new Date().toISOString().slice(0, 10);

  // Dynamic sources list and selected toggle set
  private availableSources: SourceItem[] = [];
  private activeSourceKeys: Set<string> = new Set();

  private isLoading: boolean = false;

  // Absolute bounds of available usage data from DB
  private minHistoryDate: string = '2023-11-01';
  private maxHistoryDate: string = new Date().toISOString().slice(0, 10);

  // Current selected range
  private hasCustomRange: boolean = false;
  private customStartDate: string = '2023-11-01';
  private customEndDate: string = new Date().toISOString().slice(0, 10);

  private sliderInstance: NoUiSliderAPI | null = null;

  private chart: TimelineChart;
  private spatialChart: SpatialMatrixChart;

  constructor(callbacks: TelemetryDashboardCallbacks) {
    this.callbacks = callbacks;
    this.element = document.createElement('div');
    this.element.className = 'telemetry-page-view';
    this.chart = new TimelineChart();
    this.spatialChart = new SpatialMatrixChart();

    this.render();

    onLanguageChange(() => {
      this.render();
      if (this.activeSubTab === 'timeline') this.mountChart();
      if (this.activeSubTab === 'spatial') this.mountSpatialChart();
    });

    this.initData();
  }

  public getElement(): HTMLElement {
    return this.element;
  }

  public async refresh(): Promise<void> {
    await this.initData();
  }

  private async initData(): Promise<void> {
    this.isLoading = true;
    try {
      const [gridRes, teleRes] = await Promise.all([
        apiClient.getGridInfo(),
        apiClient.getTelemetryStats(),
      ]);
      this.gridInfo = gridRes;
      this.telemetry = teleRes;

      await this.fetchTimeline();
    } catch (err) {
      console.warn('Failed to initialize telemetry dashboard:', err);
    } finally {
      this.isLoading = false;
      this.render();
      this.mountChart();
    }
  }

  private async fetchTimeline(): Promise<void> {
    this.isLoading = true;
    try {
      const params: Parameters<typeof apiClient.getTelemetryTimeline>[0] = {
        granularity: this.activeGranularity,
      };
      if (this.activePreset !== 'all' || this.hasCustomRange) {
        params.from_date = this.customStartDate;
        params.to_date = this.customEndDate;
      }

      const res = await apiClient.getTelemetryTimeline(params);

      this.timelineData = res;

      // Extract available sources dynamically from response
      if (res.available_sources && res.available_sources.length > 0) {
        this.availableSources = res.available_sources;
      } else {
        // Fallback discovery from series items if available_sources empty
        const discovered = new Set<string>();
        for (const item of res.series) {
          if (item.calls.sources) {
            Object.keys(item.calls.sources).forEach((k) => discovered.add(k));
          }
        }
        if (discovered.size > 0) {
          this.availableSources = Array.from(discovered).map((name, i) => ({
            id: i,
            name,
            description: null,
          }));
        } else if (this.availableSources.length === 0) {
          this.availableSources = [
            { id: 1, name: 'rest_api', description: 'REST API' },
            { id: 2, name: 'web_map', description: 'Web Map' },
            { id: 3, name: 'desktop_cli', description: 'Desktop / CLI' },
          ];
        }
      }

      // Initialize active toggles on first load: activate all non-zero/known sources
      if (this.activeSourceKeys.size === 0) {
        this.availableSources.forEach((s) => {
          if (s.name !== 'unknown') {
            this.activeSourceKeys.add(s.name);
          }
        });
        if (this.activeSourceKeys.size === 0 && this.availableSources.length > 0) {
          this.activeSourceKeys.add(this.availableSources[0].name);
        }
      }

      // Track earliest available date from database
      if (res.summary.from_date && this.activePreset === 'all') {
        this.minHistoryDate = res.summary.from_date;
        if (!this.hasCustomRange) {
          this.customStartDate = this.minHistoryDate;
        }
      }
    } catch (err) {
      console.warn('Failed to fetch timeline series:', err);
    } finally {
      this.isLoading = false;
      this.render();
      this.mountChart();
    }
  }

  private applyPreset(preset: TelemetryPreset): void {
    this.activePreset = preset;
    this.hasCustomRange = false;
    const now = new Date();
    this.customEndDate = now.toISOString().slice(0, 10);

    if (preset === '7d') {
      const d = new Date();
      d.setDate(d.getDate() - 7);
      this.customStartDate = d.toISOString().slice(0, 10);
    } else if (preset === '30d') {
      const d = new Date();
      d.setDate(d.getDate() - 30);
      this.customStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'ytd') {
      this.customStartDate = `${now.getFullYear()}-01-01`;
    } else if (preset === '12m') {
      const d = new Date();
      d.setFullYear(d.getFullYear() - 1);
      this.customStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'all') {
      this.customStartDate = this.minHistoryDate;
    }

    this.fetchTimeline();
  }

  private toggleSource(sourceKey: string): void {
    if (this.activeSourceKeys.has(sourceKey)) {
      this.activeSourceKeys.delete(sourceKey);
    } else {
      this.activeSourceKeys.add(sourceKey);
    }
    this.render();
    this.mountChart();
  }

  private render(): void {
    this.element.innerHTML = `
      <div class="telemetry-page-container">
        <!-- Top Navigation Bar -->
        <div class="telemetry-page-header">
          <div class="telemetry-title-block">
            <button class="btn btn-outline" id="btn-back-to-studio">
              ${t('stats.backToStudio')}
            </button>
            <div class="telemetry-page-title">${t('stats.title')}</div>
          </div>

          <div class="telemetry-subtabs">
            <button class="telemetry-subtab-btn ${this.activeSubTab === 'timeline' ? 'active' : ''}" id="tab-timeline-btn">
              ${t('stats.tabTimeline')}
            </button>
            <button class="telemetry-subtab-btn ${this.activeSubTab === 'spatial' ? 'active' : ''}" id="tab-spatial-btn">
              ${t('stats.tabSpatial')}
            </button>
            <button class="telemetry-subtab-btn ${this.activeSubTab === 'grid' ? 'active' : ''}" id="tab-grid-btn">
              ${t('stats.tabGridInfo')}
            </button>
          </div>
        </div>

        <!-- Main Body -->
        <div class="telemetry-page-content">
          ${this.activeSubTab === 'timeline'
            ? this.renderTimelineHtml()
            : (this.activeSubTab === 'spatial' ? this.renderSpatialHtml() : this.renderGridHtml())}
        </div>
      </div>
    `;

    // Top Navigation Listeners
    this.element.querySelector('#btn-back-to-studio')?.addEventListener('click', () => {
      this.callbacks.onBackToStudio();
    });

    this.element.querySelector('#tab-timeline-btn')?.addEventListener('click', () => {
      this.activeSubTab = 'timeline';
      this.render();
      this.mountChart();
    });

    this.element.querySelector('#tab-spatial-btn')?.addEventListener('click', () => {
      this.activeSubTab = 'spatial';
      this.render();
      this.mountSpatialChart();
      if (!this.spatialData) {
        this.fetchSpatialMatrix();
      }
    });

    this.element.querySelector('#tab-grid-btn')?.addEventListener('click', () => {
      this.activeSubTab = 'grid';
      this.render();
    });

    if (this.activeSubTab === 'timeline') {
      this.attachTimelineEventListeners();
    } else if (this.activeSubTab === 'spatial') {
      this.attachSpatialEventListeners();
    }
  }

  private renderTimelineHtml(): string {
    // Dynamic calculation of total calls and total points from ACTIVE SELECTED sources
    let totalCalls = 0;
    let totalPoints = 0;

    if (this.timelineData?.series) {
      for (const item of this.timelineData.series) {
        for (const srcKey of this.activeSourceKeys) {
          totalCalls += item.calls.sources?.[srcKey] ?? (item.calls[srcKey] as number) ?? 0;
          totalPoints += item.points.sources?.[srcKey] ?? (item.points[srcKey] as number) ?? 0;
        }
      }
    } else {
      totalCalls = this.timelineData?.summary?.total_calls ?? (this.telemetry?.summary.total_calls ?? 0);
      totalPoints = this.timelineData?.summary?.total_points ?? (this.telemetry?.summary.total_points ?? 0);
    }

    const avgRatio = totalCalls > 0 ? (totalPoints / totalCalls).toFixed(1) : '0';

    let periodText = 'All Logged History';
    const summary = this.timelineData?.summary;
    const startStr = summary?.from_date ?? this.customStartDate;
    const endStr = summary?.to_date ?? this.customEndDate;

    if (startStr && endStr) {
      const d1 = new Date(startStr);
      const d2 = new Date(endStr);
      const daysDiff = Math.max(0, Math.round((d2.getTime() - d1.getTime()) / 86400000) + 1);

      let countLabel = '';
      if (this.activeGranularity === 'monthly') {
        const monthsCount = Math.max(1, (d2.getFullYear() - d1.getFullYear()) * 12 + (d2.getMonth() - d1.getMonth()) + 1);
        countLabel = t('stats.monthsCount', { count: monthsCount });
      } else if (this.activeGranularity === 'weekly') {
        const weeksCount = Math.max(1, Math.ceil(daysDiff / 7));
        countLabel = t('stats.weeksCount', { count: weeksCount });
      } else {
        countLabel = t('stats.daysCount', { count: daysDiff });
      }
      periodText = `${startStr} → ${endStr} (${countLabel})`;
    }

    // Dynamic Source Toggles with dual-color indicators (Count & Coordinates)
    const sourceTogglesHtml = this.availableSources
      .map((s, idx) => {
        const isActive = this.activeSourceKeys.has(s.name);
        const cfg = this.chart.getSourceConfig(s.name, idx);
        return `
          <label class="source-toggle-card ${isActive ? 'active' : ''}" data-source-toggle="${s.name}">
            <input type="checkbox" class="source-checkbox" ${isActive ? 'checked' : ''} data-source-name="${s.name}" />
            <div class="source-swatch-duo">
              <span class="swatch-bar swatch-calls" style="background-color: ${cfg.callsColor};" title="${t('stats.metricTypeCalls')}"></span>
              <span class="swatch-bar swatch-points" style="background-color: ${cfg.pointsColor};" title="${t('stats.metricTypePoints')}"></span>
            </div>
            <span class="source-name-label">${cfg.label}</span>
          </label>
        `;
      })
      .join('');

    // Persistent Total indicator legend card
    const totalLegendCardHtml = `
      <div class="source-toggle-card persistent-total" title="${t('stats.seriesTotal')}">
        <div class="source-swatch-duo">
          <span class="swatch-bar swatch-calls dashed" style="background-color: ${TOTAL_SERIES_CONFIG.callsColor};"></span>
          <span class="swatch-bar swatch-points dashed" style="background-color: ${TOTAL_SERIES_CONFIG.pointsColor};"></span>
        </div>
        <span class="source-name-label font-bold">${t('stats.seriesTotal')}</span>
      </div>
    `;

    return `
      <!-- KPI Banner -->
      <div class="telemetry-kpi-banner">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val blue">${totalCalls.toLocaleString()}</div>
          <div class="telemetry-kpi-label">${t('stats.totalCalls')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val cyan">${totalPoints.toLocaleString()}</div>
          <div class="telemetry-kpi-label">${t('stats.totalPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val emerald">${avgRatio}</div>
          <div class="telemetry-kpi-label">${t('stats.batchIntensity')}</div>
        </div>
        <div class="telemetry-kpi-card wide">
          <div class="telemetry-kpi-val muted" style="font-size: 13.5px;">${periodText}</div>
          <div class="telemetry-kpi-label">${t('stats.activePeriod')}</div>
        </div>
      </div>

      <!-- Controls Deck -->
      <div class="telemetry-deck">
        <!-- Row 1: Presets & Granularity -->
        <div class="telemetry-deck-row">
          <div class="telemetry-btn-toggle-group">
            <button class="telemetry-toggle-btn ${this.activePreset === '7d' ? 'active' : ''}" data-preset="7d">${t('stats.preset7D')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === '30d' ? 'active' : ''}" data-preset="30d">${t('stats.preset30D')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === 'ytd' ? 'active' : ''}" data-preset="ytd">${t('stats.presetYTD')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === '12m' ? 'active' : ''}" data-preset="12m">${t('stats.preset12M')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === 'all' ? 'active' : ''}" data-preset="all">${t('stats.presetAll')}</button>
          </div>

          <div class="telemetry-btn-toggle-group">
            <button class="telemetry-toggle-btn ${this.activeGranularity === 'daily' ? 'active' : ''}" data-granularity="daily">${t('stats.granularityDaily')}</button>
            <button class="telemetry-toggle-btn ${this.activeGranularity === 'weekly' ? 'active' : ''}" data-granularity="weekly">${t('stats.granularityWeekly')}</button>
            <button class="telemetry-toggle-btn ${this.activeGranularity === 'monthly' ? 'active' : ''}" data-granularity="monthly">${t('stats.granularityMonthly')}</button>
          </div>
        </div>

        <!-- Row 2: noUiSlider Graphical Time Window Slider -->
        <div class="telemetry-range-slider-block">
          <div class="range-slider-header">
            <span class="range-slider-title">${t('stats.rangeSliderLabel')}</span>
            <div class="range-slider-readout">
              <span class="range-date-pill" id="range-start-label">${this.customStartDate}</span>
              <span class="range-arrow">→</span>
              <span class="range-date-pill" id="range-end-label">${this.customEndDate}</span>
            </div>
          </div>

          <div class="nouislider-container">
            <div id="nouislider-time-window"></div>
          </div>

          <div class="range-slider-limits">
            <span class="slider-limit-label">${t('stats.sliderEarliest', { date: this.minHistoryDate })}</span>
            <span class="slider-limit-label">${t('stats.sliderToday', { date: this.maxHistoryDate })}</span>
          </div>
        </div>

        <!-- Row 3: Dynamic Sources Toggle Boxes -->
        <div class="telemetry-deck-row source-toggles-row">
          <div class="sources-label-col">
            <span class="facet-label">${t('stats.sourcesFilterLabel')}</span>
          </div>
          <div class="sources-toggle-list">
            ${sourceTogglesHtml}
            ${totalLegendCardHtml}
          </div>
        </div>
      </div>

      <!-- Chart Viewport -->
      <div class="telemetry-chart-viewport full-page">
        ${this.isLoading ? `<div class="timeline-loading-overlay">${t('stats.loadingTimeline')}</div>` : ''}
        <div class="chart-legend-row dual-metric-legend">
          <div class="legend-badge">
            <span class="legend-box" style="background-color: #6366f1;"></span>
            <strong>${t('stats.metricTypeCalls')}</strong>
            <span class="muted font-normal">(${t('stats.callsLegend') || 'Calls'})</span>
          </div>
          <div class="legend-badge">
            <span class="legend-box" style="background-color: #38bdf8;"></span>
            <strong>${t('stats.metricTypePoints')}</strong>
            <span class="muted font-normal">(${t('stats.pointsLegend') || 'Coordinates'})</span>
          </div>
        </div>
        <div id="timeline-chart-slot"></div>
      </div>
    `;
  }

  private renderGridHtml(): string {
    return `
      <div class="grid-metadata-card">
        <div class="card-title">${t('stats.activeGridTitle')}</div>
        <div class="grid-info-table">
          <div class="grid-info-row">
            <span class="muted">${t('stats.gridFile')}</span>
            <span class="mono primary-light font-bold">${this.gridInfo?.grid_file || 'transdat.spg'}</span>
          </div>
          <div class="grid-info-row">
            <span class="muted">${t('stats.quasigeoid')}</span>
            <span>${this.gridInfo?.quasigeoid_model || 'Model Cvasigeoid 2008'}</span>
          </div>
          <div class="grid-info-row">
            <span class="muted">${t('stats.projCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.projected || 'EPSG:31700 (Stereo 70)'}</span>
          </div>
          <div class="grid-info-row">
            <span class="muted">${t('stats.geoCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.geographic || 'EPSG:4258 (ETRS89)'}</span>
          </div>
        </div>

        <div class="accuracy-callout">
          <div class="accuracy-callout-title">${t('stats.accuracyNoteTitle')}</div>
          <div class="accuracy-callout-text">${t('stats.accuracyNoteText')}</div>
        </div>
      </div>
    `;
  }

  private attachTimelineEventListeners(): void {
    // Presets
    this.element.querySelectorAll('[data-preset]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const p = (e.currentTarget as HTMLElement).getAttribute('data-preset') as TelemetryPreset;
        if (p) this.applyPreset(p);
      });
    });

    // Granularity
    this.element.querySelectorAll('[data-granularity]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const g = (e.currentTarget as HTMLElement).getAttribute('data-granularity') as Granularity;
        if (g && g !== this.activeGranularity) {
          this.activeGranularity = g;
          this.fetchTimeline();
        }
      });
    });

    // noUiSlider Initialization
    this.initNoUiSlider();

    // Source Toggles Checkboxes
    this.element.querySelectorAll<HTMLInputElement>('.source-checkbox').forEach((checkbox) => {
      checkbox.addEventListener('change', (e) => {
        const target = e.target as HTMLInputElement;
        const srcName = target.getAttribute('data-source-name');
        if (srcName) {
          this.toggleSource(srcName);
        }
      });
    });
  }

  private initNoUiSlider(): void {
    const sliderContainer = this.element.querySelector('#nouislider-time-window') as HTMLElement | null;
    if (!sliderContainer) return;

    if (this.sliderInstance) {
      try {
        this.sliderInstance.destroy();
      } catch {
        // Ignore if already destroyed
      }
      this.sliderInstance = null;
    }

    const minTs = new Date(this.minHistoryDate).getTime();
    const maxTs = new Date(this.maxHistoryDate).getTime();
    const startTs = Math.max(minTs, Math.min(new Date(this.customStartDate).getTime(), maxTs));
    const endTs = Math.max(minTs, Math.min(new Date(this.customEndDate).getTime(), maxTs));

    const slider = noUiSlider.create(sliderContainer, {
      start: [startTs, endTs],
      connect: true,
      behaviour: 'drag',
      step: 86400000,
      range: {
        min: minTs,
        max: Math.max(maxTs, minTs + 86400000),
      },
    });

    this.sliderInstance = slider;

    const labelStart = this.element.querySelector('#range-start-label');
    const labelEnd = this.element.querySelector('#range-end-label');

    slider.on('update', (values) => {
      const sTs = Number(values[0]);
      const eTs = Number(values[1]);
      const sStr = new Date(sTs).toISOString().slice(0, 10);
      const eStr = new Date(eTs).toISOString().slice(0, 10);

      if (labelStart) labelStart.textContent = sStr;
      if (labelEnd) labelEnd.textContent = eStr;
    });

    slider.on('change', (values) => {
      const sTs = Number(values[0]);
      const eTs = Number(values[1]);
      this.hasCustomRange = true;
      this.customStartDate = new Date(sTs).toISOString().slice(0, 10);
      this.customEndDate = new Date(eTs).toISOString().slice(0, 10);
      this.fetchTimeline();
    });
  }

  private mountChart(): void {
    const slot = this.element.querySelector('#timeline-chart-slot');
    if (!slot) return;
    slot.innerHTML = '';
    slot.appendChild(this.chart.getElement());

    if (this.timelineData?.series) {
      this.chart.update(this.timelineData.series, {
        activeSourceKeys: this.activeSourceKeys,
        availableSources: this.availableSources,
        showTotal: true,
      });
    }
  }

  // ==========================================
  // Spatial Matrix Visual Logic
  // ==========================================

  private async fetchSpatialMatrix(): Promise<void> {
    this.isSpatialLoading = true;
    this.render();

    try {
      const params: Parameters<typeof apiClient.getSpatialMatrix>[0] = {
        metric: this.spatialMetric,
      };

      // Date range filter
      if (this.spatialHasCustomRange) {
        params.from_date = this.spatialCustomStartDate;
        params.to_date = this.spatialCustomEndDate;
      } else {
        const now = new Date();
        if (this.spatialPreset === '7d') {
          const d = new Date();
          d.setDate(d.getDate() - 7);
          params.from_date = d.toISOString().slice(0, 10);
        } else if (this.spatialPreset === '30d') {
          const d = new Date();
          d.setDate(d.getDate() - 30);
          params.from_date = d.toISOString().slice(0, 10);
        } else if (this.spatialPreset === 'ytd') {
          params.from_date = `${now.getFullYear()}-01-01`;
        } else if (this.spatialPreset === '12m') {
          const d = new Date();
          d.setFullYear(d.getFullYear() - 1);
          params.from_date = d.toISOString().slice(0, 10);
        }
      }

      // Map active source keys to source IDs
      if (this.spatialActiveSourceKeys.size > 0 && this.availableSources.length > 0) {
        const ids = this.availableSources
          .filter((s) => this.spatialActiveSourceKeys.has(s.name))
          .map((s) => s.id);
        if (ids.length > 0 && ids.length < this.availableSources.length) {
          params.sources = ids.join(',');
        }
      }

      const res = await apiClient.getSpatialMatrix(params);
      this.spatialData = res;
    } catch (err) {
      console.warn('Failed to fetch spatial matrix telemetry:', err);
    } finally {
      this.isSpatialLoading = false;
      this.render();
      this.mountSpatialChart();
    }
  }

  private applySpatialPreset(preset: TelemetryPreset): void {
    this.spatialPreset = preset;
    this.spatialHasCustomRange = false;
    const now = new Date();
    this.spatialCustomEndDate = now.toISOString().slice(0, 10);

    if (preset === '7d') {
      const d = new Date();
      d.setDate(d.getDate() - 7);
      this.spatialCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === '30d') {
      const d = new Date();
      d.setDate(d.getDate() - 30);
      this.spatialCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'ytd') {
      this.spatialCustomStartDate = `${now.getFullYear()}-01-01`;
    } else if (preset === '12m') {
      const d = new Date();
      d.setFullYear(d.getFullYear() - 1);
      this.spatialCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'all') {
      this.spatialCustomStartDate = this.minHistoryDate;
    }

    this.fetchSpatialMatrix();
  }

  private renderSpatialHtml(): string {
    const quadrants = this.spatialData?.summary.quadrants;
    const totalVal = this.spatialData?.summary.total_metric_value || 0;
    const activeCells = this.spatialData?.summary.active_cells_count || 0;
    const maxCellHits = this.spatialData?.summary.max_cell_hits || 0;

    // Default quadrant items if not yet loaded
    const qNW = quadrants?.NW || { value: 0, percentage: 0 };
    const qNE = quadrants?.NE || { value: 0, percentage: 0 };
    const qSW = quadrants?.SW || { value: 0, percentage: 0 };
    const qSE = quadrants?.SE || { value: 0, percentage: 0 };

    // Source toggle chips
    const sourceChipsHtml = this.availableSources
      .filter((s) => s.name !== 'unknown')
      .map((s) => {
        const isActive = this.spatialActiveSourceKeys.size === 0 || this.spatialActiveSourceKeys.has(s.name);
        return `
          <button class="spatial-source-chip ${isActive ? 'active' : ''}" data-spatial-source="${s.name}">
            <span class="chip-check">${isActive ? '✓' : ''}</span>
            <span>${t(`stats.source${s.name.charAt(0).toUpperCase() + s.name.slice(1)}`) || s.name}</span>
          </button>
        `;
      })
      .join('');

    return `
      <!-- Top Deck: 4 Mathematical Quadrant Cards -->
      <div class="spatial-quadrants-deck">
        <div class="quadrant-card nw">
          <div class="quadrant-header">
            <span class="quadrant-code">NW</span>
            <span class="quadrant-name">${t('stats.quadrantNW')}</span>
          </div>
          <div class="quadrant-value mono">${qNW.value.toLocaleString()}</div>
          <div class="quadrant-share">${t('stats.shareOfTotal', { pct: qNW.percentage })}</div>
        </div>

        <div class="quadrant-card ne">
          <div class="quadrant-header">
            <span class="quadrant-code">NE</span>
            <span class="quadrant-name">${t('stats.quadrantNE')}</span>
          </div>
          <div class="quadrant-value mono">${qNE.value.toLocaleString()}</div>
          <div class="quadrant-share">${t('stats.shareOfTotal', { pct: qNE.percentage })}</div>
        </div>

        <div class="quadrant-card sw">
          <div class="quadrant-header">
            <span class="quadrant-code">SW</span>
            <span class="quadrant-name">${t('stats.quadrantSW')}</span>
          </div>
          <div class="quadrant-value mono">${qSW.value.toLocaleString()}</div>
          <div class="quadrant-share">${t('stats.shareOfTotal', { pct: qSW.percentage })}</div>
        </div>

        <div class="quadrant-card se">
          <div class="quadrant-header">
            <span class="quadrant-code">SE</span>
            <span class="quadrant-name">${t('stats.quadrantSE')}</span>
          </div>
          <div class="quadrant-value mono">${qSE.value.toLocaleString()}</div>
          <div class="quadrant-share">${t('stats.shareOfTotal', { pct: qSE.percentage })}</div>
        </div>
      </div>

      <!-- Controls Deck -->
      <div class="telemetry-deck spatial-controls-deck">
        <!-- Row 1: Time Horizon Presets -->
        <div class="telemetry-deck-row">
          <div class="telemetry-facet-group">
            <span class="facet-label">${t('stats.activePeriod')}:</span>
            <div class="btn-group">
              <button class="btn btn-sm ${this.spatialPreset === 'all' && !this.spatialHasCustomRange ? 'btn-primary' : 'btn-outline'}" data-spatial-preset="all">
                ${t('stats.presetAll')}
              </button>
              <button class="btn btn-sm ${this.spatialPreset === '12m' && !this.spatialHasCustomRange ? 'btn-primary' : 'btn-outline'}" data-spatial-preset="12m">
                ${t('stats.preset12M')}
              </button>
              <button class="btn btn-sm ${this.spatialPreset === 'ytd' && !this.spatialHasCustomRange ? 'btn-primary' : 'btn-outline'}" data-spatial-preset="ytd">
                ${t('stats.presetYTD')}
              </button>
              <button class="btn btn-sm ${this.spatialPreset === '30d' && !this.spatialHasCustomRange ? 'btn-primary' : 'btn-outline'}" data-spatial-preset="30d">
                ${t('stats.preset30D')}
              </button>
              <button class="btn btn-sm ${this.spatialPreset === '7d' && !this.spatialHasCustomRange ? 'btn-primary' : 'btn-outline'}" data-spatial-preset="7d">
                ${t('stats.preset7D')}
              </button>
            </div>
          </div>

          <!-- Total Count & Cell Area Badges -->
          <div class="spatial-overview-badges">
            <span class="stat-pill">
              <span class="muted">${t('stats.activeCellsCount', { count: activeCells.toLocaleString() })}</span>
            </span>
            <span class="stat-pill cell-area-pill">
              <span class="muted">${t('stats.cellAreaNote')}</span>
            </span>
            <span class="stat-pill highlight">
              <strong>${totalVal.toLocaleString()}</strong> ${this.spatialMetric === 'points' ? t('stats.metricPoints') : t('stats.metricCalls')}
            </span>
          </div>
        </div>

        <!-- Row 2: Vanilla Date Range Scrubber -->
        <div class="telemetry-deck-row spatial-date-scrubber-row">
          <div class="spatial-date-scrubber-label">
            <span class="facet-label">${t('stats.rangeSliderLabel')}</span>
            <span class="mono primary-light font-bold" id="spatial-range-display">
              ${this.spatialCustomStartDate} — ${this.spatialCustomEndDate}
            </span>
          </div>
          <div class="spatial-date-inputs">
            <input type="date" class="form-control form-control-sm" id="spatial-date-start"
              min="${this.minHistoryDate}" max="${this.maxHistoryDate}" value="${this.spatialCustomStartDate}" />
            <span class="muted">➔</span>
            <input type="date" class="form-control form-control-sm" id="spatial-date-end"
              min="${this.minHistoryDate}" max="${this.maxHistoryDate}" value="${this.spatialCustomEndDate}" />
          </div>
        </div>

        <!-- Row 3: Metric Toggle & Scale Mode Toggle (Zero Emojis) -->
        <div class="telemetry-deck-row facet-row">
          <div class="telemetry-facet-group">
            <span class="facet-label">${t('stats.metricToggleLabel')}</span>
            <div class="btn-group">
              <button class="btn btn-sm ${this.spatialMetric === 'points' ? 'btn-primary' : 'btn-outline'}" id="spatial-metric-points-btn">
                ${t('stats.metricPoints')}
              </button>
              <button class="btn btn-sm ${this.spatialMetric === 'calls' ? 'btn-primary' : 'btn-outline'}" id="spatial-metric-calls-btn">
                ${t('stats.metricCalls')}
              </button>
            </div>
          </div>

          <div class="telemetry-facet-group">
            <span class="facet-label">${t('stats.scaleToggleLabel')}</span>
            <div class="btn-group">
              <button class="btn btn-sm ${this.spatialScaleMode === 'log' ? 'btn-primary' : 'btn-outline'}" id="spatial-scale-log-btn">
                ${t('stats.scaleLog')}
              </button>
              <button class="btn btn-sm ${this.spatialScaleMode === 'linear' ? 'btn-primary' : 'btn-outline'}" id="spatial-scale-linear-btn">
                ${t('stats.scaleLinear')}
              </button>
            </div>
          </div>
        </div>

        <!-- Row 4: Client Sources Filter Chips -->
        <div class="telemetry-deck-row source-toggles-row">
          <div class="sources-label-col">
            <span class="facet-label">${t('stats.sourcesFilterLabel')}</span>
          </div>
          <div class="spatial-source-chips-list">
            ${sourceChipsHtml}
          </div>
        </div>
      </div>

      <!-- Spatial Matrix Viewport -->
      <div class="telemetry-chart-viewport spatial-matrix-viewport">
        ${this.isSpatialLoading ? `<div class="timeline-loading-overlay">${t('stats.loadingTimeline')}</div>` : ''}
        
        <!-- Top Info Header & 5-Step Discrete Color Legend (Strict i18n) -->
        <div class="spatial-viewport-legend-row">
          <div class="spatial-legend-bar">
            <span class="legend-bar-label">${t('stats.scaleToggleLabel')}</span>
            <div class="legend-color-step">
              <span class="step-swatch ghost"></span>
              <span class="step-label">${t('stats.levelEmpty')}</span>
            </div>
            <div class="legend-color-step">
              <span class="step-swatch level-1"></span>
              <span class="step-label">${t('stats.levelLow')}</span>
            </div>
            <div class="legend-color-step">
              <span class="step-swatch level-2"></span>
              <span class="step-label">${t('stats.levelMedium')}</span>
            </div>
            <div class="legend-color-step">
              <span class="step-swatch level-3"></span>
              <span class="step-label">${t('stats.levelHigh')}</span>
            </div>
            <div class="legend-color-step">
              <span class="step-swatch level-4"></span>
              <span class="step-label">${t('stats.levelPeak', { count: maxCellHits.toLocaleString() })}</span>
            </div>
          </div>

          <div class="spatial-legend-notes">
            <span class="spatial-ghost-note">${t('stats.ghostCellsNote')}</span>
            <span class="spatial-quadrant-note">${t('stats.quadrantLinesNote')}</span>
          </div>
        </div>

        <!-- Canvas Container Slot -->
        <div id="spatial-matrix-slot"></div>
      </div>
    `;
  }

  private attachSpatialEventListeners(): void {
    // Preset buttons
    this.element.querySelectorAll('[data-spatial-preset]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const p = (e.currentTarget as HTMLElement).getAttribute('data-spatial-preset') as TelemetryPreset;
        if (p) {
          this.applySpatialPreset(p);
        }
      });
    });

    // Date range inputs
    const startInput = this.element.querySelector('#spatial-date-start') as HTMLInputElement | null;
    const endInput = this.element.querySelector('#spatial-date-end') as HTMLInputElement | null;

    const onDateInputChange = () => {
      if (startInput && endInput) {
        if (startInput.value && endInput.value && startInput.value <= endInput.value) {
          this.spatialHasCustomRange = true;
          this.spatialCustomStartDate = startInput.value;
          this.spatialCustomEndDate = endInput.value;
          this.fetchSpatialMatrix();
        }
      }
    };

    startInput?.addEventListener('change', onDateInputChange);
    endInput?.addEventListener('change', onDateInputChange);

    // Metric toggle
    this.element.querySelector('#spatial-metric-points-btn')?.addEventListener('click', () => {
      if (this.spatialMetric !== 'points') {
        this.spatialMetric = 'points';
        this.fetchSpatialMatrix();
      }
    });

    this.element.querySelector('#spatial-metric-calls-btn')?.addEventListener('click', () => {
      if (this.spatialMetric !== 'calls') {
        this.spatialMetric = 'calls';
        this.fetchSpatialMatrix();
      }
    });

    // Scale toggle
    this.element.querySelector('#spatial-scale-log-btn')?.addEventListener('click', () => {
      if (this.spatialScaleMode !== 'log') {
        this.spatialScaleMode = 'log';
        this.render();
        this.mountSpatialChart();
      }
    });

    this.element.querySelector('#spatial-scale-linear-btn')?.addEventListener('click', () => {
      if (this.spatialScaleMode !== 'linear') {
        this.spatialScaleMode = 'linear';
        this.render();
        this.mountSpatialChart();
      }
    });

    // Source chips
    this.element.querySelectorAll('[data-spatial-source]').forEach((chip) => {
      chip.addEventListener('click', (e) => {
        const srcName = (e.currentTarget as HTMLElement).getAttribute('data-spatial-source');
        if (srcName) {
          if (this.spatialActiveSourceKeys.has(srcName)) {
            this.spatialActiveSourceKeys.delete(srcName);
          } else {
            this.spatialActiveSourceKeys.add(srcName);
          }
          this.fetchSpatialMatrix();
        }
      });
    });
  }

  private mountSpatialChart(): void {
    const slot = this.element.querySelector('#spatial-matrix-slot');
    if (!slot) return;
    slot.innerHTML = '';
    slot.appendChild(this.spatialChart.getElement());

    if (this.spatialData) {
      this.spatialChart.update(this.spatialData, {
        metric: this.spatialMetric,
        scaleMode: this.spatialScaleMode,
      });
    }
  }
}
