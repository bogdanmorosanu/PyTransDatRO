import noUiSlider, { type API as NoUiSliderAPI } from 'nouislider';
import 'nouislider/dist/nouislider.css';
import { t, onLanguageChange } from '../i18n';
import {
  apiClient,
  type GridInfoResponse,
  type TelemetryStatsResponse,
  type TimelineResponse,
  type SourceItem,
} from '../utils/apiClient';
import {
  TimelineChart,
  TOTAL_SERIES_CONFIG,
} from './TimelineChart';

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

  // View state
  private activeSubTab: 'timeline' | 'grid' = 'timeline';
  private activePreset: TelemetryPreset = 'all';
  private activeGranularity: Granularity = 'daily';

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

  constructor(callbacks: TelemetryDashboardCallbacks) {
    this.callbacks = callbacks;
    this.element = document.createElement('div');
    this.element.className = 'telemetry-page-view';
    this.chart = new TimelineChart();

    this.render();

    onLanguageChange(() => {
      this.render();
      this.mountChart();
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
              📊 ${t('stats.tabTimeline')}
            </button>
            <button class="telemetry-subtab-btn ${this.activeSubTab === 'grid' ? 'active' : ''}" id="tab-grid-btn">
              🌐 ${t('stats.tabGrid')}
            </button>
          </div>
        </div>

        <!-- Main Body -->
        <div class="telemetry-page-content">
          ${this.activeSubTab === 'timeline' ? this.renderTimelineHtml() : this.renderGridHtml()}
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

    this.element.querySelector('#tab-grid-btn')?.addEventListener('click', () => {
      this.activeSubTab = 'grid';
      this.render();
    });

    if (this.activeSubTab === 'timeline') {
      this.attachTimelineEventListeners();
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
}
