import { t, onLanguageChange } from '../i18n';
import {
  apiClient,
  type TelemetryStatsResponse,
  type TimelineResponse,
  type SourceItem,
  type SpatialMatrixResponse,
} from '../utils/apiClient';
import {
  TimelineChart,
  TIMELINE_COLOR_PRESETS,
  type TimelineMetricFilter,
} from './TimelineChart';
import { SpatialMatrixChart } from './SpatialMatrixChart';

export type TelemetryPreset = 'all' | '12m' | 'ytd' | '30d' | '7d';
export type Granularity = 'daily' | 'weekly' | 'monthly';

export interface TelemetryDashboardCallbacks {
  onBackToStudio?: () => void;
}

// Canonical lookup dictionary for source display names
export const SOURCE_LABEL_LOOKUP: Record<string, string> = {
  rest_api: 'REST API',
  web_app: 'Web App',
  web_map: 'Web App', // backward compatibility alias
  desktop_cli: 'Desktop / CLI',
};

export class TelemetryDashboard {
  private element: HTMLElement;

  private telemetry: TelemetryStatsResponse | null = null;
  private timelineData: TimelineResponse | null = null;
  private spatialData: SpatialMatrixResponse | null = null;

  // View state
  private activeSubTab: 'timeline' | 'spatial' = 'timeline';
  private activePreset: TelemetryPreset = 'all';
  private activeGranularity: Granularity = 'daily';
  private compShowPoints: boolean = true;
  private compShowCalls: boolean = true;

  public get computationMetricFilter(): TimelineMetricFilter {
    if (this.compShowPoints && this.compShowCalls) return 'both';
    if (this.compShowPoints) return 'points';
    if (this.compShowCalls) return 'calls';
    return 'both';
  }

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

  // Current selected range for computation
  private hasCustomRange: boolean = false;
  private customStartDate: string = '2023-11-01';
  private customEndDate: string = new Date().toISOString().slice(0, 10);

  private chart: TimelineChart;
  private spatialChart: SpatialMatrixChart;

  constructor(_callbacks?: TelemetryDashboardCallbacks) {
    void _callbacks;
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

  public getSourceDisplayName(sourceName: string): string {
    if (SOURCE_LABEL_LOOKUP[sourceName]) {
      return SOURCE_LABEL_LOOKUP[sourceName];
    }
    // Return exact DB text without transformation
    return sourceName;
  }

  private async initData(): Promise<void> {
    this.isLoading = true;
    try {
      const teleRes = await apiClient.getTelemetryStats();
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
            { id: 2, name: 'web_app', description: 'Web App' },
            { id: 3, name: 'desktop_cli', description: 'Desktop / CLI' },
          ];
        }
      }

      // Initialize active toggles: activate all non-unknown sources
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

      if (this.spatialActiveSourceKeys.size === 0) {
        this.availableSources.forEach((s) => {
          if (s.name !== 'unknown') {
            this.spatialActiveSourceKeys.add(s.name);
          }
        });
        if (this.spatialActiveSourceKeys.size === 0 && this.availableSources.length > 0) {
          this.spatialActiveSourceKeys.add(this.availableSources[0].name);
        }
      }

      // Track earliest available date from database
      if (res.summary.from_date && this.activePreset === 'all') {
        this.minHistoryDate = res.summary.from_date;
        if (!this.hasCustomRange) {
          this.customStartDate = this.minHistoryDate;
        }
        if (!this.spatialHasCustomRange && this.spatialPreset === 'all') {
          this.spatialCustomStartDate = this.minHistoryDate;
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
      if (this.activeSourceKeys.size > 1) {
        this.activeSourceKeys.delete(sourceKey);
      }
    } else {
      this.activeSourceKeys.add(sourceKey);
    }
    this.render();
    this.mountChart();
  }

  private toggleSpatialSource(sourceKey: string): void {
    if (this.spatialActiveSourceKeys.has(sourceKey)) {
      if (this.spatialActiveSourceKeys.size > 1) {
        this.spatialActiveSourceKeys.delete(sourceKey);
      }
    } else {
      this.spatialActiveSourceKeys.add(sourceKey);
    }
    this.fetchSpatialMatrix();
  }

  private render(): void {
    this.element.innerHTML = `
      <div class="telemetry-page-container">
        <!-- Top Navigation Bar -->
        <div class="telemetry-page-header">
          <div class="telemetry-title-block">
            <div class="telemetry-page-title">${t('stats.title')}</div>
            <div class="header-brand-nav-divider"></div>
            <nav class="telemetry-nav-tabs" role="tablist">
              <button 
                type="button" 
                class="telemetry-nav-tab-btn ${this.activeSubTab === 'timeline' ? 'active' : ''}" 
                id="tab-timeline-btn"
                role="tab"
                aria-selected="${this.activeSubTab === 'timeline'}"
              >
                <span>${t('stats.tabComputation')}</span>
                <span class="telemetry-tab-indicator"></span>
              </button>
              <button 
                type="button" 
                class="telemetry-nav-tab-btn ${this.activeSubTab === 'spatial' ? 'active' : ''}" 
                id="tab-spatial-btn"
                role="tab"
                aria-selected="${this.activeSubTab === 'spatial'}"
              >
                <span>${t('stats.tabSpatial')}</span>
                <span class="telemetry-tab-indicator"></span>
              </button>
            </nav>
          </div>
        </div>

        <!-- Main Body -->
        <div class="telemetry-page-content">
          ${this.activeSubTab === 'timeline' ? this.renderTimelineHtml() : this.renderSpatialHtml()}
        </div>
      </div>
    `;

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

    if (this.activeSubTab === 'timeline') {
      this.attachTimelineEventListeners();
    } else {
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

    // Filter out unknown source from UI
    const displaySources = this.availableSources.filter((s) => s.name !== 'unknown');

    // Dynamic Source Toggles (Individual toggle buttons styled like theme-native toggle buttons)
    const sourceTogglesHtml = displaySources
      .map((s) => {
        const isActive = this.activeSourceKeys.has(s.name);
        const displayName = this.getSourceDisplayName(s.name);
        return `
          <button type="button" class="telemetry-toggle-btn individual ${isActive ? 'active' : ''}" data-source-name="${s.name}">
            ${displayName}
          </button>
        `;
      })
      .join('');

    return `
      <!-- KPI Banner: Neutral Theme Colors, Clean Info Presentation -->
      <div class="telemetry-kpi-banner">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${totalPoints.toLocaleString()}</div>
          <div class="telemetry-kpi-label">${t('stats.totalPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${totalCalls.toLocaleString()}</div>
          <div class="telemetry-kpi-label">${t('stats.totalCalls')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${avgRatio}</div>
          <div class="telemetry-kpi-label">${t('stats.batchIntensity')}</div>
        </div>
      </div>

      <!-- Controls Deck -->
      <div class="telemetry-deck">
        <!-- Row 1: Unified Single Bar Presets (Left) & Granularity (Right) -->
        <div class="telemetry-deck-row unified-controls-row" style="justify-content: space-between;">
          <div class="telemetry-facet-group">
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.activePreset === 'all' && !this.hasCustomRange ? 'active' : ''}" data-preset="all">${t('stats.presetAll')}</button>
              <button class="telemetry-toggle-btn ${this.activePreset === '12m' && !this.hasCustomRange ? 'active' : ''}" data-preset="12m">${t('stats.preset12M')}</button>
              <button class="telemetry-toggle-btn ${this.activePreset === 'ytd' && !this.hasCustomRange ? 'active' : ''}" data-preset="ytd">${t('stats.presetYTD')}</button>
              <button class="telemetry-toggle-btn ${this.activePreset === '30d' && !this.hasCustomRange ? 'active' : ''}" data-preset="30d">${t('stats.preset30D')}</button>
              <button class="telemetry-toggle-btn ${this.activePreset === '7d' && !this.hasCustomRange ? 'active' : ''}" data-preset="7d">${t('stats.preset7D')}</button>
            </div>
          </div>

          <div class="telemetry-facet-group">
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.activeGranularity === 'daily' ? 'active' : ''}" data-granularity="daily">${t('stats.granularityDaily')}</button>
              <button class="telemetry-toggle-btn ${this.activeGranularity === 'weekly' ? 'active' : ''}" data-granularity="weekly">${t('stats.granularityWeekly')}</button>
              <button class="telemetry-toggle-btn ${this.activeGranularity === 'monthly' ? 'active' : ''}" data-granularity="monthly">${t('stats.granularityMonthly')}</button>
            </div>
          </div>
        </div>

        <!-- Row 2: Native Themed Calendar Range Controls & Range Limits Info -->
        <div class="telemetry-deck-row telemetry-date-controls-row">
          <div class="telemetry-date-inputs-cluster">
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="comp-date-start"
                min="${this.minHistoryDate}" 
                max="${this.maxHistoryDate}" 
                value="${this.customStartDate}" 
              />
            </div>
            <span class="range-arrow">→</span>
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="comp-date-end"
                min="${this.minHistoryDate}" 
                max="${this.maxHistoryDate}" 
                value="${this.customEndDate}" 
              />
            </div>
          </div>

          <div class="telemetry-range-limits-pill">
            <span class="limit-subtle">${t('stats.sliderEarliest', { date: this.minHistoryDate })}</span>
            <span class="limit-divider">•</span>
            <span class="limit-subtle">${t('stats.sliderToday', { date: this.maxHistoryDate })}</span>
          </div>
        </div>

        <!-- Row 3: Merged SOURCES (Left) & METRIC (Right) -->
        <div class="telemetry-deck-row merged-filters-row">
          <div class="filters-left-group">
            <span class="facet-label">${t('stats.sourcesFilterLabel')}</span>
            <div class="sources-toggle-list">
              ${sourceTogglesHtml}
            </div>
          </div>
          <div class="filters-right-group">
            <span class="facet-label">${t('stats.metricFilterLabel')}</span>
            <div class="telemetry-individual-toggle-group">
              <button class="telemetry-toggle-btn individual ${this.compShowPoints ? 'active' : ''}" id="comp-metric-points-btn">
                ${t('stats.metricPoints')}
              </button>
              <button class="telemetry-toggle-btn individual ${this.compShowCalls ? 'active' : ''}" id="comp-metric-calls-btn">
                ${t('stats.metricCalls')}
              </button>
            </div>
          </div>
        </div>
      </div>

      <!-- Chart Viewport with Dedicated Graph Symbology Legend -->
      <div class="telemetry-chart-viewport full-page">
        ${this.isLoading ? `<div class="timeline-loading-overlay">${t('stats.loadingTimeline')}</div>` : ''}
        
        <!-- Dedicated Graph Symbology Legend (2-Line Grouped: Line 1 = Sources, Line 2 = Total) -->
        ${this.renderTimelineLegendHtml()}

        <div id="timeline-chart-slot"></div>
      </div>
    `;
  }

  private renderTimelineLegendHtml(): string {
    const showCalls = this.compShowCalls;
    const showPoints = this.compShowPoints;

    const sourceGroups: string[] = [];

    const renderSwatchPopover = (seriesKey: string, currentColor: string): string => {
      const swatches = TIMELINE_COLOR_PRESETS.map((hex) => {
        const isActive = currentColor.toLowerCase() === hex.toLowerCase();
        return `<button type="button" class="legend-swatch-option ${isActive ? 'active' : ''}" data-series-key="${seriesKey}" data-swatch-color="${hex}" style="background-color: ${hex};" title="${hex}"></button>`;
      }).join('');
      return `
        <div class="legend-color-popover" role="menu">
          ${swatches}
        </div>
      `;
    };

    // Active source lines grouped by source name
    let idx = 0;
    for (const srcKey of this.activeSourceKeys) {
      if (srcKey === 'unknown') continue;
      const config = this.chart.getSourceConfig(srcKey, idx++);
      const displayName = this.getSourceDisplayName(srcKey);

      sourceGroups.push(`
        <div class="legend-source-group">
          <span class="legend-group-title">${displayName}:</span>
          <div class="legend-group-items">
            ${showCalls ? `
              <span class="legend-sub-item interactive-legend-item" data-legend-series="${srcKey}:calls" title="Click to change line color">
                <span class="legend-line-sample" style="background-color: ${config.callsColor};"></span>
                <span>${t('stats.metricTypeCalls')}</span>
                ${renderSwatchPopover(`${srcKey}:calls`, config.callsColor)}
              </span>
            ` : ''}
            ${showPoints ? `
              <span class="legend-sub-item interactive-legend-item" data-legend-series="${srcKey}:points" title="Click to change line color">
                <span class="legend-line-sample" style="background-color: ${config.pointsColor};"></span>
                <span>${t('stats.metricTypePoints')}</span>
                ${renderSwatchPopover(`${srcKey}:points`, config.pointsColor)}
              </span>
            ` : ''}
          </div>
        </div>
      `);
    }

    // Line 2: Dedicated Total row
    let totalHtml = '';
    if (this.activeSourceKeys.size > 0) {
      const totalColors = this.chart.getTotalSeriesColors();
      totalHtml = `
        <div class="legend-total-row">
          <span class="legend-group-title">${t('stats.seriesTotal')}:</span>
          <div class="legend-group-items">
            ${showCalls ? `
              <span class="legend-sub-item interactive-legend-item" data-legend-series="total:calls" title="Click to change line color">
                <span class="legend-line-sample dashed-sample" style="background: repeating-linear-gradient(90deg, ${totalColors.callsColor} 0, ${totalColors.callsColor} 7px, transparent 7px, transparent 12px);"></span>
                <span>${t('stats.metricTypeCalls')}</span>
                ${renderSwatchPopover('total:calls', totalColors.callsColor)}
              </span>
            ` : ''}
            ${showPoints ? `
              <span class="legend-sub-item interactive-legend-item" data-legend-series="total:points" title="Click to change line color">
                <span class="legend-line-sample dashed-sample" style="background: repeating-linear-gradient(90deg, ${totalColors.pointsColor} 0, ${totalColors.pointsColor} 7px, transparent 7px, transparent 12px);"></span>
                <span>${t('stats.metricTypePoints')}</span>
                ${renderSwatchPopover('total:points', totalColors.pointsColor)}
              </span>
            ` : ''}
          </div>
        </div>
      `;
    }

    return `
      <div class="chart-legend-container">
        <div class="legend-sources-row">
          ${sourceGroups.join('')}
        </div>
        ${totalHtml}
      </div>
    `;
  }

  private attachTimelineEventListeners(): void {
    // Preset buttons
    this.element.querySelectorAll('[data-preset]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const p = (e.currentTarget as HTMLElement).getAttribute('data-preset') as TelemetryPreset;
        if (p) this.applyPreset(p);
      });
    });

    // Granularity buttons
    this.element.querySelectorAll('[data-granularity]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const g = (e.currentTarget as HTMLElement).getAttribute('data-granularity') as Granularity;
        if (g && g !== this.activeGranularity) {
          this.activeGranularity = g;
          this.fetchTimeline();
        }
      });
    });

    // Calendar Inputs with strict clamping to database boundaries
    const startInput = this.element.querySelector('#comp-date-start') as HTMLInputElement | null;
    const endInput = this.element.querySelector('#comp-date-end') as HTMLInputElement | null;

    const onDateInputChange = () => {
      if (startInput && endInput) {
        let startVal = startInput.value;
        let endVal = endInput.value;

        if (startVal < this.minHistoryDate) {
          startVal = this.minHistoryDate;
          startInput.value = this.minHistoryDate;
        }
        if (endVal > this.maxHistoryDate) {
          endVal = this.maxHistoryDate;
          endInput.value = this.maxHistoryDate;
        }
        if (startVal > endVal) {
          endVal = startVal;
          endInput.value = startVal;
        }

        this.hasCustomRange = true;
        this.customStartDate = startVal;
        this.customEndDate = endVal;
        this.fetchTimeline();
      }
    };

    startInput?.addEventListener('change', onDateInputChange);
    endInput?.addEventListener('change', onDateInputChange);

    // Metric Multi-Toggle Buttons (Points & Calls can both be active)
    this.element.querySelector('#comp-metric-points-btn')?.addEventListener('click', () => {
      if (this.compShowPoints && this.compShowCalls) {
        this.compShowPoints = false;
      } else if (!this.compShowPoints) {
        this.compShowPoints = true;
      }
      this.render();
      this.mountChart();
    });

    this.element.querySelector('#comp-metric-calls-btn')?.addEventListener('click', () => {
      if (this.compShowCalls && this.compShowPoints) {
        this.compShowCalls = false;
      } else if (!this.compShowCalls) {
        this.compShowCalls = true;
      }
      this.render();
      this.mountChart();
    });

    // Source Toggle Buttons
    this.element.querySelectorAll<HTMLButtonElement>('[data-source-name]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const target = e.currentTarget as HTMLButtonElement;
        const srcName = target.getAttribute('data-source-name');
        if (srcName) {
          this.toggleSource(srcName);
        }
      });
    });

    // Interactive Legend Color Popover Toggles & Swatch Selection
    this.element.querySelectorAll<HTMLElement>('.interactive-legend-item').forEach((item) => {
      item.addEventListener('click', (e) => {
        e.stopPropagation();
        const popover = item.querySelector('.legend-color-popover');
        if (!popover) return;
        const wasOpen = popover.classList.contains('open');
        this.element.querySelectorAll('.legend-color-popover.open').forEach((p) => p.classList.remove('open'));
        if (!wasOpen) {
          popover.classList.add('open');
        }
      });
    });

    this.element.querySelectorAll<HTMLButtonElement>('.legend-swatch-option').forEach((swatch) => {
      swatch.addEventListener('click', (e) => {
        e.stopPropagation();
        const seriesKey = swatch.getAttribute('data-series-key');
        const color = swatch.getAttribute('data-swatch-color');
        if (seriesKey && color) {
          this.chart.setLineColor(seriesKey, color);
          this.render();
          this.mountChart();
        }
      });
    });

    this.element.addEventListener('click', () => {
      this.element.querySelectorAll('.legend-color-popover.open').forEach((p) => p.classList.remove('open'));
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
        metricFilter: this.computationMetricFilter,
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
      if (this.spatialActiveSourceKeys.size === 0 && this.availableSources.length > 0) {
        this.availableSources.forEach((s) => {
          if (s.name !== 'unknown') {
            this.spatialActiveSourceKeys.add(s.name);
          }
        });
      }

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

    const qNW = quadrants?.NW || { value: 0, percentage: 0 };
    const qNE = quadrants?.NE || { value: 0, percentage: 0 };
    const qSW = quadrants?.SW || { value: 0, percentage: 0 };
    const qSE = quadrants?.SE || { value: 0, percentage: 0 };

    // Dynamic Source Toggles matching Computation lookup and omitting unknown
    const displaySources = this.availableSources.filter((s) => s.name !== 'unknown');
    const spatialSourceTogglesHtml = displaySources
      .map((s) => {
        const isActive = this.spatialActiveSourceKeys.size === 0 || this.spatialActiveSourceKeys.has(s.name);
        const displayName = this.getSourceDisplayName(s.name);
        return `
          <button type="button" class="telemetry-toggle-btn individual ${isActive ? 'active' : ''}" data-spatial-source-name="${s.name}">
            ${displayName}
          </button>
        `;
      })
      .join('');

    return `
      <!-- Overall Stats Area: Full-Width Summary Rectangle -->
      <div class="spatial-overview-banner">
        <div class="spatial-overview-col left-col">
          <div class="overview-sublabel">${this.spatialMetric === 'points' ? t('stats.totalPoints') : t('stats.totalCalls')}</div>
          <div class="overview-main-val mono">${totalVal.toLocaleString()}</div>
        </div>
        <div class="spatial-overview-col right-col">
          <div class="overview-sublabel">${t('stats.activeCellsLabel')}</div>
          <div class="overview-main-val mono">${activeCells.toLocaleString()}</div>
        </div>
      </div>

      <!-- 4 Mathematical Quadrant Cards (Name Above, Value Below with adjacent %) -->
      <div class="spatial-quadrants-deck">
        <div class="quadrant-card nw">
          <div class="quadrant-name">${t('stats.quadrantNW')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${qNW.value.toLocaleString()}</span>
            <span class="quadrant-share mono">(${qNW.percentage}%)</span>
          </div>
        </div>

        <div class="quadrant-card ne">
          <div class="quadrant-name">${t('stats.quadrantNE')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${qNE.value.toLocaleString()}</span>
            <span class="quadrant-share mono">(${qNE.percentage}%)</span>
          </div>
        </div>

        <div class="quadrant-card sw">
          <div class="quadrant-name">${t('stats.quadrantSW')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${qSW.value.toLocaleString()}</span>
            <span class="quadrant-share mono">(${qSW.percentage}%)</span>
          </div>
        </div>

        <div class="quadrant-card se">
          <div class="quadrant-name">${t('stats.quadrantSE')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${qSE.value.toLocaleString()}</span>
            <span class="quadrant-share mono">(${qSE.percentage}%)</span>
          </div>
        </div>
      </div>

      <!-- Controls Deck (Compact, Matching Computation) -->
      <div class="telemetry-deck spatial-controls-deck">
        <!-- Row 1: Unified Presets Single Bar -->
        <div class="telemetry-deck-row unified-controls-row">
          <div class="telemetry-facet-group">
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.spatialPreset === 'all' && !this.spatialHasCustomRange ? 'active' : ''}" data-spatial-preset="all">
                ${t('stats.presetAll')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialPreset === '12m' && !this.spatialHasCustomRange ? 'active' : ''}" data-spatial-preset="12m">
                ${t('stats.preset12M')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialPreset === 'ytd' && !this.spatialHasCustomRange ? 'active' : ''}" data-spatial-preset="ytd">
                ${t('stats.presetYTD')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialPreset === '30d' && !this.spatialHasCustomRange ? 'active' : ''}" data-spatial-preset="30d">
                ${t('stats.preset30D')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialPreset === '7d' && !this.spatialHasCustomRange ? 'active' : ''}" data-spatial-preset="7d">
                ${t('stats.preset7D')}
              </button>
            </div>
          </div>
        </div>

        <!-- Row 2: Native Themed Calendar Inputs & Limits -->
        <div class="telemetry-deck-row telemetry-date-controls-row">
          <div class="telemetry-date-inputs-cluster">
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="spatial-date-start"
                min="${this.minHistoryDate}" 
                max="${this.maxHistoryDate}" 
                value="${this.spatialCustomStartDate}" 
              />
            </div>
            <span class="range-arrow">→</span>
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="spatial-date-end"
                min="${this.minHistoryDate}" 
                max="${this.maxHistoryDate}" 
                value="${this.spatialCustomEndDate}" 
              />
            </div>
          </div>

          <div class="telemetry-range-limits-pill">
            <span class="limit-subtle">${t('stats.sliderEarliest', { date: this.minHistoryDate })}</span>
            <span class="limit-divider">•</span>
            <span class="limit-subtle">${t('stats.sliderToday', { date: this.maxHistoryDate })}</span>
          </div>
        </div>

        <!-- Row 3: Merged SOURCES (Left) & METRIC (Right) -->
        <div class="telemetry-deck-row merged-filters-row">
          <div class="filters-left-group">
            <span class="facet-label">${t('stats.sourcesFilterLabel')}</span>
            <div class="sources-toggle-list">
              ${spatialSourceTogglesHtml}
            </div>
          </div>
          <div class="filters-right-group">
            <span class="facet-label">${t('stats.metricFilterLabel')}</span>
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.spatialMetric === 'points' ? 'active' : ''}" id="spatial-metric-points-btn">
                ${t('stats.metricPoints')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialMetric === 'calls' ? 'active' : ''}" id="spatial-metric-calls-btn">
                ${t('stats.metricCalls')}
              </button>
            </div>
          </div>
        </div>
      </div>

      <!-- Spatial Matrix Viewport -->
      <div class="telemetry-chart-viewport spatial-matrix-viewport">
        ${this.isSpatialLoading ? `<div class="timeline-loading-overlay">${t('stats.loadingTimeline')}</div>` : ''}
        
        <!-- Canvas Container Slot -->
        <div id="spatial-matrix-slot"></div>

        <!-- Bottom Heatmap Controls & Legend Bar -->
        <div class="spatial-matrix-bottom-bar">
          <!-- 5-Step Discrete Color Ramp Legend -->
          <div class="spatial-legend-bar">
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

          <!-- Bottom-Right Scaling Selector -->
          <div class="spatial-scaling-selector-bottom">
            <span class="facet-label">${t('stats.scaleToggleLabel')}</span>
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.spatialScaleMode === 'log' ? 'active' : ''}" id="spatial-scale-log-btn">
                ${t('stats.scaleLog')}
              </button>
              <button class="telemetry-toggle-btn ${this.spatialScaleMode === 'linear' ? 'active' : ''}" id="spatial-scale-linear-btn">
                ${t('stats.scaleLinear')}
              </button>
            </div>
          </div>
        </div>

        <!-- Stacked Footnotes: Resolution, Cell Area, Quadrant Lines -->
        <div class="spatial-notes-block">
          <div class="spatial-note-item">${t('stats.gridCellSizeNote')}</div>
          <div class="spatial-note-item">${t('stats.quadrantLinesNote')}</div>
        </div>
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

    // Date range inputs with strict clamping to database boundaries
    const startInput = this.element.querySelector('#spatial-date-start') as HTMLInputElement | null;
    const endInput = this.element.querySelector('#spatial-date-end') as HTMLInputElement | null;

    const onDateInputChange = () => {
      if (startInput && endInput) {
        let startVal = startInput.value;
        let endVal = endInput.value;

        if (startVal < this.minHistoryDate) {
          startVal = this.minHistoryDate;
          startInput.value = this.minHistoryDate;
        }
        if (endVal > this.maxHistoryDate) {
          endVal = this.maxHistoryDate;
          endInput.value = this.maxHistoryDate;
        }
        if (startVal > endVal) {
          endVal = startVal;
          endInput.value = startVal;
        }

        this.spatialHasCustomRange = true;
        this.spatialCustomStartDate = startVal;
        this.spatialCustomEndDate = endVal;
        this.fetchSpatialMatrix();
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

    // Source toggle buttons
    this.element.querySelectorAll<HTMLButtonElement>('[data-spatial-source-name]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const target = e.currentTarget as HTMLButtonElement;
        const srcName = target.getAttribute('data-spatial-source-name');
        if (srcName) {
          this.toggleSpatialSource(srcName);
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
