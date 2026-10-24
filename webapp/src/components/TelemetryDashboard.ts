import { formatTelemetryNumber, formatTelemetryDecimal, formatTelemetryPercent, formatDisplayDate } from '../utils/numberFormat';
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
import { LegacyTimelineChart, LEGACY_SERIES_DEFS } from './LegacyTimelineChart';
import type { LegacyTimelineResponse, LegacyStatsResponse } from '../utils/apiClient';

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
  private activeSubTab: 'timeline' | 'spatial' | 'legacy' = 'timeline';

  // Legacy view state
  private legacyData: LegacyTimelineResponse | null = null;
  private legacyStats: LegacyStatsResponse | null = null;
  private legacyPreset: TelemetryPreset = 'all';
  private legacyGranularity: Granularity = 'daily';
  private legacyActiveSources: Set<string> = new Set(['app', 'service']);
  private legacyActiveOps: Set<string> = new Set(['st70', 'st30']);
  private legacyShowTotal: boolean = true;
  private isLegacyLoading: boolean = false;
  private legacyMinDate: string = '2011-01-09';
  private legacyMaxDate: string = '2026-10-14';
  private legacyCustomStartDate: string = '2011-01-09';
  private legacyCustomEndDate: string = '2026-10-14';
  private legacyHasCustomRange: boolean = false;
  private legacyChart: LegacyTimelineChart;
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
  private compShowTotal: boolean = true;

  private chart: TimelineChart;
  private spatialChart: SpatialMatrixChart;

  constructor(_callbacks?: TelemetryDashboardCallbacks) {
    void _callbacks;
    this.element = document.createElement('div');
    this.element.className = 'telemetry-page-view';
    this.chart = new TimelineChart();
    this.spatialChart = new SpatialMatrixChart();
    this.legacyChart = new LegacyTimelineChart();

    this.render();

    onLanguageChange(() => {
      this.render();
      if (this.activeSubTab === 'timeline') this.mountChart();
      if (this.activeSubTab === 'spatial') this.mountSpatialChart();
      if (this.activeSubTab === 'legacy') this.mountLegacyChart();
    });

    window.addEventListener('themechange', () => {
      this.render();
      if (this.activeSubTab === 'timeline') this.mountChart();
      if (this.activeSubTab === 'spatial') this.mountSpatialChart();
      if (this.activeSubTab === 'legacy') this.mountLegacyChart();
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
        let fromDate = res.summary.from_date;
        if (fromDate.length === 7) {
          fromDate = `${fromDate}-01`;
        }
        this.minHistoryDate = fromDate;
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
              <button 
                type="button" 
                class="telemetry-nav-tab-btn ${this.activeSubTab === 'legacy' ? 'active' : ''}" 
                id="tab-legacy-btn"
                role="tab"
                aria-selected="${this.activeSubTab === 'legacy'}"
              >
                <span>${t('stats.tabLegacy')}</span>
                <span class="telemetry-tab-indicator"></span>
              </button>
            </nav>
          </div>
        </div>

        <!-- Main Body -->
        <div class="telemetry-page-content">
          ${this.activeSubTab === 'timeline' ? this.renderTimelineHtml() : (this.activeSubTab === 'spatial' ? this.renderSpatialHtml() : this.renderLegacyHtml())}
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

    this.element.querySelector('#tab-legacy-btn')?.addEventListener('click', () => {
      this.activeSubTab = 'legacy';
      this.render();
      this.mountLegacyChart();
      if (!this.legacyData) {
        this.fetchLegacyData();
      }
    });

    if (this.activeSubTab === 'timeline') {
      this.attachTimelineEventListeners();
    } else if (this.activeSubTab === 'spatial') {
      this.attachSpatialEventListeners();
    } else {
      this.attachLegacyEventListeners();
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

    const avgRatio = totalCalls > 0 ? formatTelemetryDecimal(totalPoints / totalCalls, 1) : '0';

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

    const compEffectiveFrom = this.timelineData?.summary?.effective_from_date ?? this.timelineData?.summary?.from_date ?? this.customStartDate;
    const compEffectiveTo = this.timelineData?.summary?.effective_to_date ?? this.timelineData?.summary?.to_date ?? this.customEndDate;
    const compMinFormatted = formatDisplayDate(this.minHistoryDate.length === 7 ? `${this.minHistoryDate}-01` : this.minHistoryDate);
    const compMaxFormatted = formatDisplayDate(this.maxHistoryDate);
    const compActiveStartFormatted = formatDisplayDate(compEffectiveFrom);
    const compActiveEndFormatted = formatDisplayDate(compEffectiveTo);

    return `
      <!-- KPI Banner: Neutral Theme Colors, Clean Info Presentation -->
      <div class="telemetry-kpi-banner">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(totalPoints)}</div>
          <div class="telemetry-kpi-label">${t('stats.totalPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(totalCalls)}</div>
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
        <div class="telemetry-deck-row unified-controls-row flex-between-center">
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

        <!-- Row 2: Native Themed Calendar Range Controls, Bounds (Left) & Active Period (Right) -->
        <div class="telemetry-deck-row telemetry-date-controls-row">
          <!-- Left Group: User Query Selection & Historical Limits -->
          <div class="telemetry-date-left-group">
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

            <div class="telemetry-date-limits-cluster">
              <div class="telemetry-info-subgroup">
                <span class="info-label">${t('stats.earliestAvailable')}</span>
                <span class="info-field-readonly">${compMinFormatted}</span>
              </div>
              <div class="telemetry-info-subgroup">
                <span class="info-label">${t('stats.currentDate')}</span>
                <span class="info-field-readonly">${compMaxFormatted}</span>
              </div>
            </div>
          </div>

          <!-- Vertical Visual Divider -->
          <div class="telemetry-deck-vdivider" aria-hidden="true"></div>

          <!-- Right Group: Database Effective Period Result -->
          <div class="telemetry-date-active-cluster" id="comp-effective-period-group">
            <span class="form-label telemetry-active-main-label">${t('stats.effectivePeriod')}</span>
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.periodStart')}</span>
              <span class="info-field-readonly active-period-value">${compActiveStartFormatted}</span>
            </div>
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.periodEnd')}</span>
              <span class="info-field-readonly active-period-value">${compActiveEndFormatted}</span>
            </div>
            <span class="info-help-indicator" id="comp-active-period-info-btn" aria-label="${t('stats.activePeriodTooltip')}" tabindex="0">
              <svg width="13" height="13" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                <circle cx="8" cy="8" r="7"/>
                <line x1="8" y1="7.5" x2="8" y2="11.5"/>
                <circle cx="8" cy="4.5" r="0.75" fill="currentColor"/>
              </svg>
              <span class="info-tooltip-bubble" role="tooltip">${t('stats.activePeriodTooltip')}</span>
            </span>
          </div>
        </div>

        <!-- Row 3: Merged Sources (Left) & Metrics (Right) -->
        <div class="telemetry-deck-row merged-filters-row">
          <div class="filters-left-group">
            <label class="form-label">${t('stats.sourcesFilterLabel')}</label>
            <div class="telemetry-btn-toggle-group multi-select">
              ${sourceTogglesHtml}
            </div>
          </div>
          <div class="filters-right-group">
            <label class="form-label">${t('stats.metricsFilterLabel')}</label>
            <div class="telemetry-btn-toggle-group multi-select">
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

    // Top-Left Dedicated Total Cluster (Stationary Anchor)
    let totalHtml = '';
    if (this.activeSourceKeys.size > 0) {
      const totalColors = this.chart.getTotalSeriesColors();
      totalHtml = `
        <div class="legend-total-cluster">
          <button type="button" class="telemetry-toggle-btn individual compact ${this.compShowTotal ? 'active' : ''}" id="comp-toggle-total">
            ${t('stats.seriesTotal')}
          </button>
          ${
            this.compShowTotal
              ? `
            <div class="legend-total-swatches">
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
          `
              : ''
          }
        </div>
      `;
    }

    return `
      <div class="chart-legend-container">
        ${totalHtml}
        <div class="legend-sources-cluster">
          ${sourceGroups.join('')}
        </div>
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

    // Total series toggle
    this.element.querySelector('#comp-toggle-total')?.addEventListener('click', () => {
      this.toggleCompTotal();
    });

    this.element.addEventListener('click', () => {
      this.element.querySelectorAll('.legend-color-popover.open').forEach((p) => p.classList.remove('open'));
    });
  }

  private toggleCompTotal(): void {
    this.compShowTotal = !this.compShowTotal;
    this.render();
    this.mountChart();
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
        showTotal: this.compShowTotal,
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

    const spEffectiveFrom = this.spatialData?.summary?.effective_from_date ?? this.spatialData?.summary?.from_date ?? this.spatialCustomStartDate;
    const spEffectiveTo = this.spatialData?.summary?.effective_to_date ?? this.spatialData?.summary?.to_date ?? this.spatialCustomEndDate;
    const spMinFormatted = formatDisplayDate(this.minHistoryDate.length === 7 ? `${this.minHistoryDate}-01` : this.minHistoryDate);
    const spMaxFormatted = formatDisplayDate(this.maxHistoryDate);
    const spActiveStartFormatted = formatDisplayDate(spEffectiveFrom);
    const spActiveEndFormatted = formatDisplayDate(spEffectiveTo);

    return `
      <!-- Overall Stats Area: 2 Archetype 8.1 Standalone KPI Cards -->
      <div class="telemetry-kpi-banner two-cards">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(totalVal)}</div>
          <div class="telemetry-kpi-label">${this.spatialMetric === 'points' ? t('stats.totalPoints') : t('stats.totalCalls')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(activeCells)}</div>
          <div class="telemetry-kpi-label">${t('stats.activeCellsLabel')}</div>
        </div>
      </div>

      <!-- 4 Mathematical Quadrant Cards (Name Above, Value Below with adjacent %) -->
      <div class="spatial-quadrants-deck">
        <div class="quadrant-card nw">
          <div class="quadrant-name">${t('stats.quadrantNW')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${formatTelemetryNumber(qNW.value)}</span>
            <span class="quadrant-share mono">(${formatTelemetryPercent(qNW.percentage, 1)})</span>
          </div>
        </div>

        <div class="quadrant-card ne">
          <div class="quadrant-name">${t('stats.quadrantNE')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${formatTelemetryNumber(qNE.value)}</span>
            <span class="quadrant-share mono">(${formatTelemetryPercent(qNE.percentage, 1)})</span>
          </div>
        </div>

        <div class="quadrant-card sw">
          <div class="quadrant-name">${t('stats.quadrantSW')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${formatTelemetryNumber(qSW.value)}</span>
            <span class="quadrant-share mono">(${formatTelemetryPercent(qSW.percentage, 1)})</span>
          </div>
        </div>

        <div class="quadrant-card se">
          <div class="quadrant-name">${t('stats.quadrantSE')}</div>
          <div class="quadrant-val-row">
            <span class="quadrant-value mono">${formatTelemetryNumber(qSE.value)}</span>
            <span class="quadrant-share mono">(${formatTelemetryPercent(qSE.percentage, 1)})</span>
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

        <!-- Row 2: Native Themed Calendar Inputs, Bounds (Left) & Active Period (Right) -->
        <div class="telemetry-deck-row telemetry-date-controls-row">
          <!-- Left Group: User Query Selection & Historical Limits -->
          <div class="telemetry-date-left-group">
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

            <div class="telemetry-date-limits-cluster">
              <div class="telemetry-info-subgroup">
                <span class="info-label">${t('stats.earliestAvailable')}</span>
                <span class="info-field-readonly">${spMinFormatted}</span>
              </div>
              <div class="telemetry-info-subgroup">
                <span class="info-label">${t('stats.currentDate')}</span>
                <span class="info-field-readonly">${spMaxFormatted}</span>
              </div>
            </div>
          </div>

          <!-- Vertical Visual Divider -->
          <div class="telemetry-deck-vdivider" aria-hidden="true"></div>

          <!-- Right Group: Database Effective Period Result -->
          <div class="telemetry-date-active-cluster" id="spatial-effective-period-group">
            <span class="form-label telemetry-active-main-label">${t('stats.effectivePeriod')}</span>
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.periodStart')}</span>
              <span class="info-field-readonly active-period-value">${spActiveStartFormatted}</span>
            </div>
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.periodEnd')}</span>
              <span class="info-field-readonly active-period-value">${spActiveEndFormatted}</span>
            </div>
            <span class="info-help-indicator" id="spatial-active-period-info-btn" aria-label="${t('stats.activePeriodTooltip')}" tabindex="0">
              <svg width="13" height="13" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                <circle cx="8" cy="8" r="7"/>
                <line x1="8" y1="7.5" x2="8" y2="11.5"/>
                <circle cx="8" cy="4.5" r="0.75" fill="currentColor"/>
              </svg>
              <span class="info-tooltip-bubble" role="tooltip">${t('stats.activePeriodTooltip')}</span>
            </span>
          </div>
        </div>

        <!-- Row 3: Merged SOURCES (Left) & METRIC (Right) -->
        <div class="telemetry-deck-row merged-filters-row">
          <div class="filters-left-group">
            <label class="form-label">${t('stats.sourcesFilterLabel')}</label>
            <div class="telemetry-btn-toggle-group multi-select">
              ${spatialSourceTogglesHtml}
            </div>
          </div>
          <div class="filters-right-group">
            <label class="form-label">${t('stats.metricFilterLabel')}</label>
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
              <span class="step-label">${t('stats.levelPeak', { count: formatTelemetryNumber(maxCellHits) })}</span>
            </div>
          </div>

          <!-- Bottom-Right Scaling Selector -->
          <div class="spatial-scaling-selector-bottom">
            <label class="form-label">${t('stats.scaleToggleLabel')}</label>
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

  // ==========================================================================
  // Legacy Subpage Methods
  // ==========================================================================
  private async fetchLegacyData(): Promise<void> {
    this.isLegacyLoading = true;
    this.render();
    try {
      const [statsRes, timelineRes] = await Promise.all([
        this.legacyStats ? Promise.resolve(this.legacyStats) : apiClient.getLegacyStats(),
        apiClient.getLegacyTimeline({
          granularity: this.legacyGranularity,
          from_date: (this.legacyPreset !== 'all' || this.legacyHasCustomRange) ? this.legacyCustomStartDate : undefined,
          to_date: (this.legacyPreset !== 'all' || this.legacyHasCustomRange) ? this.legacyCustomEndDate : undefined,
          source: Array.from(this.legacyActiveSources).join(','),
          operation: Array.from(this.legacyActiveOps).join(','),
        }),
      ]);

      this.legacyStats = statsRes;
      this.legacyData = timelineRes;

      if (timelineRes.bounds.min_date) {
        this.legacyMinDate = timelineRes.bounds.min_date;
      }
      if (timelineRes.bounds.max_date) {
        this.legacyMaxDate = timelineRes.bounds.max_date;
        if (!this.legacyHasCustomRange && this.legacyPreset === 'all') {
          this.legacyCustomStartDate = this.legacyMinDate;
          this.legacyCustomEndDate = this.legacyMaxDate;
        }
      }
    } catch (err) {
      console.warn('Failed to fetch legacy timeline data:', err);
    } finally {
      this.isLegacyLoading = false;
      this.render();
      this.mountLegacyChart();
    }
  }

  private applyLegacyPreset(preset: TelemetryPreset): void {
    this.legacyPreset = preset;
    this.legacyHasCustomRange = false;

    // Use maxHistoryDate from DB (2026-10-14), not the live browser clock!
    const baseDate = new Date(this.legacyMaxDate);
    this.legacyCustomEndDate = this.legacyMaxDate;

    if (preset === '7d') {
      const d = new Date(baseDate);
      d.setDate(d.getDate() - 7);
      this.legacyCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === '30d') {
      const d = new Date(baseDate);
      d.setDate(d.getDate() - 30);
      this.legacyCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'ytd') {
      this.legacyCustomStartDate = `${baseDate.getFullYear()}-01-01`;
    } else if (preset === '12m') {
      const d = new Date(baseDate);
      d.setFullYear(d.getFullYear() - 1);
      this.legacyCustomStartDate = d.toISOString().slice(0, 10);
    } else if (preset === 'all') {
      this.legacyCustomStartDate = this.legacyMinDate;
    }

    this.fetchLegacyData();
  }

  private toggleLegacySource(sourceKey: string): void {
    if (this.legacyActiveSources.has(sourceKey)) {
      if (this.legacyActiveSources.size > 1) {
        this.legacyActiveSources.delete(sourceKey);
      }
    } else {
      this.legacyActiveSources.add(sourceKey);
    }
    this.fetchLegacyData();
  }

  private toggleLegacyOp(opKey: string): void {
    if (this.legacyActiveOps.has(opKey)) {
      if (this.legacyActiveOps.size > 1) {
        this.legacyActiveOps.delete(opKey);
      }
    } else {
      this.legacyActiveOps.add(opKey);
    }
    this.fetchLegacyData();
  }

  private toggleLegacyTotal(): void {
    this.legacyShowTotal = !this.legacyShowTotal;
    this.render();
    this.mountLegacyChart();
  }

  private mountLegacyChart(): void {
    const slot = this.element.querySelector('#legacy-chart-slot');
    if (!slot) return;
    slot.innerHTML = '';
    slot.appendChild(this.legacyChart.getElement());

    // Build active keys based on selected sources and ops
    const activeKeys = new Set<string>();
    for (const src of this.legacyActiveSources) {
      for (const op of this.legacyActiveOps) {
        activeKeys.add(`${src}_${op}`);
      }
    }

    this.legacyChart.update(this.legacyData?.series || [], {
      activeKeys,
      showTotal: this.legacyShowTotal,
    });
  }

  private renderLegacyHtml(): string {
    let totalPoints = 0;
    let appPoints = 0;
    let servicePoints = 0;
    let st70Points = 0;
    let st30Points = 0;

    if (this.legacyData?.series) {
      for (const item of this.legacyData.series) {
        totalPoints += item.points.total;
        if (this.legacyActiveSources.has('app')) {
          if (this.legacyActiveOps.has('st70')) appPoints += item.points.app_st70;
          if (this.legacyActiveOps.has('st30')) appPoints += item.points.app_st30;
        }
        if (this.legacyActiveSources.has('service')) {
          if (this.legacyActiveOps.has('st70')) servicePoints += item.points.service_st70;
          if (this.legacyActiveOps.has('st30')) servicePoints += item.points.service_st30;
        }
        if (this.legacyActiveOps.has('st70')) {
          if (this.legacyActiveSources.has('app')) st70Points += item.points.app_st70;
          if (this.legacyActiveSources.has('service')) st70Points += item.points.service_st70;
        }
        if (this.legacyActiveOps.has('st30')) {
          if (this.legacyActiveSources.has('app')) st30Points += item.points.app_st30;
          if (this.legacyActiveSources.has('service')) st30Points += item.points.service_st30;
        }
      }
    } else if (this.legacyStats) {
      totalPoints = this.legacyStats.total_points;
      appPoints = this.legacyStats.app_st70 + this.legacyStats.app_st30;
      servicePoints = this.legacyStats.service_st70 + this.legacyStats.service_st30;
      st70Points = this.legacyStats.app_st70 + this.legacyStats.service_st70;
      st30Points = this.legacyStats.app_st30 + this.legacyStats.service_st30;
    }

    const hasApp = this.legacyActiveSources.has('app');
    const hasService = this.legacyActiveSources.has('service');
    const hasSt70 = this.legacyActiveOps.has('st70');
    const hasSt30 = this.legacyActiveOps.has('st30');

    return `
      <!-- KPI Banner: Points Only -->
      <div class="telemetry-kpi-banner five-cards">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(totalPoints)}</div>
          <div class="telemetry-kpi-label">${t('stats.legacyTotalPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(appPoints)}</div>
          <div class="telemetry-kpi-label">${t('stats.legacyAppPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(servicePoints)}</div>
          <div class="telemetry-kpi-label">${t('stats.legacyServicePoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(st70Points)}</div>
          <div class="telemetry-kpi-label">${t('stats.legacyStereo70Points')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val neutral-val">${formatTelemetryNumber(st30Points)}</div>
          <div class="telemetry-kpi-label">${t('stats.legacyStereo30Points')}</div>
        </div>
      </div>

      <!-- Controls Deck -->
      <div class="telemetry-deck">
        <!-- Row 1: Unified Presets & Granularity -->
        <div class="telemetry-deck-row unified-controls-row flex-between-center">
          <div class="telemetry-facet-group">
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.legacyPreset === 'all' && !this.legacyHasCustomRange ? 'active' : ''}" data-legacy-preset="all">${t('stats.presetAll')}</button>
              <button class="telemetry-toggle-btn ${this.legacyPreset === '12m' && !this.legacyHasCustomRange ? 'active' : ''}" data-legacy-preset="12m">${t('stats.preset12M')}</button>
              <button class="telemetry-toggle-btn ${this.legacyPreset === 'ytd' && !this.legacyHasCustomRange ? 'active' : ''}" data-legacy-preset="ytd">${t('stats.presetYTD')}</button>
              <button class="telemetry-toggle-btn ${this.legacyPreset === '30d' && !this.legacyHasCustomRange ? 'active' : ''}" data-legacy-preset="30d">${t('stats.preset30D')}</button>
              <button class="telemetry-toggle-btn ${this.legacyPreset === '7d' && !this.legacyHasCustomRange ? 'active' : ''}" data-legacy-preset="7d">${t('stats.preset7D')}</button>
            </div>
          </div>

          <div class="telemetry-facet-group">
            <div class="telemetry-btn-toggle-group">
              <button class="telemetry-toggle-btn ${this.legacyGranularity === 'daily' ? 'active' : ''}" data-legacy-granularity="daily">${t('stats.granularityDaily')}</button>
              <button class="telemetry-toggle-btn ${this.legacyGranularity === 'weekly' ? 'active' : ''}" data-legacy-granularity="weekly">${t('stats.granularityWeekly')}</button>
              <button class="telemetry-toggle-btn ${this.legacyGranularity === 'monthly' ? 'active' : ''}" data-legacy-granularity="monthly">${t('stats.granularityMonthly')}</button>
            </div>
          </div>
        </div>

        <!-- Row 2: Fixed Archive Calendar Range Controls & Range Limits Info -->
        <div class="telemetry-deck-row telemetry-date-controls-row">
          <div class="telemetry-date-inputs-cluster">
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="legacy-date-start"
                min="${this.legacyMinDate}" 
                max="${this.legacyMaxDate}" 
                value="${this.legacyCustomStartDate}" 
              />
            </div>
            <span class="range-arrow">→</span>
            <div class="date-input-group">
              <input 
                type="date" 
                class="telemetry-date-input" 
                id="legacy-date-end"
                min="${this.legacyMinDate}" 
                max="${this.legacyMaxDate}" 
                value="${this.legacyCustomEndDate}" 
              />
            </div>
          </div>

          <div class="telemetry-date-limits-cluster">
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.archiveStart')}</span>
              <span class="info-field-readonly">${this.legacyMinDate}</span>
            </div>
            <div class="telemetry-info-subgroup">
              <span class="info-label">${t('stats.archiveEnd')}</span>
              <span class="info-field-readonly">${this.legacyMaxDate}</span>
            </div>
          </div>
        </div>

        <!-- Row 3: Faceted Filters: Source (Left) & Operation (Right) -->
        <div class="telemetry-deck-row merged-filters-row">
          <div class="filters-left-group">
            <label class="form-label">${t('stats.sourcesFilterLabel')}</label>
            <div class="telemetry-btn-toggle-group multi-select">
              <button type="button" class="telemetry-toggle-btn individual ${hasApp ? 'active' : ''}" data-legacy-source="app">
                ${t('stats.legacySourceApp')}
              </button>
              <button type="button" class="telemetry-toggle-btn individual ${hasService ? 'active' : ''}" data-legacy-source="service">
                ${t('stats.legacySourceService')}
              </button>
            </div>
          </div>

          <div class="filters-right-group">
            <label class="form-label">${t('stats.legacyOpLabel')}</label>
            <div class="telemetry-btn-toggle-group multi-select">
              <button type="button" class="telemetry-toggle-btn individual ${hasSt70 ? 'active' : ''}" data-legacy-op="st70">
                ${t('stats.legacyOpStereo70')}
              </button>
              <button type="button" class="telemetry-toggle-btn individual ${hasSt30 ? 'active' : ''}" data-legacy-op="st30">
                ${t('stats.legacyOpStereo30')}
              </button>
            </div>
          </div>
        </div>
      </div>

      <!-- Chart Viewport -->
      <div class="telemetry-chart-viewport full-page">
        ${this.isLegacyLoading ? `<div class="timeline-loading-overlay">${t('stats.loadingLegacyTimeline')}</div>` : ''}
        ${this.renderLegacyLegendHtml()}
        <div id="legacy-chart-slot"></div>
      </div>
    `;
  }

  private renderLegacyLegendHtml(): string {
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

    const sourceGroups: string[] = [];

    for (const src of this.legacyActiveSources) {
      const srcName = src === 'app' ? t('stats.legacySourceApp') : t('stats.legacySourceService');
      const opItems: string[] = [];

      for (const op of this.legacyActiveOps) {
        const seriesKey = `${src}_${op}`;
        const def = LEGACY_SERIES_DEFS[seriesKey];
        if (!def) continue;
        const color = this.legacyChart.getLineColor(seriesKey, def.defaultColor);
        const opLabel = op === 'st70' ? t('stats.legacyOpStereo70') : t('stats.legacyOpStereo30');

        opItems.push(`
          <span class="legend-sub-item interactive-legend-item" data-legend-series="${seriesKey}" title="Click to change line color">
            <span class="legend-line-sample" style="background-color: ${color};"></span>
            <span>${opLabel}</span>
            ${renderSwatchPopover(seriesKey, color)}
          </span>
        `);
      }

      sourceGroups.push(`
        <div class="legend-source-group">
          <span class="legend-group-title">${srcName}:</span>
          <div class="legend-group-items">
            ${opItems.join('')}
          </div>
        </div>
      `);
    }

    let totalHtml = '';
    if (this.legacyActiveSources.size > 0 && this.legacyActiveOps.size > 0) {
      const totalColor = this.legacyChart.getTotalColor();
      totalHtml = `
        <div class="legend-total-cluster">
          <button type="button" class="telemetry-toggle-btn individual compact ${this.legacyShowTotal ? 'active' : ''}" id="legacy-toggle-total">
            ${t('stats.seriesTotal')}
          </button>
          ${
            this.legacyShowTotal
              ? `
            <div class="legend-total-swatches">
              <span class="legend-sub-item interactive-legend-item" data-legend-series="total" title="Click to change line color">
                <span class="legend-line-sample dashed-sample" style="background: repeating-linear-gradient(90deg, ${totalColor} 0, ${totalColor} 7px, transparent 7px, transparent 12px);"></span>
                <span>${t('stats.legacyTotalPoints')}</span>
                ${renderSwatchPopover('total', totalColor)}
              </span>
            </div>
          `
              : ''
          }
        </div>
      `;
    }

    return `
      <div class="chart-legend-container">
        ${totalHtml}
        <div class="legend-sources-row">
          ${sourceGroups.join('')}
        </div>
      </div>
    `;
  }

  private attachLegacyEventListeners(): void {
    // Preset buttons
    this.element.querySelectorAll('[data-legacy-preset]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const preset = (e.currentTarget as HTMLElement).getAttribute('data-legacy-preset') as TelemetryPreset;
        if (preset) this.applyLegacyPreset(preset);
      });
    });

    // Granularity buttons
    this.element.querySelectorAll('[data-legacy-granularity]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const gran = (e.currentTarget as HTMLElement).getAttribute('data-legacy-granularity') as Granularity;
        if (gran) {
          this.legacyGranularity = gran;
          this.fetchLegacyData();
        }
      });
    });

    // Source buttons
    this.element.querySelectorAll('[data-legacy-source]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const src = (e.currentTarget as HTMLElement).getAttribute('data-legacy-source');
        if (src) this.toggleLegacySource(src);
      });
    });

    // Op buttons
    this.element.querySelectorAll('[data-legacy-op]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const op = (e.currentTarget as HTMLElement).getAttribute('data-legacy-op');
        if (op) this.toggleLegacyOp(op);
      });
    });

    // Total overlay toggle
    this.element.querySelector('#legacy-toggle-total')?.addEventListener('click', () => {
      this.toggleLegacyTotal();
    });

    // Date range inputs
    const startInput = this.element.querySelector('#legacy-date-start') as HTMLInputElement;
    const endInput = this.element.querySelector('#legacy-date-end') as HTMLInputElement;

    const onDateChange = () => {
      if (startInput && endInput && startInput.value && endInput.value) {
        if (startInput.value <= endInput.value) {
          this.legacyCustomStartDate = startInput.value;
          this.legacyCustomEndDate = endInput.value;
          this.legacyHasCustomRange = true;
          this.fetchLegacyData();
        }
      }
    };

    startInput?.addEventListener('change', onDateChange);
    endInput?.addEventListener('change', onDateChange);

    // Interactive Legend Color Popover Toggles & Swatch Selection for Legacy
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
          this.legacyChart.setLineColor(seriesKey, color);
          this.render();
          this.mountLegacyChart();
        }
      });
    });

    this.element.addEventListener('click', () => {
      this.element.querySelectorAll('.legend-color-popover.open').forEach((p) => p.classList.remove('open'));
    });
  }
}