import { t, onLanguageChange } from '../i18n';
import { formatTelemetryNumber, formatTelemetryDecimal } from '../utils/numberFormat';
import {
  apiClient,
  type GridInfoResponse,
  type TelemetryStatsResponse,
  type TimelineResponse,
} from '../utils/apiClient';
import { TimelineChart } from './TimelineChart';

export type TelemetryPreset = 'all' | '12m' | 'ytd' | '30d' | '7d';
export type Granularity = 'daily' | 'monthly';
export type SourceFilter = 'all' | 'web_map' | 'rest_api' | 'desktop_cli';

export class StatsView {
  private container: HTMLElement;
  private backdrop: HTMLElement | null = null;
  private gridInfo: GridInfoResponse | null = null;
  private telemetry: TelemetryStatsResponse | null = null;
  private timelineData: TimelineResponse | null = null;

  // State
  private activeSubTab: 'timeline' | 'grid' = 'timeline';
  private activePreset: TelemetryPreset = 'all';
  private activeGranularity: Granularity = 'daily';
  private activeSource: SourceFilter = 'all';
  private isLoadingTimeline: boolean = false;

  private timelineChart: TimelineChart;

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'tab-content';
    this.timelineChart = new TimelineChart();

    this.render();

    onLanguageChange(() => {
      this.render();
      if (this.backdrop) this.renderModal();
    });
    this.loadData();
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public async open(): Promise<void> {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';
    root.appendChild(this.backdrop);

    this.renderModal();
    await Promise.all([this.loadData(), this.fetchTimeline()]);
    this.renderModal();
  }

  public close(): void {
    if (this.backdrop) {
      this.backdrop.remove();
      this.backdrop = null;
    }
  }

  private async loadData(): Promise<void> {
    try {
      const [gridRes, teleRes] = await Promise.all([
        apiClient.getGridInfo(),
        apiClient.getTelemetryStats(),
      ]);
      this.gridInfo = gridRes;
      this.telemetry = teleRes;
    } catch (err) {
      console.warn('Failed to load telemetry / grid metadata', err);
    } finally {
      this.render();
    }
  }

  private async fetchTimeline(): Promise<void> {
    this.isLoadingTimeline = true;
    try {
      const now = new Date();
      let fromDate: string | undefined = undefined;
      const toDate = now.toISOString().slice(0, 10);

      if (this.activePreset === '7d') {
        const d = new Date();
        d.setDate(d.getDate() - 7);
        fromDate = d.toISOString().slice(0, 10);
      } else if (this.activePreset === '30d') {
        const d = new Date();
        d.setDate(d.getDate() - 30);
        fromDate = d.toISOString().slice(0, 10);
      } else if (this.activePreset === 'ytd') {
        fromDate = `${now.getFullYear()}-01-01`;
      } else if (this.activePreset === '12m') {
        const d = new Date();
        d.setFullYear(d.getFullYear() - 1);
        fromDate = d.toISOString().slice(0, 10);
      }

      let sourceId: number | undefined = undefined;
      if (this.activeSource === 'web_map') sourceId = 2;
      else if (this.activeSource === 'rest_api') sourceId = 1;
      else if (this.activeSource === 'desktop_cli') sourceId = 3;

      const res = await apiClient.getTelemetryTimeline({
        from_date: fromDate,
        to_date: toDate,
        granularity: this.activeGranularity,
        source: sourceId,
      });

      this.timelineData = res;
    } catch (err) {
      console.warn('Failed to load telemetry timeline:', err);
    } finally {
      this.isLoadingTimeline = false;
      if (this.backdrop && this.activeSubTab === 'timeline') {
        this.updateChartAndKPIs();
      }
    }
  }

  private renderModal(): void {
    if (!this.backdrop) return;

    this.backdrop.innerHTML = `
      <div class="modal-dialog telemetry-dialog">
        <div class="modal-header telemetry-modal-header">
          <div class="telemetry-header-title-group">
            <div class="modal-title">${t('stats.title')}</div>
            <div class="telemetry-subtabs">
              <button class="telemetry-subtab-btn ${this.activeSubTab === 'timeline' ? 'active' : ''}" id="tab-btn-timeline">
                📊 ${t('stats.tabTimeline')}
              </button>
              <button class="telemetry-subtab-btn ${this.activeSubTab === 'grid' ? 'active' : ''}" id="tab-btn-grid">
                🌐 ${t('stats.tabGrid')}
              </button>
            </div>
          </div>
          <button class="btn-icon" id="stats-close-icon" aria-label="Close">✕</button>
        </div>

        <div class="modal-body telemetry-modal-body">
          ${this.activeSubTab === 'timeline' ? this.renderTimelineTabHtml() : this.renderGridTabHtml()}
        </div>

        <div class="modal-footer">
          <button class="btn btn-secondary" id="stats-close-btn">${t('stats.closeBtn') || 'Close'}</button>
        </div>
      </div>
    `;

    // Event Listeners for Tabs
    this.backdrop.querySelector('#tab-btn-timeline')?.addEventListener('click', () => {
      if (this.activeSubTab !== 'timeline') {
        this.activeSubTab = 'timeline';
        this.renderModal();
        if (!this.timelineData) this.fetchTimeline();
        else this.updateChartAndKPIs();
      }
    });

    this.backdrop.querySelector('#tab-btn-grid')?.addEventListener('click', () => {
      if (this.activeSubTab !== 'grid') {
        this.activeSubTab = 'grid';
        this.renderModal();
      }
    });

    this.backdrop.querySelector('#stats-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#stats-close-btn')?.addEventListener('click', () => this.close());
    this.backdrop.addEventListener('click', (e) => {
      if (e.target === this.backdrop) this.close();
    });

    if (this.activeSubTab === 'timeline') {
      this.attachTimelineControls();
      this.mountChart();
    }
  }

  private renderTimelineTabHtml(): string {
    const summary = this.timelineData?.summary;
    const totalCalls = summary?.total_calls ?? (this.telemetry?.summary.total_calls ?? 0);
    const totalPoints = summary?.total_points ?? (this.telemetry?.summary.total_points ?? 0);
    const avgRatio = typeof summary?.avg_points_per_call === 'number'
      ? formatTelemetryDecimal(summary.avg_points_per_call, 1)
      : (totalCalls > 0 ? formatTelemetryDecimal(totalPoints / totalCalls, 1) : '0');
    const periodText = summary?.from_date && summary?.to_date
      ? `${summary.from_date} → ${summary.to_date} (${summary.days_count ?? ''} ${t('stats.daysCount', { count: summary.days_count ?? 0 })})`
      : 'All Logged History';

    return `
      <!-- KPI Banner -->
      <div class="telemetry-kpi-banner">
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val blue" id="kpi-calls">${formatTelemetryNumber(totalCalls)}</div>
          <div class="telemetry-kpi-label">${t('stats.totalCalls')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val cyan" id="kpi-points">${formatTelemetryNumber(totalPoints)}</div>
          <div class="telemetry-kpi-label">${t('stats.totalPoints')}</div>
        </div>
        <div class="telemetry-kpi-card">
          <div class="telemetry-kpi-val emerald" id="kpi-ratio">${avgRatio}</div>
          <div class="telemetry-kpi-label">${t('stats.batchIntensity')}</div>
        </div>
        <div class="telemetry-kpi-card wide">
          <div class="telemetry-kpi-val muted date-kpi" id="kpi-period">${periodText}</div>
          <div class="telemetry-kpi-label">${t('stats.activePeriod')}</div>
        </div>
      </div>

      <!-- Controls & Filter Toolbar -->
      <div class="telemetry-toolbar">
        <div class="telemetry-controls-group">
          <!-- Presets -->
          <div class="telemetry-btn-toggle-group">
            <button class="telemetry-toggle-btn ${this.activePreset === '7d' ? 'active' : ''}" data-preset="7d">${t('stats.preset7D')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === '30d' ? 'active' : ''}" data-preset="30d">${t('stats.preset30D')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === 'ytd' ? 'active' : ''}" data-preset="ytd">${t('stats.presetYTD')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === '12m' ? 'active' : ''}" data-preset="12m">${t('stats.preset12M')}</button>
            <button class="telemetry-toggle-btn ${this.activePreset === 'all' ? 'active' : ''}" data-preset="all">${t('stats.presetAll')}</button>
          </div>

          <!-- Granularity -->
          <div class="telemetry-btn-toggle-group">
            <button class="telemetry-toggle-btn ${this.activeGranularity === 'daily' ? 'active' : ''}" data-granularity="daily">${t('stats.granularityDaily')}</button>
            <button class="telemetry-toggle-btn ${this.activeGranularity === 'monthly' ? 'active' : ''}" data-granularity="monthly">${t('stats.granularityMonthly')}</button>
          </div>
        </div>

        <!-- Sources Filter -->
        <div class="telemetry-source-chips">
          <span class="telemetry-chips-label">${t('stats.sourceFilter')}</span>
          <button class="filter-chip ${this.activeSource === 'all' ? 'active' : ''}" data-source="all">${t('stats.sourceAll')}</button>
          <button class="filter-chip ${this.activeSource === 'web_map' ? 'active' : ''}" data-source="web_map">${t('stats.sourceWebMap')}</button>
          <button class="filter-chip ${this.activeSource === 'rest_api' ? 'active' : ''}" data-source="rest_api">${t('stats.sourceRestApi')}</button>
          <button class="filter-chip ${this.activeSource === 'desktop_cli' ? 'active' : ''}" data-source="desktop_cli">${t('stats.sourceDesktopCli')}</button>
        </div>
      </div>

      <!-- Chart Area -->
      <div class="telemetry-chart-viewport">
        ${this.isLoadingTimeline ? `<div class="timeline-loading-overlay">${t('stats.loadingTimeline')}</div>` : ''}
        <div class="chart-legend-row">
          <div class="legend-badge"><span class="legend-box blue"></span> ${t('stats.callsLegend')}</div>
          <div class="legend-badge"><span class="legend-box cyan"></span> ${t('stats.pointsLegend')}</div>
        </div>
        <div id="timeline-chart-slot"></div>
      </div>
    `;
  }

  private renderGridTabHtml(): string {
    return `
      <!-- Grid Info Section -->
      <div>
        <div class="section-micro-label">
          ${t('stats.activeGridTitle')}
        </div>
        <div class="content-card">
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.gridFile')}</span>
            <span class="mono">${this.gridInfo?.grid_file || 'transdat.spg'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.quasigeoid')}</span>
            <span>${this.gridInfo?.quasigeoid_model || 'Model Cvasigeoid 2008'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.projCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.projected || 'EPSG:31700 (Stereo 70)'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.geoCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.geographic || 'EPSG:4258 (ETRS89)'}</span>
          </div>
        </div>
      </div>

      <!-- Accuracy Notice -->
      <div class="callout-box accent-primary">
        <span class="callout-title">${t('stats.accuracyNoteTitle')}</span>
        <p class="callout-body">${t('stats.accuracyNoteText')}</p>
      </div>
    `;
  }

  private attachTimelineControls(): void {
    if (!this.backdrop) return;

    // Presets
    this.backdrop.querySelectorAll('[data-preset]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const val = (e.currentTarget as HTMLElement).getAttribute('data-preset') as TelemetryPreset;
        if (val && val !== this.activePreset) {
          this.activePreset = val;
          this.renderModal();
          this.fetchTimeline();
        }
      });
    });

    // Granularity
    this.backdrop.querySelectorAll('[data-granularity]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const val = (e.currentTarget as HTMLElement).getAttribute('data-granularity') as Granularity;
        if (val && val !== this.activeGranularity) {
          this.activeGranularity = val;
          this.renderModal();
          this.fetchTimeline();
        }
      });
    });

    // Source Filter
    this.backdrop.querySelectorAll('[data-source]').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        const val = (e.currentTarget as HTMLElement).getAttribute('data-source') as SourceFilter;
        if (val && val !== this.activeSource) {
          this.activeSource = val;
          this.renderModal();
          this.fetchTimeline();
        }
      });
    });
  }

  private mountChart(): void {
    if (!this.backdrop) return;
    const slot = this.backdrop.querySelector('#timeline-chart-slot');
    if (!slot) return;
    slot.innerHTML = '';
    slot.appendChild(this.timelineChart.getElement());

    if (this.timelineData?.series) {
      const activeKeys = new Set(
        this.activeSource === 'all'
          ? ['rest_api', 'web_map', 'desktop_cli']
          : [this.activeSource]
      );
      this.timelineChart.update(this.timelineData.series, {
        activeSourceKeys: activeKeys,
        availableSources: this.timelineData.available_sources || [],
        showTotal: true,
      });
    }
  }

  private updateChartAndKPIs(): void {
    if (!this.backdrop || this.activeSubTab !== 'timeline') return;
    this.renderModal();
  }

  private render(): void {
    const totalCalls = this.telemetry?.summary.total_calls ?? 0;
    const totalPoints = this.telemetry?.summary.total_points ?? 0;

    this.container.innerHTML = `
      <div class="studio-card">
        <div class="card-title">
          <span>${t('stats.activeGridTitle')}</span>
        </div>
        <div class="content-card">
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.gridFile')}</span>
            <span class="mono">${this.gridInfo?.grid_file || 'transdat.spg'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.quasigeoid')}</span>
            <span>${this.gridInfo?.quasigeoid_model || 'Model Cvasigeoid 2008'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.projCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.projected || 'EPSG:31700 (Stereo 70)'}</span>
          </div>
          <div class="content-row-divided">
            <span style="color: var(--text-muted);">${t('stats.geoCRS')}</span>
            <span class="mono">${this.gridInfo?.crs.geographic || 'EPSG:4258 (ETRS89)'}</span>
          </div>
        </div>

        <div class="callout-box accent-primary">
          <span class="callout-title">${t('stats.accuracyNoteTitle')}</span>
          <p class="callout-body">${t('stats.accuracyNoteText')}</p>
        </div>
      </div>

      <div class="studio-card">
        <div class="card-title">
          <span>${t('stats.telemetryTitle')}</span>
        </div>
        <div class="status-grid">
          <div class="status-badge-card">
            <div class="status-badge-num">${formatTelemetryNumber(totalCalls)}</div>
            <div class="status-badge-label">${t('stats.totalCalls')}</div>
          </div>
          <div class="status-badge-card success">
            <div class="status-badge-num">${formatTelemetryNumber(totalPoints)}</div>
            <div class="status-badge-label">${t('stats.totalPoints')}</div>
          </div>
        </div>
      </div>
    `;
  }
}
