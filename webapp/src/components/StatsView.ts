import { t, onLanguageChange } from '../i18n';
import { apiClient, type GridInfoResponse, type TelemetryStatsResponse } from '../utils/apiClient';

export class StatsView {
  private container: HTMLElement;
  private gridInfo: GridInfoResponse | null = null;
  private telemetry: TelemetryStatsResponse | null = null;

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'tab-content';
    this.render();

    onLanguageChange(() => this.render());
    this.loadData();
  }

  public getElement(): HTMLElement {
    return this.container;
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

  private render(): void {
    const totalCalls = this.telemetry?.summary.total_calls ?? 0;
    const totalPoints = this.telemetry?.summary.total_points ?? 0;

    this.container.innerHTML = `
      <!-- Technical Geodetic Card -->
      <div class="studio-card">
        <div class="card-title">
          <span>🌐 ${t('stats.activeGridTitle')}</span>
        </div>
        <div style="font-size: 13px; line-height: 1.8;">
          <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
            <span style="color: var(--text-muted);">${t('stats.gridFile')}</span>
            <span class="mono" style="font-weight: 600; color: var(--color-primary-light);">${this.gridInfo?.grid_file || 'transdat.spg'}</span>
          </div>
          <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
            <span style="color: var(--text-muted);">${t('stats.quasigeoid')}</span>
            <span style="font-weight: 500;">${this.gridInfo?.quasigeoid_model || 'Model Cvasigeoid 2008'}</span>
          </div>
          <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
            <span style="color: var(--text-muted);">${t('stats.projCRS')}</span>
            <span class="mono" style="font-size: 12px;">${this.gridInfo?.crs.projected || 'EPSG:31700 (Stereo 70)'}</span>
          </div>
          <div style="display: flex; justify-content: space-between; padding: 6px 0;">
            <span style="color: var(--text-muted);">${t('stats.geoCRS')}</span>
            <span class="mono" style="font-size: 12px;">${this.gridInfo?.crs.geographic || 'EPSG:4258 (ETRS89)'}</span>
          </div>
        </div>

        <!-- Accuracy Notice -->
        <div style="margin-top: 14px; padding: 12px; background-color: var(--bg-surface-input); border-left: 3px solid var(--color-primary); border-radius: var(--radius-sm);">
          <div style="font-weight: 600; font-size: 12px; margin-bottom: 4px; color: var(--color-primary-light);">
            ⚖️ ${t('stats.accuracyNoteTitle')}
          </div>
          <div style="font-size: 12px; color: var(--text-muted); line-height: 1.6;">
            ${t('stats.accuracyNoteText')}
          </div>
        </div>
      </div>

      <!-- Telemetry Counters -->
      <div class="studio-card">
        <div class="card-title">
          <span>📈 ${t('stats.telemetryTitle')}</span>
        </div>
        <div class="status-grid">
          <div class="status-badge-card">
            <div class="status-badge-num" style="color: var(--color-primary-light);">${totalCalls.toLocaleString()}</div>
            <div class="status-badge-label">${t('stats.totalCalls')}</div>
          </div>
          <div class="status-badge-card">
            <div class="status-badge-num" style="color: var(--color-success);">${totalPoints.toLocaleString()}</div>
            <div class="status-badge-label">${t('stats.totalPoints')}</div>
          </div>
        </div>

        <!-- Daily Trends List -->
        ${
          this.telemetry && this.telemetry.daily_trend.length > 0
            ? `
          <div style="margin-top: 10px;">
            <div class="form-label">Ultimele 7 zile</div>
            <div style="font-size: 12px; font-family: var(--font-mono); max-height: 180px; overflow-y: auto;">
              ${this.telemetry.daily_trend
                .slice(0, 7)
                .map(
                  (d) => `
                <div style="display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px dashed var(--border-subtle);">
                  <span>${d.date} (${d.source})</span>
                  <span><strong>${d.points}</strong> pct / ${d.calls} req</span>
                </div>
              `
                )
                .join('')}
            </div>
          </div>
        `
            : ''
        }
      </div>
    `;
  }
}
