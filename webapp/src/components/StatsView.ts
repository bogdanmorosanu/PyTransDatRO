import { t, onLanguageChange } from '../i18n';
import { apiClient, type GridInfoResponse, type TelemetryStatsResponse } from '../utils/apiClient';

export class StatsView {
  private container: HTMLElement;
  private backdrop: HTMLElement | null = null;
  private gridInfo: GridInfoResponse | null = null;
  private telemetry: TelemetryStatsResponse | null = null;

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'tab-content';
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
    await this.loadData();
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

  private renderModal(): void {
    if (!this.backdrop) return;
    const totalCalls = this.telemetry?.summary.total_calls ?? 0;
    const totalPoints = this.telemetry?.summary.total_points ?? 0;

    this.backdrop.innerHTML = `
      <div class="modal-dialog" style="max-width: 580px;">
        <div class="modal-header">
          <div class="modal-title">${t('stats.title')}</div>
          <button class="btn-icon" id="stats-close-icon" aria-label="Close">✕</button>
        </div>
        <div class="modal-body" style="display: flex; flex-direction: column; gap: 16px;">
          
          <!-- Grid Info Section -->
          <div>
            <div style="font-size: 12px; font-weight: 700; color: var(--color-primary-light); text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">
              ${t('stats.activeGridTitle')}
            </div>
            <div style="font-size: 13px; line-height: 1.8; background-color: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: var(--radius-md); padding: 10px 14px;">
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 4px 0;">
                <span style="color: var(--text-muted);">${t('stats.gridFile')}</span>
                <span class="mono" style="font-weight: 600; color: var(--color-primary-light);">${this.gridInfo?.grid_file || 'transdat.spg'}</span>
              </div>
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 4px 0;">
                <span style="color: var(--text-muted);">${t('stats.quasigeoid')}</span>
                <span style="font-weight: 500;">${this.gridInfo?.quasigeoid_model || 'Model Cvasigeoid 2008'}</span>
              </div>
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 4px 0;">
                <span style="color: var(--text-muted);">${t('stats.projCRS')}</span>
                <span class="mono" style="font-size: 12px;">${this.gridInfo?.crs.projected || 'EPSG:31700 (Stereo 70)'}</span>
              </div>
              <div style="display: flex; justify-content: space-between; padding: 4px 0;">
                <span style="color: var(--text-muted);">${t('stats.geoCRS')}</span>
                <span class="mono" style="font-size: 12px;">${this.gridInfo?.crs.geographic || 'EPSG:4258 (ETRS89)'}</span>
              </div>
            </div>
          </div>

          <!-- Accuracy Notice -->
          <div style="padding: 10px 14px; background-color: var(--bg-surface-input); border-left: 3px solid var(--color-primary); border-radius: var(--radius-sm);">
            <div style="font-weight: 600; font-size: 12px; margin-bottom: 3px; color: var(--color-primary-light);">
              ${t('stats.accuracyNoteTitle')}
            </div>
            <div style="font-size: 12px; color: var(--text-muted); line-height: 1.5;">
              ${t('stats.accuracyNoteText')}
            </div>
          </div>

          <!-- Telemetry Counters -->
          <div>
            <div style="font-size: 12px; font-weight: 700; color: var(--color-primary-light); text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">
              ${t('stats.telemetryTitle')}
            </div>
            <div class="status-grid" style="margin-bottom: 0;">
              <div class="status-badge-card">
                <div class="status-badge-num" style="color: var(--color-primary-light);">${totalCalls.toLocaleString()}</div>
                <div class="status-badge-label">${t('stats.totalCalls')}</div>
              </div>
              <div class="status-badge-card">
                <div class="status-badge-num" style="color: var(--color-success);">${totalPoints.toLocaleString()}</div>
                <div class="status-badge-label">${t('stats.totalPoints')}</div>
              </div>
            </div>
          </div>

        </div>
        <div class="modal-footer">
          <button class="btn btn-secondary" id="stats-close-btn">${t('stats.closeBtn') || 'Close'}</button>
        </div>
      </div>
    `;

    this.backdrop.querySelector('#stats-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#stats-close-btn')?.addEventListener('click', () => this.close());
    this.backdrop.addEventListener('click', (e) => {
      if (e.target === this.backdrop) this.close();
    });
  }

  private render(): void {
    const totalCalls = this.telemetry?.summary.total_calls ?? 0;
    const totalPoints = this.telemetry?.summary.total_points ?? 0;

    this.container.innerHTML = `
      <div class="studio-card">
        <div class="card-title">
          <span>${t('stats.activeGridTitle')}</span>
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

        <div style="margin-top: 14px; padding: 12px; background-color: var(--bg-surface-input); border-left: 3px solid var(--color-primary); border-radius: var(--radius-sm);">
          <div style="font-weight: 600; font-size: 12px; margin-bottom: 4px; color: var(--color-primary-light);">
            ${t('stats.accuracyNoteTitle')}
          </div>
          <div style="font-size: 12px; color: var(--text-muted); line-height: 1.6;">
            ${t('stats.accuracyNoteText')}
          </div>
        </div>
      </div>

      <div class="studio-card">
        <div class="card-title">
          <span>${t('stats.telemetryTitle')}</span>
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
      </div>
    `;
  }
}
