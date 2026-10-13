import { t } from '../i18n';
import { apiClient, type GridInfoResponse } from '../utils/apiClient';

export class GridInfoModal {
  private backdrop: HTMLElement | null = null;
  private gridInfo: GridInfoResponse | null = null;
  private isLoading: boolean = false;

  public async open(): Promise<void> {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';
    root.appendChild(this.backdrop);

    this.render();
    await this.loadGridData();
  }

  public close(): void {
    if (this.backdrop) {
      this.backdrop.remove();
      this.backdrop = null;
    }
  }

  private async loadGridData(): Promise<void> {
    this.isLoading = true;
    this.render();
    try {
      this.gridInfo = await apiClient.getGridInfo();
    } catch (err) {
      console.warn('Failed to fetch grid info for modal:', err);
    } finally {
      this.isLoading = false;
      this.render();
    }
  }

  private render(): void {
    if (!this.backdrop) return;

    const gridFile = this.gridInfo?.grid_file || (this.isLoading ? '...' : 'transdat.spg');
    const quasigeoid = this.gridInfo?.quasigeoid_model || (this.isLoading ? '...' : 'Model Cvasigeoid 2008');
    const projCrs = this.gridInfo?.crs?.projected || (this.isLoading ? '...' : 'EPSG:31700 (Stereo 70)');
    const geoCrs = this.gridInfo?.crs?.geographic || (this.isLoading ? '...' : 'EPSG:4258 (ETRS89)');

    this.backdrop.innerHTML = `
      <div class="modal-dialog" style="max-width: 580px;">
        <div class="modal-header">
          <div class="modal-title">${t('header.gridInfo')}</div>
          <button class="btn-icon" id="gridinfo-close-icon" aria-label="Close">✕</button>
        </div>
        <div class="modal-body" style="display: flex; flex-direction: column; gap: 16px; font-size: 13px; line-height: 1.6;">
          
          <div>
            <div style="font-size: 11px; font-weight: 700; color: var(--color-primary-light); text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 8px;">
              ${t('stats.activeGridTitle') || 'Active Transformation Grid'}
            </div>
            <div style="font-size: 13px; line-height: 1.8; background-color: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: var(--radius-md); padding: 12px 16px;">
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
                <span style="color: var(--text-muted);">${t('stats.gridFile') || 'Grid File'}</span>
                <span class="mono" style="font-weight: 600; color: var(--color-primary-light);">${gridFile}</span>
              </div>
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
                <span style="color: var(--text-muted);">${t('stats.quasigeoid') || 'Quasigeoid Model'}</span>
                <span style="font-weight: 500;">${quasigeoid}</span>
              </div>
              <div style="display: flex; justify-content: space-between; border-bottom: 1px solid var(--border-subtle); padding: 6px 0;">
                <span style="color: var(--text-muted);">${t('stats.projCRS') || 'Projected CRS'}</span>
                <span class="mono" style="font-size: 12px;">${projCrs}</span>
              </div>
              <div style="display: flex; justify-content: space-between; padding: 6px 0;">
                <span style="color: var(--text-muted);">${t('stats.geoCRS') || 'Geographic CRS'}</span>
                <span class="mono" style="font-size: 12px;">${geoCrs}</span>
              </div>
            </div>
          </div>

          <!-- Official Accuracy & Constraints Notice -->
          <div style="padding: 12px 14px; background-color: var(--bg-surface-input); border-left: 3px solid var(--color-primary); border-radius: var(--radius-sm);">
            <div style="font-weight: 600; font-size: 12px; margin-bottom: 4px; color: var(--color-primary-light);">
              ${t('stats.accuracyNoteTitle') || 'Standard Conformity & Limitations'}
            </div>
            <div style="font-size: 12px; color: var(--text-muted); line-height: 1.5;">
              ${t('stats.accuracyNoteText') || 'Transformations strictly follow official ANCPI grid boundaries and algorithms. Coordinates outside grid envelope are flagged.'}
            </div>
          </div>

        </div>
        <div class="modal-footer">
          <button class="btn btn-secondary" id="gridinfo-close-btn">${t('common.close') || 'Close'}</button>
        </div>
      </div>
    `;

    this.backdrop.querySelector('#gridinfo-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#gridinfo-close-btn')?.addEventListener('click', () => this.close());
    this.backdrop.addEventListener('click', (e) => {
      if (e.target === this.backdrop) this.close();
    });
  }
}
