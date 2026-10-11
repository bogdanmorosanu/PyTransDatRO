import { t } from '../i18n';

export class GeodeticHelpModal {
  private backdrop: HTMLElement | null = null;

  public open(): void {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';

    this.backdrop.innerHTML = `
      <div class="modal-dialog" style="max-width: 600px;">
        <div class="modal-header">
          <div class="modal-title">${t('geodeticHelp.title')}</div>
          <button class="btn-icon" id="gh-close-icon" aria-label="Close">✕</button>
        </div>
        <div class="modal-body" style="display: flex; flex-direction: column; gap: 16px; font-size: 13px; line-height: 1.6;">
          
          <!-- Datums & Reference Frames -->
          <div>
            <h4 style="font-size: 13px; font-weight: 700; color: var(--color-primary-light); margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.5px;">
              ${t('geodeticHelp.datumsTitle')}
            </h4>
            <div style="display: flex; flex-direction: column; gap: 8px; color: var(--text-main);">
              <div style="background-color: var(--bg-surface-elevated); padding: 8px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-subtle);">
                <strong>Stereo 70 (EPSG:31700):</strong>
                <p style="color: var(--text-muted); margin-top: 2px;">${t('geodeticHelp.stereo70Desc')}</p>
              </div>
              <div style="background-color: var(--bg-surface-elevated); padding: 8px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-subtle);">
                <strong>ETRS89 (EPSG:4258):</strong>
                <p style="color: var(--text-muted); margin-top: 2px;">${t('geodeticHelp.etrs89Desc')}</p>
              </div>
              <div style="background-color: var(--bg-surface-elevated); padding: 8px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-subtle);">
                <strong>Vertical Reference & Quasigeoid:</strong>
                <p style="color: var(--text-muted); margin-top: 2px;">${t('geodeticHelp.verticalDesc')}</p>
              </div>
            </div>
          </div>

          <!-- Coordinate Definitions -->
          <div>
            <h4 style="font-size: 13px; font-weight: 700; color: var(--color-primary-light); margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.5px;">
              ${t('geodeticHelp.coordsTitle')}
            </h4>
            <ul style="padding-left: 18px; color: var(--text-muted); display: flex; flex-direction: column; gap: 4px;">
              <li>${t('geodeticHelp.coordsPlanimetric')}</li>
              <li>${t('geodeticHelp.coordsGeographic')}</li>
            </ul>
          </div>

          <!-- Formats & Auto-Detection -->
          <div>
            <h4 style="font-size: 13px; font-weight: 700; color: var(--color-primary-light); margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.5px;">
              ${t('geodeticHelp.formatsTitle')}
            </h4>
            <p style="color: var(--text-muted);">
              ${t('geodeticHelp.formatsDesc')}
            </p>
          </div>

          <!-- Accuracy & Standards -->
          <div>
            <h4 style="font-size: 13px; font-weight: 700; color: var(--color-primary-light); margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.5px;">
              ${t('geodeticHelp.accuracyTitle')}
            </h4>
            <p style="color: var(--text-muted);">
              ${t('geodeticHelp.accuracyDesc')}
            </p>
          </div>

        </div>
        <div class="modal-footer">
          <button class="btn btn-secondary" id="gh-close-btn">${t('geodeticHelp.closeBtn')}</button>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    this.backdrop.querySelector('#gh-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#gh-close-btn')?.addEventListener('click', () => this.close());

    this.backdrop.addEventListener('click', (e) => {
      if (e.target === this.backdrop) this.close();
    });
  }

  public close(): void {
    if (this.backdrop) {
      this.backdrop.remove();
      this.backdrop = null;
    }
  }
}
