import { t } from '../i18n';

export class VolumeWarningModal {
  private backdrop: HTMLElement | null = null;

  public open(count: number, onConfirm: (choice: 'all' | 'sample' | 'skip') => void): void {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';

    this.backdrop.innerHTML = `
      <div class="modal-dialog">
        <div class="modal-header">
          <div class="modal-title" style="color: var(--color-warning);">⚠ ${t('volumeWarning.title')}</div>
          <button class="btn-icon" id="vw-close-icon">✕</button>
        </div>
        <div class="modal-body">
          <p style="font-size: 13px; line-height: 1.6; color: var(--text-main); margin-bottom: 20px;">
            ${t('volumeWarning.message', { count })}
          </p>

          <div style="display: flex; flex-direction: column; gap: 10px;">
            <button class="btn btn-primary" id="vw-sample-btn">
              ⚡ ${t('volumeWarning.btnSample')}
            </button>
            <button class="btn btn-secondary" id="vw-all-btn">
              ⚠ ${t('volumeWarning.btnRenderAll')}
            </button>
            <button class="btn btn-outline" id="vw-skip-btn">
              ✕ ${t('volumeWarning.btnSkip')}
            </button>
          </div>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    this.backdrop.querySelector('#vw-close-icon')?.addEventListener('click', () => {
      this.close();
      onConfirm('skip');
    });

    this.backdrop.querySelector('#vw-sample-btn')?.addEventListener('click', () => {
      this.close();
      onConfirm('sample');
    });

    this.backdrop.querySelector('#vw-all-btn')?.addEventListener('click', () => {
      this.close();
      onConfirm('all');
    });

    this.backdrop.querySelector('#vw-skip-btn')?.addEventListener('click', () => {
      this.close();
      onConfirm('skip');
    });
  }

  public close(): void {
    if (this.backdrop) {
      this.backdrop.remove();
      this.backdrop = null;
    }
  }
}
