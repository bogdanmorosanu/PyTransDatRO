import { t } from '../i18n';
import { APP_CONFIG } from '../config';

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
          <div class="modal-title warning">${t('volumeWarning.title')}</div>
          <button class="btn-icon modal-close-btn" id="vw-close-icon" title="${t('common.close') || 'Close'}">✕</button>
        </div>
        <div class="modal-body">
          <p class="modal-body-text">
            ${t('volumeWarning.message', { count })}
          </p>

          <div class="flex-column-gap-md">
            <button class="btn btn-primary" id="vw-sample-btn">
              ${t('volumeWarning.btnSample', { sample: APP_CONFIG.batchSampleSize.toLocaleString() })}
            </button>
            <button class="btn btn-secondary" id="vw-all-btn">
              ${t('volumeWarning.btnRenderAll')}
            </button>
            <button class="btn btn-outline" id="vw-skip-btn">
              ${t('volumeWarning.btnSkip')}
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
