import { t } from '../i18n';
import { APP_CONFIG } from '../config';

export class VolumeWarningModal {
  private backdrop: HTMLElement | null = null;
  private keydownHandler: ((e: KeyboardEvent) => void) | null = null;

  public open(count: number, onConfirm: (choice: 'all' | 'sample' | 'skip') => void): void {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';

    this.backdrop.innerHTML = `
      <div class="modal-dialog modal-dialog-compact">
        <div class="modal-header">
          <div class="modal-title warning">${t('volumeWarning.title')}</div>
          <button class="btn-icon modal-close-btn" id="vw-close-icon" title="${t('common.close') || 'Close'}">✕</button>
        </div>
        <div class="modal-body">
          <p class="modal-body-text">
            ${t('volumeWarning.message', { count: count.toLocaleString() })}
          </p>
        </div>
        <div class="modal-footer">
          <button class="btn btn-outline btn-sm" id="vw-skip-btn">
            ${t('volumeWarning.btnSkip')}
          </button>
          <button class="btn btn-secondary btn-sm" id="vw-sample-btn">
            ${t('volumeWarning.btnSample', { sample: APP_CONFIG.batchSampleSize.toLocaleString() })}
          </button>
          <button class="btn btn-primary btn-sm" id="vw-all-btn">
            ${t('volumeWarning.btnRenderAll')}
          </button>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    const handleAction = (choice: 'all' | 'sample' | 'skip') => {
      this.close();
      onConfirm(choice);
    };

    this.backdrop.querySelector('#vw-close-icon')?.addEventListener('click', () => handleAction('skip'));
    this.backdrop.querySelector('#vw-skip-btn')?.addEventListener('click', () => handleAction('skip'));
    this.backdrop.querySelector('#vw-sample-btn')?.addEventListener('click', () => handleAction('sample'));
    this.backdrop.querySelector('#vw-all-btn')?.addEventListener('click', () => handleAction('all'));

    this.backdrop.addEventListener('click', (e) => {
      if (e.target === this.backdrop) {
        handleAction('skip');
      }
    });

    this.keydownHandler = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        handleAction('skip');
      }
    };
    window.addEventListener('keydown', this.keydownHandler);
  }

  public close(): void {
    if (this.keydownHandler) {
      window.removeEventListener('keydown', this.keydownHandler);
      this.keydownHandler = null;
    }
    if (this.backdrop) {
      this.backdrop.remove();
      this.backdrop = null;
    }
  }
}
