import { t } from '../i18n';
import { APP_CONFIG } from '../config';
import { renderRadioGroup } from '../ui';

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
          <div class="modal-title">${t('volumeWarning.title')}</div>
          <button class="btn-icon modal-close-btn" id="vw-close-icon" title="${t('common.close') || 'Close'}">✕</button>
        </div>
        <div class="modal-body">
          <div class="callout-box fhp-callout-box accent-warning vw-callout">
            <span class="callout-title fhp-callout-title">${t('volumeWarning.calloutTitle')}</span>
            <p class="callout-body fhp-callout-body">
              ${t('volumeWarning.message', { count: count.toLocaleString() })}
            </p>
          </div>
          <label class="form-label vw-prompt-label">${t('volumeWarning.promptLabel')}</label>
          <div id="vw-radio-slot"></div>
        </div>
        <div class="modal-footer">
          <button class="btn btn-primary btn-sm" id="vw-confirm-btn">
            ${t('volumeWarning.btnConfirm')}
          </button>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    let selectedChoice: 'all' | 'sample' | 'skip' = 'all';

    const radioGroup = renderRadioGroup<'all' | 'sample' | 'skip'>({
      name: 'vw-rendering-choice',
      ariaLabel: t('volumeWarning.promptLabel'),
      selectedValue: 'all',
      options: [
        {
          value: 'all',
          id: 'vw-opt-all',
          label: t('volumeWarning.optAll'),
          description: t('volumeWarning.optAllDesc', { count: count.toLocaleString() }),
          badge: `${count.toLocaleString()} pts`,
        },
        {
          value: 'sample',
          id: 'vw-opt-sample',
          label: t('volumeWarning.optSample'),
          description: t('volumeWarning.optSampleDesc', {
            sample: APP_CONFIG.batchSampleSize.toLocaleString(),
          }),
          badge: `${APP_CONFIG.batchSampleSize.toLocaleString()} pts`,
        },
        {
          value: 'skip',
          id: 'vw-opt-skip',
          label: t('volumeWarning.optSkip'),
          description: t('volumeWarning.optSkipDesc'),
        },
      ],
      onChange: (val) => {
        selectedChoice = val;
      },
    });

    const slot = this.backdrop.querySelector('#vw-radio-slot');
    slot?.appendChild(radioGroup);

    const handleAction = (choice: 'all' | 'sample' | 'skip') => {
      this.close();
      onConfirm(choice);
    };

    this.backdrop.querySelector('#vw-close-icon')?.addEventListener('click', () => handleAction('skip'));
    this.backdrop.querySelector('#vw-confirm-btn')?.addEventListener('click', () => handleAction(selectedChoice));

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
