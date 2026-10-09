import { t } from '../i18n';
import { showToast } from '../utils/toast';

export class FeedbackModal {
  private backdrop: HTMLElement | null = null;

  public open(): void {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close();

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';

    this.backdrop.innerHTML = `
      <div class="modal-dialog">
        <div class="modal-header">
          <div class="modal-title">💬 ${t('feedback.title')}</div>
          <button class="btn-icon" id="fb-close-icon">✕</button>
        </div>
        <div class="modal-body">
          <p style="color: var(--text-muted); font-size: 13px; margin-bottom: 16px;">
            ${t('feedback.desc')}
          </p>

          <div class="form-group">
            <label class="form-label">${t('feedback.typeLabel')}</label>
            <select class="form-input" id="fb-type-select">
              <option value="discrepancy">${t('feedback.typeDiscrepancy')}</option>
              <option value="feature">${t('feedback.typeFeature')}</option>
            </select>
          </div>

          <div class="form-group">
            <label class="form-label">${t('feedback.coordsLabel')}</label>
            <input type="text" class="form-input mono" id="fb-coords-input" placeholder="e.g. N: 500000, E: 500000, H: 100" />
          </div>

          <div class="form-group">
            <label class="form-label">${t('feedback.messageLabel')}</label>
            <textarea class="form-input" id="fb-message-input" rows="4" placeholder="Descrie diferența observată sau formatul dorit..."></textarea>
          </div>

          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
            <div class="form-group">
              <label class="form-label">${t('feedback.nameLabel')}</label>
              <input type="text" class="form-input" id="fb-name-input" />
            </div>
            <div class="form-group">
              <label class="form-label">${t('feedback.emailLabel')}</label>
              <input type="email" class="form-input" id="fb-email-input" />
            </div>
          </div>
        </div>
        <div class="modal-footer">
          <button class="btn btn-secondary" id="fb-cancel-btn">Anulează</button>
          <button class="btn btn-primary" id="fb-submit-btn">${t('feedback.submitBtn')}</button>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    this.backdrop.querySelector('#fb-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#fb-cancel-btn')?.addEventListener('click', () => this.close());

    this.backdrop.querySelector('#fb-submit-btn')?.addEventListener('click', () => {
      showToast(t('feedback.successToast'), 'success');
      this.close();
    });

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
