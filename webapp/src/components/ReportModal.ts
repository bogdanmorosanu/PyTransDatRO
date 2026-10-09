import { t } from '../i18n';
import type { BatchReportData } from './BatchConverter';

export class ReportModal {
  private backdrop: HTMLElement | null = null;

  public open(data: BatchReportData): void {
    const root = document.getElementById('modal-root');
    if (!root) return;

    this.close(); // Clean up if open

    this.backdrop = document.createElement('div');
    this.backdrop.className = 'modal-backdrop';

    this.backdrop.innerHTML = `
      <div class="modal-dialog">
        <div class="modal-header">
          <div class="modal-title">📊 ${t('report.title')}</div>
          <button class="btn-icon" id="modal-close-icon">✕</button>
        </div>
        <div class="modal-body">
          <div class="status-grid">
            <div class="status-badge-card">
              <div class="status-badge-num">${data.total}</div>
              <div class="status-badge-label">${t('report.totalPoints')}</div>
            </div>
            <div class="status-badge-card success">
              <div class="status-badge-num" style="color: var(--color-success);">${data.valid}</div>
              <div class="status-badge-label">${t('report.validPoints')}</div>
            </div>
            <div class="status-badge-card ${data.outOfGrid > 0 ? 'warning' : ''}">
              <div class="status-badge-num" style="color: ${data.outOfGrid > 0 ? 'var(--color-warning)' : 'var(--text-dim)'};">${data.outOfGrid}</div>
              <div class="status-badge-label">${t('report.outOfGrid')}</div>
            </div>
            <div class="status-badge-card ${data.invalid > 0 ? 'danger' : ''}">
              <div class="status-badge-num" style="color: ${data.invalid > 0 ? 'var(--color-danger)' : 'var(--text-dim)'};">${data.invalid}</div>
              <div class="status-badge-label">${t('report.invalidRows')}</div>
            </div>
          </div>

          ${
            data.outOfGrid > 0
              ? `
            <div style="padding: 10px 14px; background-color: var(--color-warning-bg); border: 1px solid var(--color-warning); border-radius: var(--radius-sm); font-size: 12px; color: var(--color-warning); margin-bottom: 14px;">
              ⚠ <strong>${data.outOfGrid}</strong> coordonate sunt în afara grilei oficiale de distorsiuni TransDatRO. Coordonatele acestora au fost propagate fără corecție spline.
            </div>
          `
              : ''
          }

          ${
            data.invalid > 0
              ? `
            <div style="padding: 10px 14px; background-color: var(--color-danger-bg); border: 1px solid var(--color-danger); border-radius: var(--radius-sm); font-size: 12px; color: var(--color-danger);">
              ✕ <strong>${data.invalid}</strong> linii conțin erori de sintaxă sau un număr incorect de coloane.
            </div>
          `
              : ''
          }
        </div>
        <div class="modal-footer">
          <button class="btn btn-primary" id="modal-close-btn">${t('report.closeBtn')}</button>
        </div>
      </div>
    `;

    root.appendChild(this.backdrop);

    this.backdrop.querySelector('#modal-close-icon')?.addEventListener('click', () => this.close());
    this.backdrop.querySelector('#modal-close-btn')?.addEventListener('click', () => this.close());
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
