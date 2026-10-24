import { t } from '../i18n';
import { APP_CONFIG } from '../config';
import { apiClient } from '../utils/apiClient';
import { showToast } from '../utils/toast';

export type FeedbackReportType = 'discrepancy' | 'grid_update' | 'bug' | 'feature' | 'api_down';

export class FeedbackPanel {
  private element: HTMLElement | null = null;
  private currentType: FeedbackReportType = 'discrepancy';
  private activeGridInfo: string = 'RO_TransDatRO_ETRS89_Stereo70.spg';
  private isMinimized: boolean = false;
  private isSubmitting: boolean = false;
  private renderTimestamp: number = 0;

  // Math challenge state
  private mathNum1: number = 0;
  private mathNum2: number = 0;
  private mathAnswer: number = 0;

  // Dragging state
  private isDragging: boolean = false;
  private dragStartX: number = 0;
  private dragStartY: number = 0;
  private initialLeft: number = 0;
  private initialTop: number = 0;

  constructor() {}

  public isOpen(): boolean {
    return this.element !== null;
  }

  public toggle(initialType?: FeedbackReportType): void {
    if (this.isOpen()) {
      if (initialType && this.currentType !== initialType) {
        this.setType(initialType);
      } else {
        this.close();
      }
    } else {
      this.open(initialType);
    }
  }

  public async open(initialType?: FeedbackReportType): Promise<void> {
    if (initialType) {
      this.currentType = initialType;
    }

    if (this.element) {
      this.element.style.display = 'flex';
      this.setType(this.currentType);
      return;
    }

    this.element = document.createElement('div');
    this.element.className = 'floating-feedback-panel';
    document.body.appendChild(this.element);

    this.render();
    await this.fetchGridInfo();
  }

  public close(): void {
    if (this.element) {
      this.element.remove();
      this.element = null;
    }
  }

  public setType(type: FeedbackReportType): void {
    this.currentType = type;
    if (this.element) {
      const select = this.element.querySelector('#fb-type-select') as HTMLSelectElement | null;
      if (select && select.value !== type) {
        select.value = type;
      }
      this.populateTemplateForCurrentType();
    }
  }

  private async fetchGridInfo(): Promise<void> {
    try {
      const info = await apiClient.getGridInfo();
      if (info && info.grid_file) {
        this.activeGridInfo = info.grid_file;
        // If current type is grid_update, refresh the template to display the loaded grid
        if (this.currentType === 'grid_update') {
          this.populateTemplateForCurrentType();
        }
      }
    } catch {
      // Fallback grid string already set
    }
  }

  private generateMathChallenge(): void {
    this.mathNum1 = Math.floor(Math.random() * 8) + 2; // 2..9
    this.mathNum2 = Math.floor(Math.random() * 8) + 1; // 1..8
    this.mathAnswer = this.mathNum1 + this.mathNum2;
  }

  private getTemplateText(type: FeedbackReportType): string {
    const rawTemplate = t(`feedback.templates.${type}`) || '';
    const nowIso = new Date().toISOString();
    const clientUrl = window.location.href;
    const browserEnv = `${navigator.userAgent} (Screen: ${window.screen.width}x${window.screen.height})`;

    return rawTemplate
      .replace('{currentGrid}', this.activeGridInfo)
      .replace('{timestamp}', nowIso)
      .replace('{clientUrl}', clientUrl)
      .replace('{browserEnv}', browserEnv);
  }

  private populateTemplateForCurrentType(): void {
    if (!this.element) return;
    const textarea = this.element.querySelector('#fb-message-input') as HTMLTextAreaElement | null;
    if (textarea) {
      textarea.value = this.getTemplateText(this.currentType);
    }
  }

  private toggleMinimize(): void {
    this.isMinimized = !this.isMinimized;
    if (this.element) {
      if (this.isMinimized) {
        this.element.classList.add('minimized');
      } else {
        this.element.classList.remove('minimized');
      }
      this.render();
    }
  }

  private copyToClipboard(): void {
    if (!this.element) return;
    const textarea = this.element.querySelector('#fb-message-input') as HTMLTextAreaElement | null;
    const text = textarea?.value || '';
    navigator.clipboard.writeText(text).then(() => {
      showToast(t('feedback.copiedToast') || 'Report copied to clipboard!', 'info');
    });
  }

  private openGithubIssue(): void {
    if (!this.element) return;
    const textarea = this.element.querySelector('#fb-message-input') as HTMLTextAreaElement | null;
    const text = textarea?.value || '';
    const title = encodeURIComponent(`[Feedback: ${this.currentType.toUpperCase()}]`);
    const body = encodeURIComponent(text);
    const url = `https://github.com/bogdanmorosanu/PyTransDatRO/issues/new?title=${title}&body=${body}`;
    window.open(url, '_blank', 'noopener,noreferrer');
  }

  private async handleSubmit(): Promise<void> {
    if (!this.element || this.isSubmitting) return;

    // 1. Minimum elapsed time check to stop automated headless scripts
    const elapsed = Date.now() - this.renderTimestamp;
    if (elapsed < 1200) {
      // Too fast for human typing/interaction
      showToast(t('feedback.mathError') || 'Please take a moment before submitting.', 'warning');
      return;
    }

    // 2. Validate human math check
    const mathInput = this.element.querySelector('#fb-math-input') as HTMLInputElement | null;
    const userAnswer = parseInt(mathInput?.value.trim() || '', 10);
    if (isNaN(userAnswer) || userAnswer !== this.mathAnswer) {
      showToast(t('feedback.mathError') || 'Incorrect human verification answer. Please try again.', 'warning');
      if (mathInput) {
        mathInput.focus();
        mathInput.select();
      }
      return;
    }

    // 3. Honeypot check (hidden input)
    const honeypot = (this.element.querySelector('#fb-hp-verify') as HTMLInputElement | null)?.value || '';

    // 4. Description content
    const textarea = this.element.querySelector('#fb-message-input') as HTMLTextAreaElement | null;
    const description = textarea?.value.trim() || '';

    if (description.length < 10) {
      showToast('Please provide a more detailed description before submitting.', 'warning');
      textarea?.focus();
      return;
    }

    this.isSubmitting = true;
    const submitBtn = this.element.querySelector('#fb-submit-btn') as HTMLButtonElement | null;
    if (submitBtn) {
      submitBtn.disabled = true;
      submitBtn.textContent = 'Submitting...';
    }

    try {
      const response = await apiClient.submitFeedback({
        report_type: this.currentType,
        description,
        client_timestamp: new Date().toISOString(),
        app_version: `v${APP_CONFIG.appVersion.replace(/^v/, '')}`,
        active_grid: this.activeGridInfo,
        honeypot,
      });

      showToast(response.message || t('feedback.successToast'), 'success');
      this.close();
    } catch (err) {
      console.warn('Feedback submit error:', err);
      showToast(t('feedback.errorToast') || 'Submission failed. You can copy the report or open a GitHub issue.', 'warning');
    } finally {
      this.isSubmitting = false;
      if (submitBtn) {
        submitBtn.disabled = false;
        submitBtn.textContent = t('feedback.submitBtn');
      }
    }
  }

  private render(): void {
    if (!this.element) return;

    this.renderTimestamp = Date.now();
    this.generateMathChallenge();

    if (this.isMinimized) {
      this.element.innerHTML = `
        <div class="drawer-header fgp-header" id="ffp-drag-handle">
          <div class="drawer-title-cluster fgp-title-cluster">
            <span class="drawer-title fgp-title">${t('feedback.title')}</span>
          </div>
          <div class="drawer-header-actions fgp-header-actions">
            <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="ffp-minimize-btn" title="${t('common.restore')}">□</button>
            <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="ffp-close-btn" title="${t('common.close')}">✕</button>
          </div>
        </div>
      `;
      this.attachHeaderListeners();
      this.attachDragListeners();
      return;
    }

    this.element.innerHTML = `
      <!-- Draggable Header -->
      <div class="drawer-header fgp-header" id="ffp-drag-handle">
        <div class="drawer-title-cluster fgp-title-cluster">
          <span class="drawer-title fgp-title">${t('feedback.title')}</span>
        </div>
        <div class="drawer-header-actions fgp-header-actions">
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="ffp-minimize-btn" title="${t('common.minimize')}">─</button>
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="ffp-close-btn" title="${t('common.close')}">✕</button>
        </div>
      </div>

      <!-- Main Body -->
      <div class="drawer-body ffp-body">
        <!-- Section 1: Feedback Content Card -->
        <div class="content-card">
          <div class="content-card-header">
            <div class="content-card-title">${t('feedback.title')}</div>
            <div class="content-card-subtitle">${t('feedback.desc')}</div>
          </div>

          <!-- Report Type Selector -->
          <div class="form-group ffp-field-group">
            <label class="form-label">${t('feedback.typeLabel')}</label>
            <select class="form-select form-select-sm" id="fb-type-select">
              <option value="discrepancy" ${this.currentType === 'discrepancy' ? 'selected' : ''}>${t('feedback.typeDiscrepancy')}</option>
              <option value="grid_update" ${this.currentType === 'grid_update' ? 'selected' : ''}>${t('feedback.typeGridUpdate')}</option>
              <option value="bug" ${this.currentType === 'bug' ? 'selected' : ''}>${t('feedback.typeBug')}</option>
              <option value="feature" ${this.currentType === 'feature' ? 'selected' : ''}>${t('feedback.typeFeature')}</option>
              <option value="api_down" ${this.currentType === 'api_down' ? 'selected' : ''}>${t('feedback.typeApiDown')}</option>
            </select>
          </div>

          <!-- Detailed Description / Template Area -->
          <div class="form-group ffp-message-group">
            <div class="callout-box">
              <span class="callout-title">${t('feedback.messageLabel')}</span>
              <p class="callout-body">${t('feedback.templateHint')}</p>
            </div>
            <textarea class="form-input mono ffp-textarea" id="fb-message-input" rows="7"></textarea>
          </div>

          <!-- Dedicated Command Actions Row (Aligned Right) -->
          <div class="ffp-actions-row">
            <button type="button" class="btn btn-outline ffp-aux-btn" id="fb-reset-template-btn" title="Revert description to default template for this report type">
              ${t('feedback.resetTemplateBtn')}
            </button>
            <button type="button" class="btn btn-outline ffp-aux-btn" id="fb-copy-btn" title="Copy report text to clipboard">
              ${t('feedback.copyBtn')}
            </button>
            <button type="button" class="btn btn-outline ffp-aux-btn" id="fb-github-btn" title="Open as GitHub issue in a new tab">
              ${t('feedback.openGithubBtn')}
            </button>
          </div>
        </div>

        <!-- Section 2: Human Verification Card -->
        <div class="content-card">
          <div class="content-card-header">
            <div class="content-card-title">${t('feedback.mathLabel')}</div>
            <div class="content-card-subtitle">${t('feedback.verificationSubtitle')}</div>
          </div>

          <!-- Zero-friction Bot Verification Row -->
          <div class="ffp-verify-content">
            <!-- Invisible Honeypot -->
            <input type="text" id="fb-hp-verify" name="website_hp" class="ffp-honeypot" tabindex="-1" autocomplete="off" />

            <!-- Lightweight Human Math Challenge -->
            <div class="ffp-math-cluster">
              <span class="ffp-math-expr mono">${this.mathNum1} + ${this.mathNum2}</span>
              <span class="ffp-math-equals">=</span>
              <input type="number" class="form-input mono ffp-math-input" id="fb-math-input" placeholder="${t('feedback.mathPlaceholder')}" min="0" max="99" />
            </div>
          </div>
        </div>
      </div>

      <!-- Footer -->
      <div class="ffp-footer">
        <button class="btn btn-secondary btn-sm" id="fb-cancel-btn">${t('common.close')}</button>
        <button class="btn btn-primary btn-sm" id="fb-submit-btn">${t('feedback.submitBtn')}</button>
      </div>

      <!-- Resizer Handle -->
      <div class="fgp-resize-handle" id="ffp-resize-handle" title="Resize window">
        <svg width="10" height="10" viewBox="0 0 10 10" fill="currentColor">
          <path d="M9 1v8H1l8-8zM9 5v4H5l4-4z"/>
        </svg>
      </div>
    `;

    this.populateTemplateForCurrentType();
    this.attachHeaderListeners();
    this.attachContentListeners();
    this.attachDragListeners();
    this.attachResizeListeners();
  }

  private attachHeaderListeners(): void {
    if (!this.element) return;
    this.element.querySelector('#ffp-minimize-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.toggleMinimize();
    });
    this.element.querySelector('#ffp-close-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.close();
    });
  }

  private attachContentListeners(): void {
    if (!this.element) return;

    const select = this.element.querySelector('#fb-type-select') as HTMLSelectElement | null;
    select?.addEventListener('change', () => {
      this.currentType = select.value as FeedbackReportType;
      this.populateTemplateForCurrentType();
    });

    this.element.querySelector('#fb-reset-template-btn')?.addEventListener('click', () => {
      this.populateTemplateForCurrentType();
    });

    this.element.querySelector('#fb-copy-btn')?.addEventListener('click', () => {
      this.copyToClipboard();
    });

    this.element.querySelector('#fb-github-btn')?.addEventListener('click', () => {
      this.openGithubIssue();
    });

    this.element.querySelector('#fb-cancel-btn')?.addEventListener('click', () => {
      this.close();
    });

    this.element.querySelector('#fb-submit-btn')?.addEventListener('click', () => {
      this.handleSubmit();
    });
  }

  private attachDragListeners(): void {
    if (!this.element) return;
    const handle = this.element.querySelector('#ffp-drag-handle') as HTMLElement | null;
    if (!handle) return;

    handle.addEventListener('pointerdown', (e: PointerEvent) => {
      if ((e.target as HTMLElement).closest('button')) return;
      if (window.innerWidth <= 768) return;

      this.isDragging = true;
      this.dragStartX = e.clientX;
      this.dragStartY = e.clientY;

      const rect = this.element!.getBoundingClientRect();
      this.initialLeft = rect.left;
      this.initialTop = rect.top;

      try {
        handle.setPointerCapture(e.pointerId);
      } catch {}

      const onPointerMove = (moveEvent: PointerEvent) => {
        if (!this.isDragging || !this.element) return;
        const dx = moveEvent.clientX - this.dragStartX;
        const dy = moveEvent.clientY - this.dragStartY;

        let newLeft = this.initialLeft + dx;
        let newTop = this.initialTop + dy;

        const maxLeft = Math.max(0, window.innerWidth - (this.element.offsetWidth || 300));
        const maxTop = Math.max(0, window.innerHeight - 50);

        newLeft = Math.max(0, Math.min(newLeft, maxLeft));
        newTop = Math.max(0, Math.min(newTop, maxTop));

        this.element.style.left = `${newLeft}px`;
        this.element.style.top = `${newTop}px`;
        this.element.style.right = 'auto';
        this.element.style.bottom = 'auto';
      };

      const onPointerUp = (upEvent: PointerEvent) => {
        this.isDragging = false;
        try {
          handle.releasePointerCapture(upEvent.pointerId);
        } catch {}
        handle.removeEventListener('pointermove', onPointerMove);
        handle.removeEventListener('pointerup', onPointerUp);
        handle.removeEventListener('pointercancel', onPointerUp);
      };

      handle.addEventListener('pointermove', onPointerMove);
      handle.addEventListener('pointerup', onPointerUp);
      handle.addEventListener('pointercancel', onPointerUp);
    });
  }

  private attachResizeListeners(): void {
    if (!this.element) return;
    const resizer = this.element.querySelector('#ffp-resize-handle') as HTMLElement | null;
    if (!resizer) return;

    resizer.addEventListener('pointerdown', (e: PointerEvent) => {
      e.preventDefault();
      e.stopPropagation();

      const startX = e.clientX;
      const startY = e.clientY;
      const startWidth = this.element!.offsetWidth;
      const startHeight = this.element!.offsetHeight;

      try {
        resizer.setPointerCapture(e.pointerId);
      } catch {}

      const onPointerMove = (moveEvent: PointerEvent) => {
        if (!this.element) return;
        moveEvent.preventDefault();
        moveEvent.stopPropagation();

        const dw = moveEvent.clientX - startX;
        const dh = moveEvent.clientY - startY;

        const minWidth = 360;
        const maxWidth = Math.max(minWidth, window.innerWidth - 40);
        const minHeight = 280;
        const maxHeight = Math.max(minHeight, window.innerHeight - 80);

        const newWidth = Math.max(minWidth, Math.min(startWidth + dw, maxWidth));
        const newHeight = Math.max(minHeight, Math.min(startHeight + dh, maxHeight));

        this.element.style.width = `${newWidth}px`;
        this.element.style.height = `${newHeight}px`;
        this.element.style.maxHeight = 'none';
      };

      const onPointerUp = (upEvent: PointerEvent) => {
        try {
          resizer.releasePointerCapture(upEvent.pointerId);
        } catch {}
        resizer.removeEventListener('pointermove', onPointerMove);
        resizer.removeEventListener('pointerup', onPointerUp);
        resizer.removeEventListener('pointercancel', onPointerUp);
      };

      resizer.addEventListener('pointermove', onPointerMove);
      resizer.addEventListener('pointerup', onPointerUp);
      resizer.addEventListener('pointercancel', onPointerUp);
    });
  }
}
