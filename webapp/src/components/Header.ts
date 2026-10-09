import { t, getLanguage, toggleLanguage, onLanguageChange } from '../i18n';

export interface HeaderCallbacks {
  onOpenFeedback: () => void;
  onToggleMobileView: () => void;
}

export class Header {
  private element: HTMLElement;
  private callbacks: HeaderCallbacks;

  constructor(callbacks: HeaderCallbacks) {
    this.callbacks = callbacks;
    this.element = document.createElement('header');
    this.element.className = 'app-header';
    this.render();

    onLanguageChange(() => this.render());
  }

  public getElement(): HTMLElement {
    return this.element;
  }

  private render(): void {
    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    const currentLang = getLanguage();

    this.element.innerHTML = `
      <div class="header-brand">
        <div class="brand-icon">RO</div>
        <div class="brand-title">
          <span>${t('app.title')}</span>
          <span class="brand-badge">${t('app.version')}</span>
        </div>
      </div>
      <div class="header-actions">
        <button class="btn btn-secondary btn-feedback" id="header-feedback-btn">
          <span>💬</span>
          <span>${t('app.feedbackBtn')}</span>
        </button>
        <a href="/docs" target="_blank" rel="noopener noreferrer" class="btn btn-outline" title="${t('app.devDocsBtn')}">
          <span>⚡</span>
          <span>API</span>
        </a>
        <button class="btn-icon" id="header-lang-btn" title="${t('app.langToggle')}">
          <span style="font-size: 12px; font-weight: 700;">${currentLang.toUpperCase()}</span>
        </button>
        <button class="btn-icon" id="header-theme-btn" title="${t('app.themeToggle')}">
          <span>${isDark ? '☀️' : '🌙'}</span>
        </button>
        <button class="btn btn-secondary mobile-only-btn" id="header-mobile-view-btn" style="display: none;">
          <span>🗺️</span>
        </button>
      </div>
    `;

    // Event Listeners
    this.element.querySelector('#header-feedback-btn')?.addEventListener('click', () => {
      this.callbacks.onOpenFeedback();
    });

    this.element.querySelector('#header-lang-btn')?.addEventListener('click', () => {
      toggleLanguage();
    });

    this.element.querySelector('#header-theme-btn')?.addEventListener('click', () => {
      const current = document.documentElement.getAttribute('data-theme');
      const next = current === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      localStorage.setItem('pytransdat_theme', next);
      this.render();
    });

    this.element.querySelector('#header-mobile-view-btn')?.addEventListener('click', () => {
      this.callbacks.onToggleMobileView();
    });
  }
}
