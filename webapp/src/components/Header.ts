import { t, getLanguage, toggleLanguage, onLanguageChange } from '../i18n';
import logoSvg from '../assets/logo.svg';

export interface HeaderCallbacks {
  onNavigate: (page: 'studio' | 'telemetry') => void;
  onOpenFeedback: () => void;
  onOpenGridInfo: () => void;
  onOpenGeodeticHelp: () => void;
  onToggleMobileView: () => void;
}

export class Header {
  private element: HTMLElement;
  private callbacks: HeaderCallbacks;
  private activePage: 'studio' | 'telemetry' = 'studio';
  private isApiOnline: boolean = true;
  private healthCheckInterval: number | null = null;

  constructor(callbacks: HeaderCallbacks) {
    this.callbacks = callbacks;
    this.element = document.createElement('header');
    this.element.className = 'app-header';
    this.render();

    onLanguageChange(() => this.render());

    // Initial check and periodic background health poll (every 30 seconds)
    this.checkApiHealth();
    this.healthCheckInterval = window.setInterval(() => this.checkApiHealth(), 30000);

    // Global listener for closing dropdown on click outside or escape key
    document.addEventListener('click', (e) => {
      const details = this.element.querySelector('#header-help-details') as HTMLDetailsElement | null;
      if (details && details.open && !details.contains(e.target as Node)) {
        details.open = false;
      }
    });

    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        const details = this.element.querySelector('#header-help-details') as HTMLDetailsElement | null;
        if (details && details.open) {
          details.open = false;
        }
      }
    });
  }

  public getElement(): HTMLElement {
    return this.element;
  }

  public setActivePage(page: 'studio' | 'telemetry'): void {
    if (this.activePage !== page) {
      this.activePage = page;
      this.render();
    }
  }

  public destroy(): void {
    if (this.healthCheckInterval !== null) {
      clearInterval(this.healthCheckInterval);
      this.healthCheckInterval = null;
    }
  }

  private async checkApiHealth(): Promise<void> {
    try {
      const res = await fetch('/health', { method: 'GET', cache: 'no-cache' });
      const online = res.ok && res.status === 200;
      if (this.isApiOnline !== online) {
        this.isApiOnline = online;
        this.render();
      }
    } catch {
      if (this.isApiOnline !== false) {
        this.isApiOnline = false;
        this.render();
      }
    }
  }

  public toggleDropdown(): void {
    const details = this.element.querySelector('#header-help-details') as HTMLDetailsElement | null;
    if (details) {
      details.open = !details.open;
    }
  }

  public closeDropdown(): void {
    const details = this.element.querySelector('#header-help-details') as HTMLDetailsElement | null;
    if (details) {
      details.open = false;
    }
  }

  private render(): void {
    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    const currentLang = getLanguage();
    const wasHelpOpen = (this.element.querySelector('#header-help-details') as HTMLDetailsElement | null)?.open ?? false;

    this.element.innerHTML = `
      <!-- Left Cluster: Brand & Primary App Navigation -->
      <div class="header-zone-left">
        <!-- Brand Cluster: Logo + Stacked Title/Version/Status -->
        <a href="/" class="header-brand-link" title="PyTransDatRO Home">
          <img src="${logoSvg}" alt="PyTransDatRO Logo" class="brand-logo-img" width="36" height="36" />
          <div class="brand-text-column">
            <span class="brand-title-text">${t('app.title')}</span>
            <div class="brand-meta-subrow">
              <span class="brand-stage-badge">${t('app.version')}</span>
              <span 
                class="health-dot ${this.isApiOnline ? 'online' : 'offline'}" 
                title="${this.isApiOnline ? 'API Connected & Healthy' : 'API Offline - Click to report issue'}"
              ></span>
            </div>
          </div>
        </a>

        ${
          !this.isApiOnline
            ? `
          <button class="header-offline-alert" id="header-offline-alert-btn" title="${t('header.reportOffline')}">
            <svg width="13" height="13" viewBox="0 0 16 16" fill="currentColor">
              <path d="M8 1a7 7 0 100 14A7 7 0 008 1zm0 3a1 1 0 011 1v4a1 1 0 11-2 0V5a1 1 0 011-1zm0 8a1.25 1.25 0 110-2.5 1.25 1.25 0 010 2.5z"/>
            </svg>
            <span>${t('header.reportOffline')}</span>
          </button>
        `
            : ''
        }

        <!-- Subtle Vertical Separator between Brand & Navigation -->
        <div class="header-brand-nav-divider"></div>

        <!-- Primary View Switcher: GitHub-Style Underline Tabs -->
        <nav class="header-nav-tabs" role="tablist">
          <button 
            type="button"
            class="header-tab-btn ${this.activePage === 'studio' ? 'active' : ''}" 
            id="nav-btn-studio" 
            role="tab" 
            aria-selected="${this.activePage === 'studio'}"
          >
            <span>Transform</span>
            <span class="header-tab-indicator"></span>
          </button>
          <button 
            type="button"
            class="header-tab-btn ${this.activePage === 'telemetry' ? 'active' : ''}" 
            id="nav-btn-telemetry" 
            role="tab" 
            aria-selected="${this.activePage === 'telemetry'}"
          >
            <span>Usage</span>
            <span class="header-tab-indicator"></span>
          </button>
        </nav>
      </div>

      <!-- Right Cluster: Attribution, Resources, Feedback & Utilities -->
      <div class="header-zone-actions">
        <!-- by PyTransDatRO link -->
        <a 
          href="https://github.com/bogdanmorosanu/PyTransDatRO" 
          target="_blank" 
          rel="noopener noreferrer" 
          class="header-repo-link" 
          title="PyTransDatRO Open Source Repository"
        >
          <span>by PyTransDatRO</span>
          <svg width="10" height="10" viewBox="0 0 12 12" fill="none" class="ext-icon">
            <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
          </svg>
        </a>

        <!-- Help Dropdown using semantic <details><summary> for guaranteed browser expansion -->
        <details class="header-dropdown-wrap" id="header-help-details" ${wasHelpOpen ? 'open' : ''}>
          <summary class="btn btn-outline header-dropdown-trigger" id="header-help-trigger">
            <span>${t('header.helpMenu')}</span>
            <svg class="chevron-icon" width="12" height="12" viewBox="0 0 12 12" fill="none">
              <path d="M2.5 4.5L6 8L9.5 4.5" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
            </svg>
          </summary>

          <div class="header-dropdown-menu" id="header-help-dropdown">
            <div class="dropdown-group-title">${t('header.pytransdatroGroup')}</div>
            <button type="button" class="dropdown-item" id="menu-help-app-btn">
              <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">
                <path d="M5.5 5.5a2.5 2.5 0 0 1 4.5 1.5c0 1.5-2 2-2 3v.5M8 13.5h.01"/>
              </svg>
              <span>${t('header.helpApp')}</span>
            </button>
            <button type="button" class="dropdown-item" id="menu-grid-info-btn">
              <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5">
                <path d="M2 4h12M2 8h12M2 12h12M4 2v12M8 2v12M12 2v12"/>
              </svg>
              <span>${t('header.gridInfo')}</span>
            </button>

            <div class="dropdown-divider"></div>
            
            <div class="dropdown-group-title">${t('header.officialSoftwareTitle')}</div>
            <a href="https://rompos.ro/index.php/download" target="_blank" rel="noopener noreferrer" class="dropdown-item">
              <span>${t('header.romposPortal')}</span>
              <svg class="ext-icon" width="11" height="11" viewBox="0 0 12 12" fill="none">
                <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
              </svg>
            </a>
            <a href="https://rompos.ro/index.php/download/download/2-software/64-romgeo-table-convert-gui" target="_blank" rel="noopener noreferrer" class="dropdown-item">
              <span>${t('header.romgeoGui')}</span>
              <svg class="ext-icon" width="11" height="11" viewBox="0 0 12 12" fill="none">
                <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
              </svg>
            </a>
            <a href="https://rompos.ro/index.php/download/download/2-software/65-romgeo-grid" target="_blank" rel="noopener noreferrer" class="dropdown-item">
              <span>${t('header.romgeoGrid')}</span>
              <svg class="ext-icon" width="11" height="11" viewBox="0 0 12 12" fill="none">
                <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
              </svg>
            </a>
          </div>
        </details>

        <!-- API Link Button -->
        <a href="/docs" target="_blank" rel="noopener noreferrer" class="btn btn-outline header-btn-link" title="${t('app.devDocsBtn')}">
          <span>${t('header.apiDocs')}</span>
          <svg class="ext-icon" width="10" height="10" viewBox="0 0 12 12" fill="none">
            <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
          </svg>
        </a>

        <!-- Feedback Button (Consistent with outline action buttons) -->
        <button type="button" class="btn btn-outline" id="header-feedback-btn">
          <span>${t('header.feedback')}</span>
        </button>

        <div class="header-divider-v"></div>

        <!-- Utilities: Lang & Theme Toggles -->
        <div class="header-utilities-dock">
          <button type="button" class="btn-icon header-util-btn" id="header-lang-btn" title="${t('app.langToggle')}">
            <span style="font-size: 11px; font-weight: 700; font-family: 'JetBrains Mono', monospace;">${currentLang.toUpperCase()}</span>
          </button>
          
          <button type="button" class="btn-icon header-util-btn" id="header-theme-btn" title="${t('app.themeToggle')}">
            ${
              isDark
                ? `<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                     <circle cx="12" cy="12" r="5"></circle>
                     <line x1="12" y1="1" x2="12" y2="3"></line>
                     <line x1="12" y1="21" x2="12" y2="23"></line>
                     <line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line>
                     <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line>
                     <line x1="1" y1="12" x2="3" y2="12"></line>
                     <line x1="21" y1="12" x2="23" y2="12"></line>
                     <line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line>
                     <line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line>
                   </svg>`
                : `<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                     <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>
                   </svg>`
            }
          </button>

          <button type="button" class="btn btn-secondary mobile-only-btn" id="header-mobile-view-btn" style="display: none;">
            <span>${t('app.mobileViewToggle')}</span>
          </button>
        </div>
      </div>
    `;

    // Attach Event Listeners
    this.element.querySelector('#nav-btn-studio')?.addEventListener('click', () => {
      this.callbacks.onNavigate('studio');
    });

    this.element.querySelector('#nav-btn-telemetry')?.addEventListener('click', () => {
      this.callbacks.onNavigate('telemetry');
    });

    this.element.querySelector('#menu-help-app-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.closeDropdown();
      this.callbacks.onOpenGeodeticHelp();
    });

    this.element.querySelector('#menu-grid-info-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.closeDropdown();
      this.callbacks.onOpenGridInfo();
    });

    this.element.querySelectorAll('#header-help-dropdown a.dropdown-item').forEach((link) => {
      link.addEventListener('click', () => {
        this.closeDropdown();
      });
    });

    this.element.querySelector('#header-feedback-btn')?.addEventListener('click', () => {
      this.callbacks.onOpenFeedback();
    });

    this.element.querySelector('#header-offline-alert-btn')?.addEventListener('click', () => {
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
