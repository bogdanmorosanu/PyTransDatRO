import { t, getLanguage, setLanguage, onLanguageChange } from '../i18n';
import { getLogoInlineSvg, syncFaviconWithTheme } from '../assets/logos';
import { renderDropdownMenu } from '../ui';
import { APP_CONFIG } from '../config';

export interface HeaderCallbacks {
  onNavigate: (page: 'studio' | 'telemetry') => void;
  onResetToFreshState?: () => void;
  onApiHealthChanged?: (online: boolean) => void;
  onOpenFeedback: (initialType?: 'discrepancy' | 'grid_update' | 'bug' | 'feature' | 'api_down') => void;
  onOpenGridInfo: () => void;
  onOpenGeodeticHelp: () => void;
  onToggleMobileView?: () => void;
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

    // Global listener for closing dropdowns on click outside or escape key
    document.addEventListener('click', (e) => {
      const dropdowns = this.element.querySelectorAll<HTMLDetailsElement>('details.header-dropdown-wrap');
      dropdowns.forEach((details) => {
        if (details.open && !details.contains(e.target as Node)) {
          details.open = false;
        }
      });
    });

    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        const dropdowns = this.element.querySelectorAll<HTMLDetailsElement>('details.header-dropdown-wrap');
        dropdowns.forEach((details) => {
          if (details.open) {
            details.open = false;
          }
        });
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
        this.callbacks.onApiHealthChanged?.(this.isApiOnline);
        this.render();
      }
    } catch {
      if (this.isApiOnline !== false) {
        this.isApiOnline = false;
        this.callbacks.onApiHealthChanged?.(this.isApiOnline);
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
    const dropdowns = this.element.querySelectorAll<HTMLDetailsElement>('details.header-dropdown-wrap');
    dropdowns.forEach((details) => {
      details.open = false;
    });
  }

  private render(): void {
    const isDark = document.documentElement.getAttribute('data-theme') !== 'light';
    const currentLang = getLanguage();
    const wasHelpOpen = (this.element.querySelector('#header-help-details') as HTMLDetailsElement | null)?.open ?? false;
    const wasLangOpen = (this.element.querySelector('#header-lang-details') as HTMLDetailsElement | null)?.open ?? false;

    this.element.innerHTML = `
      <!-- Left Cluster: Brand, Attribution & Primary App Navigation -->
      <div class="header-zone-left">
        <!-- Brand Cluster: Official C2 Map Projection Graticule Logo + Stacked Title/Version/Status -->
        <a href="/" class="header-brand-link" title="PyTransDatRO Home">
          <span class="brand-logo-wrap" id="header-brand-logo-wrap">${getLogoInlineSvg(36)}</span>
          <div class="brand-text-column">
            <span class="brand-title-text">${APP_CONFIG.appName}</span>
            <div class="brand-meta-subrow">
              <span class="brand-stage-badge">v${APP_CONFIG.appVersion.replace(/^v/, '')}</span>
              <span 
                class="health-dot ${this.isApiOnline ? 'online' : 'offline'}" 
                id="header-health-dot"
                title="${this.isApiOnline ? t('header.apiHealthyTooltip') : t('header.apiOfflineTooltip')}"
              ></span>
            </div>
          </div>
        </a>

        <!-- 2-Line Attribution: "powered by / PyTransDatRO" (Left Zone, no oblique ↗ arrow) -->
        <a 
          href="https://github.com/bogdanmorosanu/PyTransDatRO" 
          target="_blank" 
          rel="noopener noreferrer" 
          class="header-repo-link" 
          title="PyTransDatRO Open Source Repository"
        >
          <span class="repo-link-prefix">${t('header.poweredByLabel')}</span>
          <span class="repo-link-name">${t('header.poweredByBrand')}</span>
        </a>

        <!-- Subtle Vertical Separator between Brand & Navigation -->
        <div class="header-brand-nav-divider"></div>

        <!-- Primary View Switcher: Underline Tabs -->
        <nav class="header-nav-tabs" role="tablist">
          <button 
            type="button"
            class="header-tab-btn ${this.activePage === 'studio' ? 'active' : ''}" 
            id="nav-btn-studio" 
            role="tab" 
            aria-selected="${this.activePage === 'studio'}"
          >
            <span>${t('header.viewTransform')}</span>
            <span class="header-tab-indicator"></span>
          </button>
          <button 
            type="button"
            class="header-tab-btn ${this.activePage === 'telemetry' ? 'active' : ''}" 
            id="nav-btn-telemetry" 
            role="tab" 
            aria-selected="${this.activePage === 'telemetry'}"
          >
            <span>${t('header.viewUsageData')}</span>
            <span class="header-tab-indicator"></span>
          </button>
        </nav>
      </div>

      <!-- Right Cluster: Help Menu, Language Selector Menu & Theme Action -->
      <div class="header-zone-actions">
        <!-- Reusable Help Menu Slot -->
        <div id="header-help-slot"></div>

        <!-- Utilities: Language Selector Menu & Theme Action Button -->
        <div class="header-utilities-dock">
          <div id="header-lang-slot"></div>
          
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
        </div>
      </div>
    `;

    // Render Reusable Help Menu (Standard Menu Item Style)
    const helpDropdown = renderDropdownMenu({
      id: 'header-help-details',
      triggerLabel: t('header.helpMenu'),
      triggerVariant: 'menu-item',
      onOpenChange: (isOpen) => {
        if (isOpen) {
          const langDetails = this.element.querySelector('#header-lang-details') as HTMLDetailsElement | null;
          if (langDetails) langDetails.open = false;
        }
      },
      items: [
        { type: 'group', title: t('header.pytransdatroGroup') },
        {
          type: 'action',
          id: 'menu-help-app-btn',
          label: t('header.helpApp'),
          onClick: () => this.callbacks.onOpenGeodeticHelp(),
        },
        {
          type: 'action',
          id: 'menu-grid-info-btn',
          label: t('header.gridInfo'),
          onClick: () => this.callbacks.onOpenGridInfo(),
        },
        {
          type: 'action',
          id: 'menu-feedback-btn',
          label: t('header.feedback'),
          onClick: () => this.callbacks.onOpenFeedback('discrepancy'),
        },
        { type: 'divider', dashed: true },
        {
          type: 'link',
          label: t('header.apiDocs'),
          href: '/docs',
          isExternal: true,
        },
        { type: 'divider' },
        { type: 'group', title: t('header.officialSoftwareTitle') },
        {
          type: 'link',
          label: t('header.romposPortal'),
          href: 'https://rompos.ro/index.php/download',
          isExternal: true,
        },
        {
          type: 'link',
          label: t('header.romgeoGui'),
          href: 'https://rompos.ro/index.php/download/download/2-software/64-romgeo-table-convert-gui',
          isExternal: true,
        },
        {
          type: 'link',
          label: t('header.romgeoGrid'),
          href: 'https://rompos.ro/index.php/download/download/2-software/65-romgeo-grid',
          isExternal: true,
        },
        {
          type: 'link',
          label: t('header.romgeoWebApi'),
          href: 'https://api.romgeo.ro/api/v1/',
          isExternal: true,
        },
      ],
    });
    if (wasHelpOpen) {
      helpDropdown.open = true;
    }
    this.element.querySelector('#header-help-slot')?.replaceWith(helpDropdown);

    // Render Language Selector Dropdown Menu (Status-Displaying Menu Control)
    const langDropdown = renderDropdownMenu({
      id: 'header-lang-details',
      triggerLabel: currentLang.toUpperCase(),
      triggerVariant: 'menu-item',
      menuWidth: '152px',
      ariaLabel: t('app.langToggle'),
      onOpenChange: (isOpen) => {
        if (isOpen) {
          const helpDetails = this.element.querySelector('#header-help-details') as HTMLDetailsElement | null;
          if (helpDetails) helpDetails.open = false;
        }
      },
      items: [
        {
          type: 'action',
          id: 'menu-lang-ro-btn',
          label: t('header.langRo'),
          active: currentLang === 'ro',
          onClick: () => setLanguage('ro'),
        },
        {
          type: 'action',
          id: 'menu-lang-en-btn',
          label: t('header.langEn'),
          active: currentLang === 'en',
          onClick: () => setLanguage('en'),
        },
      ],
    });
    if (wasLangOpen) {
      langDropdown.open = true;
    }
    this.element.querySelector('#header-lang-slot')?.replaceWith(langDropdown);

    // Attach Event Listeners
    this.element.querySelector('#header-brand-logo-wrap')?.addEventListener('click', (e) => {
      e.preventDefault();
      e.stopPropagation();
      if (this.callbacks.onResetToFreshState) {
        this.callbacks.onResetToFreshState();
      } else {
        this.callbacks.onNavigate('studio');
      }
    });

    this.element.querySelector('.brand-title-text')?.addEventListener('click', (e) => {
      e.preventDefault();
      e.stopPropagation();
      this.callbacks.onNavigate('studio');
    });

    this.element.querySelector('.header-brand-link')?.addEventListener('click', (e) => {
      e.preventDefault();
      this.callbacks.onNavigate('studio');
    });

    this.element.querySelector('#nav-btn-studio')?.addEventListener('click', () => {
      this.callbacks.onNavigate('studio');
    });

    this.element.querySelector('#nav-btn-telemetry')?.addEventListener('click', () => {
      this.callbacks.onNavigate('telemetry');
    });

    this.element.querySelector('#header-theme-btn')?.addEventListener('click', () => {
      const current = document.documentElement.getAttribute('data-theme');
      const next = current === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      localStorage.setItem('pytransdat_theme', next);
      syncFaviconWithTheme(next);
      window.dispatchEvent(new CustomEvent('themechange', { detail: { theme: next } }));
      this.render();
    });
  }
}
