import './style/variables.css';
import './style/base.css';
import './style/layout.css';
import './style/forms.css';
import './style/buttons.css';
import './style/badges.css';
import './style/menus.css';
import './style/panels.css';
import './style/telemetry.css';
import './style/map.css';

import { t, onLanguageChange } from './i18n';
import { syncFaviconWithTheme } from './assets/logos';
import { Header } from './components/Header';
import { PointConverter, type AngularUnitType } from './components/PointConverter';
import { MapView } from './components/MapView';
import { APP_CONFIG, E2E_ENABLED } from './config';
import { BatchConverter } from './components/BatchConverter';
import { FileUpload } from './components/FileUpload';
import { GeodeticHelpPanel } from './components/GeodeticHelpPanel';
import { GridMetadataPanel } from './components/GridMetadataPanel';
import { FeedbackPanel } from './components/FeedbackPanel';
import { TelemetryDashboard } from './components/TelemetryDashboard';
import { VolumeWarningModal } from './components/VolumeWarningModal';
import { apiClient, type TransformationOp } from './utils/apiClient';
import type { CoordinateOrder } from './utils/pointParser';
import { showToast } from './utils/toast';
import { renderSelect, renderCopyBadge } from './ui';

type ActiveTab = 'point' | 'batch' | 'file';
type AppPage = 'studio' | 'telemetry';

class App {
  private activeTab: ActiveTab = 'point';
  private currentOp: TransformationOp = 'ETRS89ToStereo70';
  private currentAngularUnit: AngularUnitType = 'degrees';
  private currentOrder: CoordinateOrder = 'NE';

  private rootElement: HTMLElement;
  private workbenchLayout!: HTMLElement;
  private telemetryDashboard!: TelemetryDashboard;
  private utilityPanel!: HTMLElement;
  private sharedConfigBar!: HTMLElement;
  private tabsBar!: HTMLElement;
  private tabContentContainer!: HTMLElement;
  private isWorkbenchCollapsed: boolean = false;
  private dockToggleBtn!: HTMLButtonElement;
  private bottomSheetState: 'peek' | 'half' | 'full' = 'half';
  private bottomSheetHandleBar!: HTMLElement;

  private header!: Header;
  private globalOfflineBanner!: HTMLElement;
  private isApiOnline: boolean = true;
  private pointConverter!: PointConverter;
  private batchConverter!: BatchConverter;
  private fileUpload!: FileUpload;
  private gridMetadataPanel = new GridMetadataPanel();
  private geodeticHelpPanel = new GeodeticHelpPanel();
  private mapView!: MapView;

  private feedbackPanel = new FeedbackPanel();
  private volumeWarningModal = new VolumeWarningModal();

  constructor() {
    document.title = APP_CONFIG.appName;
    this.rootElement = document.getElementById('app-root')!;
    this.readInitialQueryParams();
    this.initTheme();
    this.initLayout();
    this.mapView.initMap();
    if (E2E_ENABLED) {
      window.__e2e = {
        mapView: () => this.mapView.getSnapshot(),
        fireMapClick: (lat: number, lng: number) => this.mapView.fireMapClick(lat, lng),
      };
    }
    this.mapView.setAngularUnit(this.currentAngularUnit);

    onLanguageChange(() => {
      this.renderSharedConfigBar();
      this.renderTabs();
      if (!this.isApiOnline) {
        this.updateOfflineBanner(false);
      }
    });
  }

  private readInitialQueryParams(): void {
    const params = new URLSearchParams(window.location.search);
    if (params.has('op')) {
      const op = params.get('op') as TransformationOp;
      if (op === 'Stereo70ToETRS89' || op === 'ETRS89ToStereo70') {
        this.currentOp = op;
      }
    }
    if (params.has('unit')) {
      const unit = params.get('unit') as AngularUnitType;
      if (unit === 'degrees' || unit === 'dms' || unit === 'radians') {
        this.currentAngularUnit = unit;
      }
    }
    if (params.has('order')) {
      const order = params.get('order') as CoordinateOrder;
      if (order === 'NE' || order === 'EN') {
        this.currentOrder = order;
      }
    }
  }

  private initTheme(): void {
    const savedTheme = (localStorage.getItem('pytransdat_theme') || 'dark') as 'dark' | 'light';
    document.documentElement.setAttribute('data-theme', savedTheme);
    document.documentElement.removeAttribute('data-btn-color');
    syncFaviconWithTheme(savedTheme === 'light' ? 'light' : 'dark');
  }

  private initLayout(): void {
    // 1. Header
    this.header = new Header({
      onNavigate: (page) => this.navigateToPage(page),
      onResetToFreshState: () => this.resetToFreshState(),
      onApiHealthChanged: (online) => this.updateOfflineBanner(online),
      onOpenFeedback: (initialType) => this.feedbackPanel.toggle(initialType),
      onOpenGridInfo: () => this.gridMetadataPanel.toggle(),
      onOpenGeodeticHelp: () => this.geodeticHelpPanel.toggle(),
    });
    this.rootElement.appendChild(this.header.getElement());

    // 1b. Global Offline Alert Banner (directly beneath header)
    this.globalOfflineBanner = document.createElement('div');
    this.globalOfflineBanner.id = 'global-offline-banner';
    this.globalOfflineBanner.className = 'global-offline-banner';
    this.globalOfflineBanner.style.display = 'none';
    this.rootElement.appendChild(this.globalOfflineBanner);

    // 2. Workbench Layout
    this.workbenchLayout = document.createElement('main');
    this.workbenchLayout.className = 'workbench-layout';
    const workbench = this.workbenchLayout;

    // 2b. Telemetry Dashboard (Full View)
    this.telemetryDashboard = new TelemetryDashboard({
      onBackToStudio: () => this.navigateToPage('studio'),
    });
    this.rootElement.appendChild(this.telemetryDashboard.getElement());
    this.telemetryDashboard.getElement().style.display = 'none';

    // Left Panel
    this.utilityPanel = document.createElement('div');
    this.utilityPanel.className = 'utility-panel sheet-half';

    // Shared Configuration Bar
    this.sharedConfigBar = document.createElement('div');
    this.sharedConfigBar.className = 'shared-config-bar';
    this.utilityPanel.appendChild(this.sharedConfigBar);
    this.renderSharedConfigBar();

    // Folder Bookmark Tabs Bar
    this.tabsBar = document.createElement('nav');
    this.tabsBar.className = 'tabs-bar';
    this.renderTabs();
    this.utilityPanel.appendChild(this.tabsBar);

    // Tab Content Area
    this.tabContentContainer = document.createElement('div');
    this.tabContentContainer.style.flex = '1';
    this.tabContentContainer.style.display = 'flex';
    this.tabContentContainer.style.flexDirection = 'column';
    this.utilityPanel.appendChild(this.tabContentContainer);

    workbench.appendChild(this.utilityPanel);

    // Right Panel: Map
    this.mapView = new MapView({
      onMapPointClicked: async (lat, lon) => {
        // Non-disruptive map click: only capture when Point tab is actively selected
        if (this.activeTab !== 'point') {
          return;
        }

        try {
          const res = await apiClient.transformPoint({
            op: 'ETRS89ToStereo70',
            coos: [lat, lon],
            unit: 'degrees',
          });

          const isUntransformable = res.warning === 'Out of grid' || res.warning === 'No data on grid';

          if (isUntransformable) {
            this.pointConverter.clearCoordinates();
            this.pointConverter.setQualityStatus(res.warning === 'No data on grid' ? 'no_data' : 'out_of_grid');
            showToast(res.warning === 'No data on grid' ? t('quality.noGridData') : t('point.outOfGridWarn'), 'warning');
            return;
          }

          if (this.currentOp === 'ETRS89ToStereo70') {
            this.pointConverter.setCoordinatesFromMap(lat, lon);
          } else {
            // Stereo 70 ➔ ETRS89
            if (res.coos && res.coos.length >= 2) {
              this.pointConverter.setStereoCoordinatesFromMap(res.coos[0], res.coos[1]);
            }
          }
        } catch (err) {
          console.error('Failed to calculate coordinates from map click', err);
          this.pointConverter.clearCoordinates();
        }
      },
      onOpenVolumeWarning: (count, callback) => {
        this.volumeWarningModal.open(count, callback);
      },
    });

    workbench.appendChild(this.mapView.getElement());
    this.rootElement.appendChild(workbench);

    // 3. Child Components with 3-Layer Map Decoupled Centering & Zoom to Extents
    this.pointConverter = new PointConverter({
      onPointTransformed: (event) => {
        this.mapView.updateSinglePoint(event);
      },
      onCenterOnMap: ({ lat, lon, outOfGrid }) => {
        if (!outOfGrid) {
          const ok = this.mapView.centerOnPoint(lat, lon);
          if (!ok) {
            showToast(t('point.outOfGridWarn'), 'warning');
          }
        } else {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
      onZoomToExtents: ({ lat, lon, outOfGrid }) => {
        if (!outOfGrid) {
          const ok = this.mapView.zoomToPoint(lat, lon, 14);
          if (!ok) {
            showToast(t('point.outOfGridWarn'), 'warning');
          }
        } else {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
    });

    this.batchConverter = new BatchConverter({
      onPointsPlotted: (points) => {
        this.mapView.plotBatchPoints(points);
      },
      onCenterOnMap: (points) => {
        const ok = this.mapView.centerOnBatch(points);
        if (!ok) {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
      onZoomToExtents: (points) => {
        const ok = this.mapView.zoomToTextExtents(points);
        if (!ok) {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
    });

    this.fileUpload = new FileUpload({
      onPointsPlotted: (points) => {
        this.mapView.plotFilePoints(points);
      },
      onCenterOnMap: (points) => {
        const ok = this.mapView.centerOnFiles(points);
        if (!ok) {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
      onZoomToExtents: (points) => {
        const ok = this.mapView.zoomToFilesExtents(points);
        if (!ok) {
          showToast(t('point.outOfGridWarn'), 'warning');
        }
      },
      onClearFiles: () => {
        this.mapView.clearFilePoints();
      },
    });

    this.pointConverter.setDirection(this.currentOp);
    this.pointConverter.setAngularUnit(this.currentAngularUnit);
    this.pointConverter.setCoordinateOrder(this.currentOrder);
    this.batchConverter.setDirection(this.currentOp);
    this.batchConverter.setAngularUnit(this.currentAngularUnit);
    this.batchConverter.setCoordinateOrder(this.currentOrder);
    this.fileUpload.setDirection(this.currentOp);
    this.fileUpload.setAngularUnit(this.currentAngularUnit);
    this.fileUpload.setCoordinateOrder(this.currentOrder);
    this.mapView.setCoordinateOrder(this.currentOrder);

    this.updateActiveTabContent();
    this.initDockableWorkbench();
    this.initBottomSheet();
  }

  private setDirection(op: TransformationOp): void {
    if (this.currentOp === op) return;
    this.currentOp = op;
    const isEtrs = op === 'ETRS89ToStereo70';
    const slider = this.sharedConfigBar.querySelector('.sliding-toggle[data-direction]');
    if (slider) {
      slider.setAttribute('data-direction', isEtrs ? 'etrs' : 's70');
      slider.querySelectorAll('.sliding-toggle-btn').forEach((btn) => {
        const btnOp = btn.getAttribute('data-op');
        const active = btnOp === op;
        btn.classList.toggle('active', active);
        btn.setAttribute('aria-checked', String(active));
      });
    }
    this.pointConverter.setDirection(op);
    this.batchConverter.setDirection(op);
    this.fileUpload.setDirection(op);
  }

  private setAngularUnit(unit: AngularUnitType): void {
    if (this.currentAngularUnit === unit) return;
    this.currentAngularUnit = unit;
    this.updateFormatExampleBadge();
    this.pointConverter.setAngularUnit(unit);
    this.batchConverter.setAngularUnit(unit);
    this.fileUpload.setAngularUnit(unit);
    this.mapView.setAngularUnit(unit);
  }

  private updateFormatExampleBadge(): void {
    let sampleText = '45.84123450';
    if (this.currentAngularUnit === 'dms') {
      sampleText = '45°50\'28.44420"';
    } else if (this.currentAngularUnit === 'radians') {
      sampleText = '0.7997813210';
    }
    const currentBadge = this.sharedConfigBar.querySelector('#angle-format-example-badge');
    const newBadge = renderCopyBadge({
      id: 'angle-format-example-badge',
      text: sampleText,
      tooltip: t('point.clickToSelect'),
      onCopied: () => showToast(t('point.copiedToast'), 'success'),
    });
    currentBadge?.replaceWith(newBadge);
  }

  private setCoordinateOrder(order: CoordinateOrder): void {
    if (this.currentOrder === order) return;
    this.currentOrder = order;
    const orderSlider = this.sharedConfigBar.querySelector('#shared-order-control');
    if (orderSlider) {
      orderSlider.setAttribute('data-order', order);
      orderSlider.querySelectorAll('.sliding-toggle-btn').forEach((btn) => {
        const btnOrder = btn.getAttribute('data-order-val');
        const active = btnOrder === order;
        btn.classList.toggle('active', active);
        btn.setAttribute('aria-checked', String(active));
      });
    }
    this.pointConverter.setCoordinateOrder(order);
    this.batchConverter.setCoordinateOrder(order);
    this.fileUpload.setCoordinateOrder(order);
    this.mapView.setCoordinateOrder(order);
  }

  private renderSharedConfigBar(): void {
    const isEtrs = this.currentOp === 'ETRS89ToStereo70';
    let sampleText = '45.84123450';
    if (this.currentAngularUnit === 'dms') {
      sampleText = '45°50\'28.44420"';
    } else if (this.currentAngularUnit === 'radians') {
      sampleText = '0.7997813210';
    }

    this.sharedConfigBar.innerHTML = `
      <!-- Discrete Section Landmark (Archetype 5) -->
      <div class="section-micro-label">${t('config.sectionTitle')}</div>

      <!-- Row 1: Sliding Direction Toggle (ETRS89 to Stereo70, Stereo70 to ETRS89) -->
      <div class="config-row">
        <label class="form-label">${t('point.directionLabel')}</label>
        <div class="sliding-toggle" data-direction="${isEtrs ? 'etrs' : 's70'}" role="radiogroup" aria-label="${t('point.directionLabel')}">
          <div class="sliding-toggle-glider"></div>
          <button type="button" class="sliding-toggle-btn ${isEtrs ? 'active' : ''}" data-op="ETRS89ToStereo70" role="radio" aria-checked="${isEtrs}">
            ${t('point.opETRSToS70')}
          </button>
          <button type="button" class="sliding-toggle-btn ${!isEtrs ? 'active' : ''}" data-op="Stereo70ToETRS89" role="radio" aria-checked="${!isEtrs}">
            ${t('point.opS70ToETRS')}
          </button>
        </div>
      </div>

      <!-- Row 2: Coordinate Order Sliding Toggle -->
      <div class="config-row">
        <label class="form-label">${t('common.coordOrder')}</label>
        <div class="sliding-toggle" id="shared-order-control" data-order="${this.currentOrder}" role="radiogroup" aria-label="${t('common.coordOrder')}">
          <div class="sliding-toggle-glider"></div>
          <button type="button" class="sliding-toggle-btn ${this.currentOrder === 'NE' ? 'active' : ''}" data-order-val="NE" role="radio" aria-checked="${this.currentOrder === 'NE'}">
            NE(H) / LatLon(h)
          </button>
          <button type="button" class="sliding-toggle-btn ${this.currentOrder === 'EN' ? 'active' : ''}" data-order-val="EN" role="radio" aria-checked="${this.currentOrder === 'EN'}">
            EN(H) / LonLat(h)
          </button>
        </div>
      </div>

      <!-- Row 3: Angular Unit + Example -->
      <div class="config-row">
        <label class="form-label">${t('point.unitLabel')}</label>
        <div class="config-unit-cluster">
          <div id="angular-unit-slot"></div>
          <div class="config-example-subgroup">
            <span class="info-label">${t('point.exampleLabel')}</span>
            <div id="format-example-slot"></div>
          </div>
        </div>
      </div>
    `;

    // Direction Toggle Listeners
    this.sharedConfigBar.querySelectorAll('.sliding-toggle[data-direction] .sliding-toggle-btn').forEach((btn) => {
      btn.addEventListener('click', () => {
        const op = btn.getAttribute('data-op') as TransformationOp;
        if (op) this.setDirection(op);
      });
    });

    // Row 3: Reusable renderSelect for Angular Unit
    const select = renderSelect<AngularUnitType>({
      id: 'angular-unit-select',
      ariaLabel: t('point.unitLabel'),
      className: 'config-unit-select',
      selectedValue: this.currentAngularUnit,
      options: [
        { value: 'degrees', label: t('point.unitDeg') },
        { value: 'dms', label: t('point.unitDMS') },
        { value: 'radians', label: t('point.unitRad') },
      ],
      onChange: (unit) => this.setAngularUnit(unit),
    });
    this.sharedConfigBar.querySelector('#angular-unit-slot')?.replaceWith(select);

    // Tier 2: Coordinate Order Sliding Toggle Listeners
    this.sharedConfigBar.querySelectorAll('#shared-order-control .sliding-toggle-btn').forEach((btn) => {
      btn.addEventListener('click', () => {
        const order = btn.getAttribute('data-order-val') as CoordinateOrder;
        if (order) this.setCoordinateOrder(order);
      });
    });

    // Tier 3: Reusable renderCopyBadge
    const copyBadge = renderCopyBadge({
      id: 'angle-format-example-badge',
      text: sampleText,
      tooltip: t('point.clickToSelect'),
      onCopied: () => showToast(t('point.copiedToast'), 'success'),
    });
    this.sharedConfigBar.querySelector('#format-example-slot')?.replaceWith(copyBadge);
  }

  private renderTabs(): void {
    this.tabsBar.innerHTML = `
      <button class="tab-btn ${this.activeTab === 'point' ? 'active' : ''}" id="tab-point" data-tab="point">
        <span>${t('tabs.point')}</span>
        <span class="tab-indicator"></span>
      </button>
      <button class="tab-btn ${this.activeTab === 'batch' ? 'active' : ''}" id="tab-text" data-tab="batch">
        <span>${t('tabs.batch')}</span>
        <span class="tab-indicator"></span>
      </button>
      <button class="tab-btn ${this.activeTab === 'file' ? 'active' : ''}" id="tab-files" data-tab="file">
        <span>${t('tabs.file')}</span>
        <span class="tab-indicator"></span>
      </button>
    `;

    this.tabsBar.querySelectorAll('.tab-btn').forEach((btn) => {
      btn.addEventListener('click', () => {
        const tab = btn.getAttribute('data-tab') as ActiveTab;
        if (tab) this.switchTab(tab);
      });
    });
  }

  private switchTab(tab: ActiveTab): void {
    if (this.activeTab === tab) return;
    this.activeTab = tab;
    this.renderTabs();
    this.updateActiveTabContent();
  }

  private updateActiveTabContent(): void {
    this.tabContentContainer.innerHTML = '';
    switch (this.activeTab) {
      case 'point':
        this.tabContentContainer.appendChild(this.pointConverter.getElement());
        break;
      case 'batch':
        this.tabContentContainer.appendChild(this.batchConverter.getElement());
        break;
      case 'file':
        this.tabContentContainer.appendChild(this.fileUpload.getElement());
        break;
    }
  }

  private initDockableWorkbench(): void {
    const saved = localStorage.getItem('pytransdat_workbench_collapsed') === 'true';
    if (saved && window.innerWidth > 768) {
      this.isWorkbenchCollapsed = true;
      this.utilityPanel.classList.add('is-collapsed');
      // Invalidate Leaflet map size on initial load so map fills 100% of container
      requestAnimationFrame(() => this.mapView.invalidateSize());
    }

    this.dockToggleBtn = document.createElement('button');
    this.dockToggleBtn.type = 'button';
    this.dockToggleBtn.className = `workbench-dock-toggle ${this.isWorkbenchCollapsed ? 'is-collapsed' : ''}`;
    this.dockToggleBtn.setAttribute('aria-label', this.isWorkbenchCollapsed ? 'Expand workbench dock' : 'Collapse workbench dock');
    this.dockToggleBtn.setAttribute('aria-expanded', String(!this.isWorkbenchCollapsed));
    this.dockToggleBtn.title = this.isWorkbenchCollapsed ? 'Expand workbench (Ctrl+\\ or Alt+W)' : 'Collapse workbench (Ctrl+\\ or Alt+W)';
    this.dockToggleBtn.innerHTML = `
      <span class="workbench-dock-toggle-icon">${this.isWorkbenchCollapsed ? '▶' : '◀'}</span>
    `;
    this.workbenchLayout.appendChild(this.dockToggleBtn);

    this.dockToggleBtn.addEventListener('click', () => {
      this.toggleWorkbenchDock();
    });

    this.utilityPanel.addEventListener('transitionend', (e) => {
      if (e.target === this.utilityPanel && (e.propertyName === 'margin-left' || e.propertyName === 'transform')) {
        this.mapView.invalidateSize();
      }
    });

    window.addEventListener('keydown', (e: KeyboardEvent) => {
      if ((e.ctrlKey && e.key === '\\') || (e.altKey && (e.key === 'w' || e.key === 'W'))) {
        e.preventDefault();
        this.toggleWorkbenchDock();
      }
    });
  }

  private toggleWorkbenchDock(): void {
    if (window.innerWidth <= 768) return;
    this.isWorkbenchCollapsed = !this.isWorkbenchCollapsed;
    this.utilityPanel.classList.toggle('is-collapsed', this.isWorkbenchCollapsed);
    this.dockToggleBtn.classList.toggle('is-collapsed', this.isWorkbenchCollapsed);
    this.dockToggleBtn.setAttribute('aria-label', this.isWorkbenchCollapsed ? 'Expand workbench dock' : 'Collapse workbench dock');
    this.dockToggleBtn.setAttribute('aria-expanded', String(!this.isWorkbenchCollapsed));
    this.dockToggleBtn.title = this.isWorkbenchCollapsed ? 'Expand workbench (Ctrl+\\ or Alt+W)' : 'Collapse workbench (Ctrl+\\ or Alt+W)';
    const icon = this.dockToggleBtn.querySelector('.workbench-dock-toggle-icon');
    if (icon) {
      icon.textContent = this.isWorkbenchCollapsed ? '▶' : '◀';
    }
    localStorage.setItem('pytransdat_workbench_collapsed', String(this.isWorkbenchCollapsed));
    this.mapView.invalidateSize();
  }

  private initBottomSheet(): void {
    this.bottomSheetHandleBar = document.createElement('div');
    this.bottomSheetHandleBar.className = 'bottom-sheet-handle-bar';
    this.bottomSheetHandleBar.setAttribute('role', 'button');
    this.bottomSheetHandleBar.setAttribute('aria-label', 'Drag or tap bottom sheet');
    this.bottomSheetHandleBar.innerHTML = `<div class="bottom-sheet-handle-pill"></div>`;
    this.utilityPanel.insertBefore(this.bottomSheetHandleBar, this.utilityPanel.firstChild);

    this.setBottomSheetState('half');

    let touchStartY = 0;
    let touchCurrentY = 0;
    let isTouching = false;
    let didDrag = false;

    // Tap handle to cycle states: half -> full -> peek -> half
    this.bottomSheetHandleBar.addEventListener('click', () => {
      if (window.innerWidth > 768) return;
      if (didDrag) {
        didDrag = false;
        return;
      }
      if (this.bottomSheetState === 'half') {
        this.setBottomSheetState('full');
      } else if (this.bottomSheetState === 'full') {
        this.setBottomSheetState('peek');
      } else {
        this.setBottomSheetState('half');
      }
    });

    // Touch swipe gestures
    this.bottomSheetHandleBar.addEventListener(
      'touchstart',
      (e: TouchEvent) => {
        if (window.innerWidth > 768) return;
        if (e.touches.length === 1) {
          isTouching = true;
          didDrag = false;
          touchStartY = e.touches[0].clientY;
          touchCurrentY = touchStartY;
          this.utilityPanel.style.transition = 'none';
        }
      },
      { passive: true },
    );

    this.bottomSheetHandleBar.addEventListener(
      'touchmove',
      (e: TouchEvent) => {
        if (!isTouching || e.touches.length !== 1) return;
        touchCurrentY = e.touches[0].clientY;
        const deltaY = touchCurrentY - touchStartY;
        if (Math.abs(deltaY) > 6) {
          didDrag = true;
        }

        const sheetHeight = this.utilityPanel.offsetHeight || window.innerHeight * 0.85;
        let baseTranslateY = 0;
        if (this.bottomSheetState === 'half') baseTranslateY = sheetHeight * 0.45;
        else if (this.bottomSheetState === 'peek') baseTranslateY = sheetHeight - 72;
        else baseTranslateY = 0;

        const currentTranslate = Math.max(0, Math.min(sheetHeight - 72, baseTranslateY + deltaY));
        this.utilityPanel.style.transform = `translateY(${currentTranslate}px)`;
      },
      { passive: true },
    );

    const finishTouch = () => {
      if (!isTouching) return;
      isTouching = false;

      const deltaY = touchCurrentY - touchStartY;
      const threshold = 40;
      const sheetHeight = this.utilityPanel.offsetHeight || window.innerHeight * 0.85;

      let targetState = this.bottomSheetState;

      if (this.bottomSheetState === 'half') {
        if (deltaY < -threshold) {
          targetState = 'full';
        } else if (deltaY > threshold) {
          targetState = 'peek';
        }
      } else if (this.bottomSheetState === 'peek') {
        if (deltaY < -sheetHeight * 0.45) {
          targetState = 'full';
        } else if (deltaY < -threshold) {
          targetState = 'half';
        }
      } else if (this.bottomSheetState === 'full') {
        if (deltaY > sheetHeight * 0.45) {
          targetState = 'peek';
        } else if (deltaY > threshold) {
          targetState = 'half';
        }
      }

      this.utilityPanel.style.transition = '';
      this.utilityPanel.style.transform = '';
      this.setBottomSheetState(targetState);

      // Invalidate Leaflet view on gesture finish
      this.mapView?.invalidateSize();

      // Suppress synthesized click following drag
      if (didDrag) {
        setTimeout(() => {
          didDrag = false;
        }, 200);
      }
    };

    this.bottomSheetHandleBar.addEventListener('touchend', finishTouch);
    this.bottomSheetHandleBar.addEventListener('touchcancel', finishTouch);

    // Auto-expand to sheet-full when user focuses any input or textarea on mobile
    this.utilityPanel.addEventListener('focusin', (e: FocusEvent) => {
      if (window.innerWidth <= 768 && this.bottomSheetState !== 'full') {
        const target = e.target as HTMLElement | null;
        if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || target.tagName === 'SELECT')) {
          this.setBottomSheetState('full');
        }
      }
    });
  }

  private setBottomSheetState(state: 'peek' | 'half' | 'full'): void {
    this.bottomSheetState = state;
    this.utilityPanel.classList.remove('sheet-peek', 'sheet-half', 'sheet-full');
    this.utilityPanel.classList.add(`sheet-${state}`);
    this.mapView?.invalidateSize();
  }

  private updateOfflineBanner(online: boolean): void {
    this.isApiOnline = online;
    if (online) {
      this.globalOfflineBanner.style.display = 'none';
      this.globalOfflineBanner.innerHTML = '';
    } else {
      this.globalOfflineBanner.style.display = 'flex';
      this.globalOfflineBanner.innerHTML = `
        <div class="global-offline-banner-content">
          <svg width="15" height="15" viewBox="0 0 16 16" fill="currentColor" class="banner-warn-icon">
            <path d="M8 1a7 7 0 100 14A7 7 0 008 1zm0 3a1 1 0 011 1v4a1 1 0 11-2 0V5a1 1 0 011-1zm0 8a1.25 1.25 0 110-2.5 1.25 1.25 0 010 2.5z"/>
          </svg>
          <span class="banner-warn-text">${t('header.apiOfflineMessage')}</span>
          <button type="button" class="banner-report-btn" id="global-offline-report-btn">
            ${t('header.reportIssueLink')} ↗
          </button>
        </div>
      `;
      this.globalOfflineBanner.querySelector('#global-offline-report-btn')?.addEventListener('click', () => {
        this.feedbackPanel.toggle('api_down');
      });
    }
    setTimeout(() => this.mapView?.invalidateSize(), 60);
  }

  public resetToFreshState(): void {
    this.pointConverter.clearCoordinates();
    this.batchConverter.clear();
    this.fileUpload.clearQueue();

    this.mapView.clearSinglePoint();
    this.mapView.clearBatchPoints();
    this.mapView.clearFilePoints();
    this.mapView.toggleMeasureMode(false);
    this.mapView.zoomToRomaniaExtents();

    this.currentOp = 'ETRS89ToStereo70';
    this.currentOrder = 'NE';
    this.currentAngularUnit = 'degrees';
    this.pointConverter.setDirection('ETRS89ToStereo70');
    this.pointConverter.setCoordinateOrder('NE');
    this.pointConverter.setAngularUnit('degrees');
    this.batchConverter.setDirection('ETRS89ToStereo70');
    this.batchConverter.setCoordinateOrder('NE');
    this.batchConverter.setAngularUnit('degrees');
    this.fileUpload.setDirection('ETRS89ToStereo70');
    this.fileUpload.setCoordinateOrder('NE');
    this.fileUpload.setAngularUnit('degrees');
    this.mapView.setCoordinateOrder('NE');
    this.mapView.setAngularUnit('degrees');
    this.renderSharedConfigBar();

    if (this.feedbackPanel.isOpen()) this.feedbackPanel.close();
    if (this.gridMetadataPanel.isOpen()) this.gridMetadataPanel.close();
    if (this.geodeticHelpPanel.isOpen()) this.geodeticHelpPanel.close();

    this.setBottomSheetState('half');
    this.switchTab('point');
    this.navigateToPage('studio');
  }

  private navigateToPage(page: AppPage): void {
    this.header.setActivePage(page);
    if (page === 'telemetry') {
      this.workbenchLayout.style.display = 'none';
      this.telemetryDashboard.getElement().style.display = 'block';
      this.telemetryDashboard.refresh();
    } else {
      this.telemetryDashboard.getElement().style.display = 'none';
      this.workbenchLayout.style.display = 'flex';
      setTimeout(() => this.mapView.invalidateSize(), 60);
    }
  }
}

// Bootstrap
document.addEventListener('DOMContentLoaded', () => {
  new App();
});
