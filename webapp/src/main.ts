import './style/variables.css';
import './style/base.css';
import './style/components.css';
import './style/map.css';

import { t, onLanguageChange } from './i18n';
import { Header } from './components/Header';
import { PointConverter, type AngularUnitType } from './components/PointConverter';
import { MapView } from './components/MapView';
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

  private header!: Header;
  private pointConverter!: PointConverter;
  private batchConverter!: BatchConverter;
  private fileUpload!: FileUpload;
  private gridMetadataPanel = new GridMetadataPanel();
  private geodeticHelpPanel = new GeodeticHelpPanel();
  private mapView!: MapView;

  private feedbackPanel = new FeedbackPanel();
  private volumeWarningModal = new VolumeWarningModal();

  constructor() {
    this.rootElement = document.getElementById('app-root')!;
    this.initTheme();
    this.initLayout();
    this.mapView.initMap();
    this.mapView.setAngularUnit(this.currentAngularUnit);

    onLanguageChange(() => {
      this.renderSharedConfigBar();
      this.renderTabs();
    });
  }

  private initTheme(): void {
    const savedTheme = localStorage.getItem('pytransdat_theme') || 'dark';
    document.documentElement.setAttribute('data-theme', savedTheme);
  }

  private initLayout(): void {
    // 1. Header
    this.header = new Header({
      onNavigate: (page) => this.navigateToPage(page),
      onOpenFeedback: (initialType) => this.feedbackPanel.toggle(initialType),
      onOpenGridInfo: () => this.gridMetadataPanel.toggle(),
      onOpenGeodeticHelp: () => this.geodeticHelpPanel.toggle(),
      onToggleMobileView: () => this.toggleMobileView(),
    });
    this.rootElement.appendChild(this.header.getElement());

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
    this.utilityPanel.className = 'utility-panel mobile-active';

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

        if (this.currentOp === 'ETRS89ToStereo70') {
          this.pointConverter.setCoordinatesFromMap(lat, lon);
        } else {
          // Stereo 70 ➔ ETRS89: Reverse-calculate clicked lat/lon into Stereo 70 meters
          try {
            const res = await apiClient.transformPoint({
              op: 'ETRS89ToStereo70',
              coos: [lat, lon],
              unit: 'degrees',
            });
            if (res.warning === 'Out of grid') {
              showToast(t('point.outOfGridWarn'), 'warning');
            } else if (res.coos && res.coos.length >= 2) {
              this.pointConverter.setStereoCoordinatesFromMap(res.coos[0], res.coos[1]);
            }
          } catch (err) {
            console.error('Failed to reverse calculate Stereo 70 coordinates', err);
          }
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

    this.pointConverter.setCoordinateOrder(this.currentOrder);
    this.batchConverter.setCoordinateOrder(this.currentOrder);
    this.fileUpload.setCoordinateOrder(this.currentOrder);

    this.updateActiveTabContent();
  }

  private setDirection(op: TransformationOp): void {
    if (this.currentOp === op) return;
    this.currentOp = op;
    this.renderSharedConfigBar();
    this.pointConverter.setDirection(op);
    this.batchConverter.setDirection(op);
    this.fileUpload.setDirection(op);
  }

  private setAngularUnit(unit: AngularUnitType): void {
    if (this.currentAngularUnit === unit) return;
    this.currentAngularUnit = unit;
    this.renderSharedConfigBar();
    this.pointConverter.setAngularUnit(unit);
    this.batchConverter.setAngularUnit(unit);
    this.fileUpload.setAngularUnit(unit);
    this.mapView.setAngularUnit(unit);
  }

  private setCoordinateOrder(order: CoordinateOrder): void {
    if (this.currentOrder === order) return;
    this.currentOrder = order;
    this.renderSharedConfigBar();
    this.pointConverter.setCoordinateOrder(order);
    this.batchConverter.setCoordinateOrder(order);
    this.fileUpload.setCoordinateOrder(order);
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
      <!-- Tier 1: Sliding Direction Toggle (ETRS89 ➔ Stereo 70 first, Stereo 70 ➔ ETRS89 second) -->
      <div class="sliding-toggle" data-direction="${isEtrs ? 'etrs' : 's70'}">
        <div class="sliding-toggle-glider"></div>
        <button class="sliding-toggle-btn ${isEtrs ? 'active' : ''}" data-op="ETRS89ToStereo70">
          ${t('point.opETRSToS70')}
        </button>
        <button class="sliding-toggle-btn ${!isEtrs ? 'active' : ''}" data-op="Stereo70ToETRS89">
          ${t('point.opS70ToETRS')}
        </button>
      </div>

      <!-- Tier 2: 50/50 Split: Angular Unit | Coordinate Order -->
      <div class="shared-config-row">
        <div class="config-field-group">
          <span class="config-field-label">${t('point.unitLabel')}</span>
          <select class="form-select" id="angular-unit-select">
            <option value="degrees" ${this.currentAngularUnit === 'degrees' ? 'selected' : ''}>${t('point.unitDeg')}</option>
            <option value="dms" ${this.currentAngularUnit === 'dms' ? 'selected' : ''}>${t('point.unitDMS')}</option>
            <option value="radians" ${this.currentAngularUnit === 'radians' ? 'selected' : ''}>${t('point.unitRad')}</option>
          </select>
        </div>

        <div class="coord-order-picker">
          <span class="config-field-label">${t('common.coordOrder')}</span>
          <div class="segmented-control" style="font-size: 11px;">
            <button class="segmented-btn ${this.currentOrder === 'NE' ? 'active' : ''}" id="shared-order-ne">NE(H)</button>
            <button class="segmented-btn ${this.currentOrder === 'EN' ? 'active' : ''}" id="shared-order-en">EN(H)</button>
          </div>
        </div>
      </div>

      <!-- Tier 3: Angle Format Example with Isolated Copyable Badge -->
      <div class="format-example-row">
        <span class="format-example-label">${t('point.exampleFormatLabel')}</span>
        <span class="format-example-badge" id="format-example-badge" title="${t('point.clickToSelect')}">${sampleText}</span>
      </div>
    `;

    // Direction Toggle Listeners
    this.sharedConfigBar.querySelectorAll('.sliding-toggle-btn').forEach((btn) => {
      btn.addEventListener('click', () => {
        const op = btn.getAttribute('data-op') as TransformationOp;
        if (op) this.setDirection(op);
      });
    });

    // Angular Unit Select Listener
    const select = this.sharedConfigBar.querySelector('#angular-unit-select') as HTMLSelectElement;
    if (select) {
      select.addEventListener('change', () => {
        const unit = select.value as AngularUnitType;
        if (unit) this.setAngularUnit(unit);
      });
    }

    // Coordinate Order Listeners
    this.sharedConfigBar.querySelector('#shared-order-ne')?.addEventListener('click', () => {
      this.setCoordinateOrder('NE');
    });
    this.sharedConfigBar.querySelector('#shared-order-en')?.addEventListener('click', () => {
      this.setCoordinateOrder('EN');
    });

    // Angle Format Example Badge Selection Listener (selects ONLY the number string)
    const exampleBadge = this.sharedConfigBar.querySelector('#format-example-badge') as HTMLElement;
    if (exampleBadge) {
      exampleBadge.addEventListener('click', () => {
        const range = document.createRange();
        range.selectNodeContents(exampleBadge);
        const sel = window.getSelection();
        sel?.removeAllRanges();
        sel?.addRange(range);
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(sampleText);
        }
        showToast(t('point.copiedToast'), 'success');
      });
    }
  }

  private renderTabs(): void {
    this.tabsBar.innerHTML = `
      <span class="tabs-prefix">${t('tabs.prefix')}</span>
      <button class="tab-btn ${this.activeTab === 'point' ? 'active' : ''}" data-tab="point">
        <span>${t('tabs.point')}</span>
      </button>
      <button class="tab-btn ${this.activeTab === 'batch' ? 'active' : ''}" data-tab="batch">
        <span>${t('tabs.batch')}</span>
      </button>
      <button class="tab-btn ${this.activeTab === 'file' ? 'active' : ''}" data-tab="file">
        <span>${t('tabs.file')}</span>
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

  private toggleMobileView(): void {
    const isUtilActive = this.utilityPanel.classList.contains('mobile-active');
    const mapPanel = this.mapView.getElement();

    if (isUtilActive) {
      this.utilityPanel.classList.remove('mobile-active');
      mapPanel.classList.add('mobile-active');
      setTimeout(() => this.mapView.invalidateSize(), 60);
    } else {
      this.utilityPanel.classList.add('mobile-active');
      mapPanel.classList.remove('mobile-active');
    }
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
