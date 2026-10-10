import './style/variables.css';
import './style/base.css';
import './style/components.css';
import './style/map.css';

import { t, onLanguageChange } from './i18n';
import { Header } from './components/Header';
import { PointConverter } from './components/PointConverter';
import { MapView } from './components/MapView';
import { BatchConverter } from './components/BatchConverter';
import { FileUpload } from './components/FileUpload';
import { StatsView } from './components/StatsView';
import { ReportModal } from './components/ReportModal';
import { FeedbackModal } from './components/FeedbackModal';
import { VolumeWarningModal } from './components/VolumeWarningModal';

type ActiveTab = 'point' | 'batch' | 'file' | 'stats';

class App {
  private activeTab: ActiveTab = 'point';
  private rootElement: HTMLElement;
  private utilityPanel!: HTMLElement;
  private tabsBar!: HTMLElement;
  private tabContentContainer!: HTMLElement;

  private header!: Header;
  private pointConverter!: PointConverter;
  private batchConverter!: BatchConverter;
  private fileUpload!: FileUpload;
  private statsView!: StatsView;
  private mapView!: MapView;

  private reportModal = new ReportModal();
  private feedbackModal = new FeedbackModal();
  private volumeWarningModal = new VolumeWarningModal();

  constructor() {
    this.rootElement = document.getElementById('app-root')!;
    this.initTheme();
    this.initLayout();
    this.mapView.initMap();

    onLanguageChange(() => this.renderTabs());
  }

  private initTheme(): void {
    const savedTheme = localStorage.getItem('pytransdat_theme') || 'dark';
    document.documentElement.setAttribute('data-theme', savedTheme);
  }

  private initLayout(): void {
    // 1. Header
    this.header = new Header({
      onOpenFeedback: () => this.feedbackModal.open(),
      onToggleMobileView: () => this.toggleMobileView(),
    });
    this.rootElement.appendChild(this.header.getElement());

    // 2. Workbench Layout
    const workbench = document.createElement('main');
    workbench.className = 'workbench-layout';

    // Left Panel
    this.utilityPanel = document.createElement('div');
    this.utilityPanel.className = 'utility-panel mobile-active';

    // Tabs Bar
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
      onMapPointClicked: (lat, lon) => {
        // Switch to point tab and set coords
        this.switchTab('point');
        this.pointConverter.setCoordinatesFromMap(lat, lon);
      },
      onOpenVolumeWarning: (count, callback) => {
        this.volumeWarningModal.open(count, callback);
      },
    });

    workbench.appendChild(this.mapView.getElement());
    this.rootElement.appendChild(workbench);

    // 3. Child Components
    this.pointConverter = new PointConverter({
      onPointTransformed: ({ lat, lon, outOfGrid }) => {
        this.mapView.updateSinglePoint(lat, lon, outOfGrid);
      },
    });

    this.batchConverter = new BatchConverter({
      onPointsPlotted: (points) => {
        this.mapView.plotBatchPoints(points);
      },
      onOpenReport: (report) => {
        this.reportModal.open(report);
      },
    });

    this.fileUpload = new FileUpload();
    this.statsView = new StatsView();

    this.updateActiveTabContent();
  }

  private renderTabs(): void {
    this.tabsBar.innerHTML = `
      <button class="tab-btn ${this.activeTab === 'point' ? 'active' : ''}" data-tab="point">
        <span>${t('tabs.point')}</span>
      </button>
      <button class="tab-btn ${this.activeTab === 'batch' ? 'active' : ''}" data-tab="batch">
        <span>${t('tabs.batch')}</span>
      </button>
      <button class="tab-btn ${this.activeTab === 'file' ? 'active' : ''}" data-tab="file">
        <span>${t('tabs.file')}</span>
      </button>
      <button class="tab-btn ${this.activeTab === 'stats' ? 'active' : ''}" data-tab="stats">
        <span>${t('tabs.stats')}</span>
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
      case 'stats':
        this.tabContentContainer.appendChild(this.statsView.getElement());
        break;
    }
  }

  private toggleMobileView(): void {
    const isUtilActive = this.utilityPanel.classList.contains('mobile-active');
    const mapPanel = this.mapView.getElement();

    if (isUtilActive) {
      this.utilityPanel.classList.remove('mobile-active');
      mapPanel.classList.add('mobile-active');
    } else {
      this.utilityPanel.classList.add('mobile-active');
      mapPanel.classList.remove('mobile-active');
    }
  }
}

// Bootstrap
document.addEventListener('DOMContentLoaded', () => {
  new App();
});
