import L from 'leaflet';
import 'leaflet/dist/leaflet.css';
import 'leaflet.markercluster';
import 'leaflet.markercluster/dist/MarkerCluster.css';
import 'leaflet.markercluster/dist/MarkerCluster.Default.css';

import { t, getLanguage, onLanguageChange } from '../i18n';
import { APP_CONFIG } from '../config';
import { formatDegToDMS, degToRad } from '../utils/dmsParser';
import type { AngularUnitType, PointTransformedEvent } from './PointConverter';
import type { CoordinateOrder } from '../utils/pointParser';
import type { BoundsBox, MapViewSnapshot } from './mapViewSnapshot';

export const ROMANIA_GRID_BOUNDS = {
  minLat: 43.5,
  maxLat: 48.3,
  minLon: 20.2,
  maxLon: 29.8,
};

export function isInsideRomaniaGrid(lat: number, lon: number): boolean {
  return (
    !isNaN(lat) &&
    !isNaN(lon) &&
    lat >= ROMANIA_GRID_BOUNDS.minLat &&
    lat <= ROMANIA_GRID_BOUNDS.maxLat &&
    lon >= ROMANIA_GRID_BOUNDS.minLon &&
    lon <= ROMANIA_GRID_BOUNDS.maxLon
  );
}

export const LAYER_COLOR_PRESETS = [
  { name: 'Cyan', color: '#38bdf8' },
  { name: 'Cobalt', color: '#2563eb' },
  { name: 'Emerald', color: '#10b981' },
  { name: 'Amber', color: '#f59e0b' },
  { name: 'Violet', color: '#8b5cf6' },
  { name: 'Rose', color: '#f43f5e' },
];

export interface MapPointItem {
  id: string | null;
  lat: number;
  lon: number;
  elev?: number;
  altitude?: number;
  ellipsoidalHeight?: number;
  stereoN?: number;
  stereoE?: number;
  zeta?: number | null;
  warning?: string | null;
  fileName?: string;
  angularUnit?: AngularUnitType;
}

export interface MapViewCallbacks {
  onMapPointClicked: (lat: number, lon: number) => void;
  onOpenVolumeWarning: (
    count: number,
    onConfirm: (choice: 'all' | 'sample' | 'skip') => void
  ) => void;
}

export class MapView {
  private container: HTMLElement;
  private mapElement: HTMLElement;
  private coordsBadge!: HTMLElement;
  private callbacks: MapViewCallbacks;

  private map: L.Map | null = null;
  private singlePointMarker: L.Marker | null = null;
  private pointLayer: L.LayerGroup | null = null;
  private textLayer: any = null; // L.MarkerClusterGroup
  private filesLayer: any = null; // L.MarkerClusterGroup

  // Measurement tool state: 'idle' | 'measuring' | 'completed'
  private measureLayer: L.LayerGroup = L.layerGroup();
  private isMeasuringMode: boolean = false;
  private measureState: 'idle' | 'measuring' | 'completed' = 'idle';
  private measureStartLatLng: L.LatLng | null = null;
  private measureCurrentLine: L.Polyline | null = null;
  private measureTooltip: L.Tooltip | null = null;

  // Basemaps (OpenStreetMap selected by default)
  private cartoDark: L.TileLayer | null = null;
  private cartoPositron: L.TileLayer | null = null;
  private osm: L.TileLayer | null = null;
  private satellite: L.TileLayer | null = null;
  private activeBasemapKey: 'dark' | 'positron' | 'osm' | 'satellite' = 'osm';

  // Custom Docks
  private controlsDock!: HTMLElement;
  private layersDock!: HTMLElement;
  private topRightDock!: HTMLElement;
  private legendWidget!: HTMLElement;
  private scaleBarWidget!: HTMLElement;

  // Configuration & State
  private currentAngularUnit: AngularUnitType = 'degrees';
  private currentOrder: CoordinateOrder = 'NE';
  private lastCursorLatLng: L.LatLng | null = null;
  private isLegendOpen: boolean = localStorage.getItem('pytransdat_legend_open') !== 'false';

  // Layer state & colors
  private layerColors: { point: string; text: string; files: string } = {
    point: '#38bdf8',
    text: '#2563eb',
    files: '#10b981',
  };
  private layerVisibility = {
    point: true,
    text: true,
    files: true,
  };
  private lastTextPoints: MapPointItem[] = [];
  private lastFilePoints: MapPointItem[] = [];

  constructor(callbacks: MapViewCallbacks) {
    this.callbacks = callbacks;

    this.loadSavedColors();
    this.loadSavedBasemap();

    this.container = document.createElement('div');
    this.container.className = 'map-panel';

    this.mapElement = document.createElement('div');
    this.mapElement.id = 'map';
    this.container.appendChild(this.mapElement);

    this.initCustomDocks();

    onLanguageChange(() => {
      this.refreshDocksContent();
      this.updateCursorBadge();
    });

    // Sync default basemap with theme
    const observer = new MutationObserver((mutations) => {
      for (const m of mutations) {
        if (m.type === 'attributes' && m.attributeName === 'data-theme') {
          const isDarkTheme = document.documentElement.getAttribute('data-theme') === 'dark';
          if (this.activeBasemapKey === 'dark' || this.activeBasemapKey === 'positron') {
            this.switchBasemap(isDarkTheme ? 'dark' : 'positron');
          }
        }
      }
    });
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  }

  private loadSavedBasemap(): void {
    try {
      const saved = localStorage.getItem('pytransdat_basemap');
      if (saved && ['dark', 'positron', 'osm', 'satellite'].includes(saved)) {
        this.activeBasemapKey = saved as any;
      }
    } catch {
      // Default osm
    }
  }

  private loadSavedColors(): void {
    try {
      const saved = localStorage.getItem('pytransdat_map_colors');
      if (saved) {
        const parsed = JSON.parse(saved);
        if (parsed.point) this.layerColors.point = parsed.point;
        if (parsed.text) this.layerColors.text = parsed.text;
        if (parsed.files) this.layerColors.files = parsed.files;
      }
    } catch {
      // Use defaults
    }
    this.applyLayerColors('point');
  }

  private saveColors(): void {
    try {
      localStorage.setItem('pytransdat_map_colors', JSON.stringify(this.layerColors));
    } catch {
      // Ignore
    }
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  /** Returns a plain-data snapshot of the map state (no live Leaflet objects). */
  public getSnapshot(): MapViewSnapshot {
    if (!this.map) {
      throw new Error('MapView.getSnapshot() called before initMap()');
    }
    const toBox = (b: L.LatLngBounds): BoundsBox => ({
      south: b.getSouth(),
      west: b.getWest(),
      north: b.getNorth(),
      east: b.getEast(),
    });
    const center = this.map.getCenter();
    const maxBounds = this.map.options.maxBounds;
    const point = this.pointLayer ? this.pointLayer.getLayers().length : 0;
    const text = this.textLayer ? this.textLayer.getLayers().length : 0;
    const files = this.filesLayer ? this.filesLayer.getLayers().length : 0;
    return {
      center: { lat: center.lat, lng: center.lng },
      zoom: this.map.getZoom(),
      bounds: toBox(this.map.getBounds()),
      maxBounds: maxBounds ? toBox(L.latLngBounds(maxBounds as L.LatLngBoundsLiteral)) : null,
      markerCount: point + text + files,
      markersByLayer: { point, text, files },
    };
  }

  public fireMapClick(lat: number, lng: number): void {
    if (this.isMeasuringMode) {
      this.handleMeasureClick(L.latLng(lat, lng));
      return;
    }
    this.callbacks.onMapPointClicked(lat, lng);
  }

  public invalidateSize(): void {
    if (this.map) {
      this.map.invalidateSize();
    }
  }

  public setAngularUnit(unit: AngularUnitType): void {
    this.currentAngularUnit = unit;
    this.updateCursorBadge();
  }

  public setCoordinateOrder(order: CoordinateOrder): void {
    if (this.currentOrder === order) return;
    this.currentOrder = order;
    this.updateCursorBadge();
  }

  public initMap(): void {
    if (this.map) return;

    // Center and zoom initialized from centralized APP_CONFIG
    this.map = L.map(this.mapElement, {
      center: APP_CONFIG.defaultMapCenter,
      zoom: APP_CONFIG.defaultMapZoom,
      minZoom: 3,
      maxZoom: 18,
      zoomControl: false, // Custom Controls Dock in top-left
    });

    // Basemaps initialization
    this.cartoDark = L.tileLayer(
      'https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
      {
        attribution: '&copy; CARTO',
        subdomains: 'abcd',
        maxZoom: 19,
      }
    );

    this.cartoPositron = L.tileLayer(
      'https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
      {
        attribution: '&copy; CARTO',
        subdomains: 'abcd',
        maxZoom: 19,
      }
    );

    this.osm = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; OpenStreetMap',
      maxZoom: 19,
    });

    this.satellite = L.tileLayer(
      'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
      {
        attribution: '&copy; Esri',
        maxZoom: 18,
      }
    );

    this.switchBasemap(this.activeBasemapKey);

    // Initialize data layers + measurement layer with color-synchronized cluster icons
    this.pointLayer = L.layerGroup().addTo(this.map);
    this.textLayer = this.createClusterGroup('text').addTo(this.map);
    this.filesLayer = this.createClusterGroup('files').addTo(this.map);
    this.measureLayer.addTo(this.map);

    // GIS Graphic Scale Bar in Top-Right
    this.updateGisScale();
    this.map.on('zoomend moveend viewreset', () => {
      this.updateGisScale();
    });

    // Safe Popup Copy handler using event listener (no inline onclick escaping issues)
    this.map.on('popupopen', (e: L.PopupEvent) => {
      const container = e.popup.getElement();
      if (!container) return;
      const copyBtn = container.querySelector<HTMLButtonElement>('.geodetic-popup-copy-btn');
      if (copyBtn && !copyBtn.dataset.bound) {
        copyBtn.dataset.bound = 'true';
        copyBtn.addEventListener('click', (ev) => {
          ev.preventDefault();
          ev.stopPropagation();
          const raw = copyBtn.getAttribute('data-copy-text');
          if (!raw) return;
          const textToCopy = decodeURIComponent(raw);
          navigator.clipboard.writeText(textToCopy).then(() => {
            const origHtml = copyBtn.innerHTML;
            copyBtn.innerHTML = `<span>${t('map.copiedAttributesFeedback')}</span>`;
            copyBtn.classList.add('copied');
            setTimeout(() => {
              copyBtn.innerHTML = origHtml;
              copyBtn.classList.remove('copied');
            }, 1500);
          }).catch((err) => {
            console.error('Clipboard copy failed', err);
          });
        });
      }
    });

    // Map Event Listeners
    this.map.on('mousemove', (e: L.LeafletMouseEvent) => {
      this.lastCursorLatLng = e.latlng;
      this.updateCursorBadge();
      if (this.isMeasuringMode && this.measureState === 'measuring' && this.measureStartLatLng) {
        this.updateMeasuringLine(e.latlng);
      }
    });

    this.map.on('click', (e: L.LeafletMouseEvent) => {
      if (this.isMeasuringMode) {
        this.handleMeasureClick(e.latlng);
        return;
      }
      this.callbacks.onMapPointClicked(e.latlng.lat, e.latlng.lng);
    });

    // Zoom level CAD point labels visibility toggle
    this.map.on('zoomend', () => {
      if (!this.map) return;
      if (this.map.getZoom() >= 15) {
        this.mapElement.classList.add('zoom-labels-visible');
      } else {
        this.mapElement.classList.remove('zoom-labels-visible');
      }
    });

    // Esc key cancels measurement mode
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape' && this.isMeasuringMode) {
        this.toggleMeasureMode(false);
      }
    });

    this.refreshDocksContent();
  }

  private createClusterGroup(layerKey: 'text' | 'files'): any {
    return (L as any).markerClusterGroup({
      showCoverageOnHover: false,
      maxClusterRadius: 50,
      spiderfyOnMaxZoom: true,
      iconCreateFunction: (cluster: any) => {
        const count = cluster.getChildCount();
        const color = this.layerColors[layerKey];
        return L.divIcon({
          html: `<div class="custom-cluster-circle" style="background-color: ${color};"><span style="color: #ffffff;">${count}</span></div>`,
          className: 'custom-marker-cluster',
          iconSize: new L.Point(36, 36),
        });
      },
    });
  }

  private initCustomDocks(): void {
    // 1. Top-Left Controls Dock
    this.controlsDock = document.createElement('div');
    this.controlsDock.className = 'map-controls-dock';
    this.container.appendChild(this.controlsDock);

    // 2. Top-Right Stack (Legend + Graphic Scale Bar)
    this.topRightDock = document.createElement('div');
    this.topRightDock.className = 'map-top-right-dock';

    this.legendWidget = document.createElement('div');
    this.legendWidget.className = 'map-legend-dock';
    this.topRightDock.appendChild(this.legendWidget);

    this.scaleBarWidget = document.createElement('div');
    this.scaleBarWidget.className = 'map-scale-dock';
    this.topRightDock.appendChild(this.scaleBarWidget);

    this.container.appendChild(this.topRightDock);

    // 3. Bottom-Left Layers Button & Panel
    this.layersDock = document.createElement('div');
    this.layersDock.className = 'map-layers-dock';
    this.container.appendChild(this.layersDock);

    // 4. Bottom-Center Coordinates Badge
    this.coordsBadge = document.createElement('div');
    this.coordsBadge.className = 'map-coords-badge';
    this.container.appendChild(this.coordsBadge);

    this.renderControlsDock();
    this.renderLegendWidget();
    this.renderLayersDock();
    this.updateCursorBadge();

    // One-time global outside-click listener for layers panel and swatches
    document.addEventListener('click', (e) => {
      if (!this.layersDock.contains(e.target as Node)) {
        const panel = this.layersDock.querySelector('#map-layers-panel');
        panel?.classList.remove('open');
        this.layersDock.querySelectorAll('.color-swatch-popover').forEach((p) => p.classList.remove('open'));
      }
    });

    // Synchronize fullscreen state when entering/exiting (including Escape or browser shortcuts)
    document.addEventListener('fullscreenchange', () => {
      const isFs = !!document.fullscreenElement;
      const fsBtn = this.controlsDock.querySelector('#map-fullscreen-btn');
      fsBtn?.classList.toggle('active', isFs);
      this.invalidateSize();
      setTimeout(() => this.invalidateSize(), 150);
    });
  }

  private refreshDocksContent(): void {
    this.renderControlsDock();
    this.renderLegendWidget();
    this.renderLayersDock();
    this.updateGisScale();
  }

  private updateCursorBadge(): void {
    if (!this.lastCursorLatLng) {
      this.coordsBadge.textContent =
        this.currentOrder === 'EN'
          ? `Long: 24.96680000  Lat: 45.94320000`
          : `Lat: 45.94320000  Long: 24.96680000`;
      return;
    }

    const lat = this.lastCursorLatLng.lat;
    const lon = this.lastCursorLatLng.lng;

    let latStr: string;
    let lonStr: string;

    if (this.currentAngularUnit === 'dms') {
      latStr = formatDegToDMS(lat);
      lonStr = formatDegToDMS(lon);
    } else if (this.currentAngularUnit === 'radians') {
      latStr = degToRad(lat).toFixed(10);
      lonStr = degToRad(lon).toFixed(10);
    } else {
      latStr = lat.toFixed(8);
      lonStr = lon.toFixed(8);
    }

    this.coordsBadge.textContent =
      this.currentOrder === 'EN'
        ? `Long: ${lonStr}  Lat: ${latStr}`
        : `Lat: ${latStr}  Long: ${lonStr}`;
  }

  private renderControlsDock(): void {
    this.controlsDock.innerHTML = `
      <div class="map-control-btn-group">
        <button class="map-dock-btn" id="map-zoom-in-btn" title="Zoom in">
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="12" y1="5" x2="12" y2="19"></line>
            <line x1="5" y1="12" x2="19" y2="12"></line>
          </svg>
        </button>
        <button class="map-dock-btn" id="map-zoom-out-btn" title="Zoom out">
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="5" y1="12" x2="19" y2="12"></line>
          </svg>
        </button>
      </div>

      <div class="map-control-btn-group">
        <button class="map-dock-btn" id="map-home-btn" title="${t('map.zoomToRomania')}">
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M3 9l9-7 9 7v11a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"></path>
            <polyline points="9 22 9 12 15 12 15 22"></polyline>
          </svg>
        </button>
        <button class="map-dock-btn ${this.isMeasuringMode ? 'active' : ''}" id="map-measure-btn" title="${t('map.measure')}">
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M21.3 8.7l-6-6a2 2 0 0 0-2.8 0L2.7 12.5a2 2 0 0 0 0 2.8l6 6a2 2 0 0 0 2.8 0l9.8-9.8a2 2 0 0 0 0-2.8z"></path>
            <line x1="6.5" y1="11.5" x2="8.5" y2="13.5"></line>
            <line x1="10" y1="8" x2="12" y2="10"></line>
            <line x1="13.5" y1="4.5" x2="15.5" y2="6.5"></line>
          </svg>
        </button>
        <button class="map-dock-btn" id="map-fullscreen-btn" title="${t('map.fullscreen')}">
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3"></path>
          </svg>
        </button>
      </div>
    `;

    this.controlsDock.querySelector('#map-zoom-in-btn')?.addEventListener('click', () => {
      this.map?.zoomIn();
    });

    this.controlsDock.querySelector('#map-zoom-out-btn')?.addEventListener('click', () => {
      this.map?.zoomOut();
    });

    this.controlsDock.querySelector('#map-home-btn')?.addEventListener('click', () => {
      this.zoomToRomaniaExtents();
    });

    this.controlsDock.querySelector('#map-measure-btn')?.addEventListener('click', () => {
      this.toggleMeasureMode();
    });

    this.controlsDock.querySelector('#map-fullscreen-btn')?.addEventListener('click', () => {
      this.toggleFullscreen();
    });
  }

  private renderLegendWidget(): void {
    this.legendWidget.innerHTML = `
      <button class="map-dock-toggle-btn ${this.isLegendOpen ? 'active' : ''}" id="map-legend-toggle" title="${t('map.legendTitle')}">
        <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <line x1="8" y1="6" x2="21" y2="6"></line>
          <line x1="8" y1="12" x2="21" y2="12"></line>
          <line x1="8" y1="18" x2="21" y2="18"></line>
          <line x1="3" y1="6" x2="3.01" y2="6"></line>
          <line x1="3" y1="12" x2="3.01" y2="12"></line>
          <line x1="3" y1="18" x2="3.01" y2="18"></line>
        </svg>
      </button>

      <div class="map-legend-panel ${this.isLegendOpen ? 'open' : ''}" id="map-legend-panel">
        <div class="map-legend-header">
          <span class="map-legend-title">${t('map.legendTitle')}</span>
          <button class="map-legend-close-btn" id="map-legend-close" title="${t('common.close')}">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">
              <line x1="18" y1="6" x2="6" y2="18"></line>
              <line x1="6" y1="6" x2="18" y2="18"></line>
            </svg>
          </button>
        </div>
        <div class="map-legend-items">
          <div class="map-legend-item ${this.layerVisibility.point ? '' : 'dimmed'}">
            <span class="legend-color-dot" style="background-color: ${this.layerColors.point};"></span>
            <span class="legend-label">${t('map.singlePointLayer')}</span>
          </div>
          <div class="map-legend-item ${this.layerVisibility.text ? '' : 'dimmed'}">
            <span class="legend-color-dot" style="background-color: ${this.layerColors.text};"></span>
            <span class="legend-label">${t('map.textLayer')}</span>
          </div>
          <div class="map-legend-item ${this.layerVisibility.files ? '' : 'dimmed'}">
            <span class="legend-color-dot" style="background-color: ${this.layerColors.files};"></span>
            <span class="legend-label">${t('map.filesLayer')}</span>
          </div>
        </div>
      </div>
    `;

    const toggleBtn = this.legendWidget.querySelector('#map-legend-toggle');
    const panel = this.legendWidget.querySelector('#map-legend-panel');
    const closeBtn = this.legendWidget.querySelector('#map-legend-close');

    toggleBtn?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.isLegendOpen = !this.isLegendOpen;
      localStorage.setItem('pytransdat_legend_open', String(this.isLegendOpen));
      panel?.classList.toggle('open', this.isLegendOpen);
      toggleBtn.classList.toggle('active', this.isLegendOpen);
    });

    closeBtn?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.isLegendOpen = false;
      localStorage.setItem('pytransdat_legend_open', 'false');
      panel?.classList.remove('open');
      toggleBtn?.classList.remove('active');
    });
  }

  private updateGisScale(): void {
    if (!this.map || !this.scaleBarWidget) return;

    const lat = this.map.getCenter().lat;
    const zoom = this.map.getZoom();
    const metersPerPixel = (40075016.686 * Math.cos((lat * Math.PI) / 180)) / Math.pow(2, zoom + 8);

    const targetWidthPx = 100;
    const targetMeters = targetWidthPx * metersPerPixel;

    const candidates = [
      1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000,
    ];

    let chosen = candidates[0];
    for (const c of candidates) {
      if (c <= targetMeters * 1.3) {
        chosen = c;
      } else {
        break;
      }
    }

    const actualWidthPx = Math.max(64, Math.round(chosen / metersPerPixel));

    let maxLabel: string;
    let unitLabel: string;

    if (chosen >= 1000) {
      const kmMax = chosen / 1000;
      const isRo = getLanguage() === 'ro';
      maxLabel = kmMax % 1 === 0 ? kmMax.toFixed(0) : (isRo ? kmMax.toFixed(1).replace('.', ',') : kmMax.toFixed(1));
      unitLabel = 'km';
    } else {
      maxLabel = chosen.toFixed(0);
      unitLabel = 'm';
    }

    this.scaleBarWidget.innerHTML = `
      <div class="gis-scale-container" style="width: ${actualWidthPx}px;">
        <div class="gis-scale-ticks">
          <span class="gis-tick tick-start"></span>
          <span class="gis-tick tick-mid"></span>
          <span class="gis-tick tick-end"></span>
        </div>
        <div class="gis-scale-bar">
          <div class="gis-scale-seg seg-first"></div>
          <div class="gis-scale-seg seg-second"></div>
        </div>
        <div class="gis-scale-labels">
          <span class="gis-lbl lbl-start">0</span>
          <span class="gis-lbl lbl-end">${maxLabel} ${unitLabel}</span>
        </div>
      </div>
    `;
  }

  private renderLayersDock(): void {
    const isPanelOpen = this.layersDock.querySelector('.map-layers-panel')?.classList.contains('open');

    this.layersDock.innerHTML = `
      <button class="map-dock-toggle-btn" id="map-layers-toggle" title="${t('map.layersTitle')}">
        <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <polygon points="12 2 2 7 12 12 22 7 12 2"></polygon>
          <polyline points="2 17 12 22 22 17"></polyline>
          <polyline points="2 12 12 17 22 12"></polyline>
        </svg>
      </button>

      <div class="map-layers-panel ${isPanelOpen ? 'open' : ''}" id="map-layers-panel">
        <!-- 1. Data Layers with Counters & Color Pickers -->
        <div>
          <div class="map-layers-section-title">${t('map.dataLayers')}</div>

          <!-- Single Point -->
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="checkbox" class="map-layer-checkbox" id="layer-chk-point" ${this.layerVisibility.point ? 'checked' : ''} />
              <span>${t('map.singlePointLayer')}</span>
              <span class="layer-count-badge">(${this.singlePointMarker ? '1 pt' : '0'})</span>
            </label>
            <div style="position: relative;">
              <button class="layer-color-trigger" style="background-color: ${this.layerColors.point};" data-color-for="point" title="${t('map.colorPicker')}"></button>
              <div class="color-swatch-popover" id="swatch-popover-point">
                ${LAYER_COLOR_PRESETS.map((p) => `
                  <div class="color-swatch-option" style="background-color: ${p.color};" data-layer="point" data-color="${p.color}"></div>
                `).join('')}
              </div>
            </div>
          </div>

          <!-- Batch Text -->
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="checkbox" class="map-layer-checkbox" id="layer-chk-text" ${this.layerVisibility.text ? 'checked' : ''} />
              <span>${t('map.textLayer')}</span>
              <span class="layer-count-badge">(${this.lastTextPoints.length > 0 ? `${this.lastTextPoints.length} pts` : '0'})</span>
            </label>
            <div style="position: relative;">
              <button class="layer-color-trigger" style="background-color: ${this.layerColors.text};" data-color-for="text" title="${t('map.colorPicker')}"></button>
              <div class="color-swatch-popover" id="swatch-popover-text">
                ${LAYER_COLOR_PRESETS.map((p) => `
                  <div class="color-swatch-option" style="background-color: ${p.color};" data-layer="text" data-color="${p.color}"></div>
                `).join('')}
              </div>
            </div>
          </div>

          <!-- Uploaded Files -->
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="checkbox" class="map-layer-checkbox" id="layer-chk-files" ${this.layerVisibility.files ? 'checked' : ''} />
              <span>${t('map.filesLayer')}</span>
              <span class="layer-count-badge">(${this.lastFilePoints.length > 0 ? `${this.lastFilePoints.length} pts` : '0'})</span>
            </label>
            <div style="position: relative;">
              <button class="layer-color-trigger" style="background-color: ${this.layerColors.files};" data-color-for="files" title="${t('map.colorPicker')}"></button>
              <div class="color-swatch-popover" id="swatch-popover-files">
                ${LAYER_COLOR_PRESETS.map((p) => `
                  <div class="color-swatch-option" style="background-color: ${p.color};" data-layer="files" data-color="${p.color}"></div>
                `).join('')}
              </div>
            </div>
          </div>
        </div>

        <!-- 2. Shift Grids (Coming Soon) -->
        <div>
          <div class="map-layers-section-title">${t('map.shiftGridsSection')}</div>
          <div class="map-layer-item disabled-layer">
            <span class="map-layer-left">
              <span class="disabled-bullet">○</span>
              <span>${t('map.shiftsEast')}</span>
            </span>
            <span class="layer-badge-soon">${t('map.comingSoon')}</span>
          </div>
          <div class="map-layer-item disabled-layer">
            <span class="map-layer-left">
              <span class="disabled-bullet">○</span>
              <span>${t('map.shiftsNorth')}</span>
            </span>
            <span class="layer-badge-soon">${t('map.comingSoon')}</span>
          </div>
          <div class="map-layer-item disabled-layer">
            <span class="map-layer-left">
              <span class="disabled-bullet">○</span>
              <span>${t('map.shiftsVector')}</span>
            </span>
            <span class="layer-badge-soon">${t('map.comingSoon')}</span>
          </div>
          <div class="map-layer-item disabled-layer">
            <span class="map-layer-left">
              <span class="disabled-bullet">○</span>
              <span>${t('map.shiftsHeight')}</span>
            </span>
            <span class="layer-badge-soon">${t('map.comingSoon')}</span>
          </div>
        </div>

        <!-- 3. Base Maps (Styled like all other layers with Radio Checkmarks) -->
        <div>
          <div class="map-layers-section-title">${t('map.baseLayers')}</div>
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="radio" name="basemap-choice" class="map-layer-checkbox" value="osm" ${this.activeBasemapKey === 'osm' ? 'checked' : ''} />
              <span>${t('map.baseOSM')}</span>
            </label>
          </div>
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="radio" name="basemap-choice" class="map-layer-checkbox" value="dark" ${this.activeBasemapKey === 'dark' ? 'checked' : ''} />
              <span>${t('map.baseDark')}</span>
            </label>
          </div>
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="radio" name="basemap-choice" class="map-layer-checkbox" value="positron" ${this.activeBasemapKey === 'positron' ? 'checked' : ''} />
              <span>${t('map.basePositron')}</span>
            </label>
          </div>
          <div class="map-layer-item">
            <label class="map-layer-left">
              <input type="radio" name="basemap-choice" class="map-layer-checkbox" value="satellite" ${this.activeBasemapKey === 'satellite' ? 'checked' : ''} />
              <span>${t('map.baseSatellite')}</span>
            </label>
          </div>
        </div>
      </div>
    `;

    // Toggle panel open/close
    const toggleBtn = this.layersDock.querySelector('#map-layers-toggle');
    const panel = this.layersDock.querySelector('#map-layers-panel');
    toggleBtn?.addEventListener('click', (e) => {
      e.stopPropagation();
      panel?.classList.toggle('open');
    });

    // Basemaps radio listeners
    this.layersDock.querySelectorAll('input[name="basemap-choice"]').forEach((radio) => {
      radio.addEventListener('change', (e: any) => {
        if (e.target.checked) {
          const key = e.target.value as any;
          if (key) this.switchBasemap(key);
        }
      });
    });

    // Layer visibility checkboxes
    this.layersDock.querySelector('#layer-chk-point')?.addEventListener('change', (e: any) => {
      this.layerVisibility.point = e.target.checked;
      if (this.pointLayer && this.map) {
        if (e.target.checked) this.map.addLayer(this.pointLayer);
        else this.map.removeLayer(this.pointLayer);
      }
      this.renderLegendWidget();
    });

    this.layersDock.querySelector('#layer-chk-text')?.addEventListener('change', (e: any) => {
      this.layerVisibility.text = e.target.checked;
      if (this.textLayer && this.map) {
        if (e.target.checked) this.map.addLayer(this.textLayer);
        else this.map.removeLayer(this.textLayer);
      }
      this.renderLegendWidget();
    });

    this.layersDock.querySelector('#layer-chk-files')?.addEventListener('change', (e: any) => {
      this.layerVisibility.files = e.target.checked;
      if (this.filesLayer && this.map) {
        if (e.target.checked) this.map.addLayer(this.filesLayer);
        else this.map.removeLayer(this.filesLayer);
      }
      this.renderLegendWidget();
    });

    // Color Swatch triggers
    this.layersDock.querySelectorAll('.layer-color-trigger').forEach((trigger) => {
      trigger.addEventListener('click', (e) => {
        e.stopPropagation();
        const layerKey = trigger.getAttribute('data-color-for');
        const popover = this.layersDock.querySelector(`#swatch-popover-${layerKey}`);
        this.layersDock.querySelectorAll('.color-swatch-popover').forEach((p) => {
          if (p !== popover) p.classList.remove('open');
        });
        popover?.classList.toggle('open');
      });
    });

    // Swatch selection
    this.layersDock.querySelectorAll('.color-swatch-option').forEach((opt) => {
      opt.addEventListener('click', (e) => {
        e.stopPropagation();
        const layerKey = opt.getAttribute('data-layer') as 'point' | 'text' | 'files';
        const color = opt.getAttribute('data-color');
        if (layerKey && color) {
          this.layerColors[layerKey] = color;
          this.saveColors();
          this.layersDock.querySelectorAll('.color-swatch-popover').forEach((p) => p.classList.remove('open'));
          this.renderLayersDock();
          this.renderLegendWidget();
          this.applyLayerColors(layerKey);
        }
      });
    });
  }

  private switchBasemap(key: 'dark' | 'positron' | 'osm' | 'satellite'): void {
    if (!this.map) return;
    this.activeBasemapKey = key;
    try {
      localStorage.setItem('pytransdat_basemap', key);
    } catch {
      // Ignore storage errors
    }

    const layers = [this.cartoDark, this.cartoPositron, this.osm, this.satellite];
    layers.forEach((l) => {
      if (l && this.map!.hasLayer(l)) this.map!.removeLayer(l);
    });

    let target: L.TileLayer | null = null;
    if (key === 'dark') target = this.cartoDark;
    else if (key === 'positron') target = this.cartoPositron;
    else if (key === 'osm') target = this.osm;
    else if (key === 'satellite') target = this.satellite;

    if (target) target.addTo(this.map);
    this.renderLayersDock();
  }

  private applyLayerColors(layerKey: 'point' | 'text' | 'files'): void {
    if (layerKey === 'point') {
      document.documentElement.style.setProperty('--point-marker-color', this.layerColors.point);
      document.documentElement.style.setProperty('--point-pulse-color', `${this.layerColors.point}66`);
    } else if (layerKey === 'text') {
      if (this.lastTextPoints.length > 0) {
        this.recreateClusterLayer('text', this.lastTextPoints);
      }
    } else if (layerKey === 'files') {
      if (this.lastFilePoints.length > 0) {
        this.recreateClusterLayer('files', this.lastFilePoints);
      }
    }
  }

  private recreateClusterLayer(layerKey: 'text' | 'files', points: MapPointItem[]): void {
    if (!this.map) return;
    const oldLayer = layerKey === 'text' ? this.textLayer : this.filesLayer;
    if (oldLayer && this.map.hasLayer(oldLayer)) {
      this.map.removeLayer(oldLayer);
    }

    const newLayer = this.createClusterGroup(layerKey);
    if (layerKey === 'text') this.textLayer = newLayer;
    else this.filesLayer = newLayer;

    this.renderPointsIntoLayer(newLayer, points, this.layerColors[layerKey], layerKey === 'files');

    if (this.layerVisibility[layerKey]) {
      newLayer.addTo(this.map);
    }
  }

  public zoomToRomaniaExtents(): void {
    if (!this.map) return;
    this.map.fitBounds(
      [
        [43.5, 20.2],
        [48.3, 29.8],
      ],
      { padding: [30, 30], animate: true }
    );
  }

  private toggleFullscreen(): void {
    if (!document.fullscreenElement) {
      this.container.requestFullscreen?.().catch(() => {});
    } else {
      document.exitFullscreen?.().catch(() => {});
    }
    setTimeout(() => this.invalidateSize(), 100);
  }

  // Simplified Distance Measurement Tool (Measures one distance at a time)
  public toggleMeasureMode(forceState?: boolean): void {
    this.isMeasuringMode = forceState !== undefined ? forceState : !this.isMeasuringMode;
    const btn = this.controlsDock.querySelector('#map-measure-btn');

    if (this.isMeasuringMode) {
      btn?.classList.add('active');
      this.mapElement.classList.add('measuring-mode');
      this.clearMeasurement();
      this.measureState = 'measuring';
    } else {
      btn?.classList.remove('active');
      this.mapElement.classList.remove('measuring-mode');
      this.clearMeasurement();
      this.measureState = 'idle';
    }
  }

  private handleMeasureClick(latlng: L.LatLng): void {
    // If previously completed, subsequent click resets previous line and starts a new measurement cleanly
    if (this.measureState === 'completed') {
      this.clearMeasurement();
      this.measureState = 'measuring';
    }

    if (!this.measureStartLatLng) {
      // Step 1: Start point
      this.measureStartLatLng = latlng;
      const startMarker = L.circleMarker(latlng, {
        radius: 5,
        color: '#ffffff',
        fillColor: '#38bdf8',
        weight: 2,
        fillOpacity: 1,
      });
      this.measureLayer.addLayer(startMarker);
    } else {
      // Step 2: Lock measurement
      const dist = this.measureStartLatLng.distanceTo(latlng);
      const isRo = getLanguage() === 'ro';
      const distStr = dist < 1000
        ? `${isRo ? dist.toFixed(1).replace('.', ',') : dist.toFixed(1)} m`
        : `${isRo ? (dist / 1000).toFixed(3).replace('.', ',') : (dist / 1000).toFixed(3)} km`;

      const endMarker = L.circleMarker(latlng, {
        radius: 5,
        color: '#ffffff',
        fillColor: '#38bdf8',
        weight: 2,
        fillOpacity: 1,
      });
      this.measureLayer.addLayer(endMarker);

      const finalLine = L.polyline([this.measureStartLatLng, latlng], {
        color: '#38bdf8',
        weight: 2.5,
        dashArray: '6, 6',
      });
      this.measureLayer.addLayer(finalLine);

      if (!this.measureTooltip && this.map) {
        this.measureTooltip = L.tooltip({
          permanent: true,
          className: 'measure-tooltip',
          direction: 'right',
          offset: [10, 0],
        })
          .setLatLng(latlng)
          .addTo(this.map);
      }

      if (this.measureTooltip) {
        this.measureTooltip.setLatLng(latlng).setContent(`${t('map.distanceLabel')}: ${distStr}`);
      }

      this.measureState = 'completed';
      this.measureStartLatLng = null;
      if (this.measureCurrentLine) {
        this.measureLayer.removeLayer(this.measureCurrentLine);
        this.measureCurrentLine = null;
      }
    }
  }

  private updateMeasuringLine(currentLatLng: L.LatLng): void {
    if (!this.map || !this.measureStartLatLng) return;

    const pts = [this.measureStartLatLng, currentLatLng];
    if (!this.measureCurrentLine) {
      this.measureCurrentLine = L.polyline(pts, {
        color: '#38bdf8',
        weight: 2,
        dashArray: '4, 4',
      });
      this.measureLayer.addLayer(this.measureCurrentLine);
    } else {
      this.measureCurrentLine.setLatLngs(pts);
    }

    const dist = this.measureStartLatLng.distanceTo(currentLatLng);
    const isRo = getLanguage() === 'ro';
    const distStr = dist < 1000
      ? `${isRo ? dist.toFixed(1).replace('.', ',') : dist.toFixed(1)} m`
      : `${isRo ? (dist / 1000).toFixed(3).replace('.', ',') : (dist / 1000).toFixed(3)} km`;

    if (!this.measureTooltip) {
      this.measureTooltip = L.tooltip({
        permanent: true,
        className: 'measure-tooltip',
        direction: 'right',
        offset: [10, 0],
      })
        .setLatLng(currentLatLng)
        .setContent(`${t('map.distanceLabel')}: ${distStr}`)
        .addTo(this.map);
    } else {
      this.measureTooltip.setLatLng(currentLatLng).setContent(`${t('map.distanceLabel')}: ${distStr}`);
    }
  }

  private clearMeasurement(): void {
    this.measureLayer.clearLayers();
    this.measureStartLatLng = null;
    this.measureCurrentLine = null;
    if (this.measureTooltip && this.map) {
      this.map.removeLayer(this.measureTooltip);
      this.measureTooltip = null;
    }
  }

  public updateSinglePoint(coords: PointTransformedEvent): void {
    if (!this.map || !this.pointLayer) return;

    this.pointLayer.clearLayers();
    this.singlePointMarker = null;

    const { lat, lon, outOfGrid, pointId, stereoN, stereoE, elev, altitude, ellipsoidalHeight, zeta } = coords;

    if (isNaN(lat) || isNaN(lon) || lat < -90 || lat > 90 || lon < -180 || lon > 180) {
      this.renderLayersDock();
      return;
    }

    this.applyLayerColors('point');

    const customIcon = L.divIcon({
      className: `geodetic-marker ${outOfGrid ? 'out-of-grid' : ''}`,
      html: `
        <div class="geodetic-pin-pulse"></div>
        <div class="geodetic-pin-dot"></div>
      `,
      iconSize: [24, 24],
      iconAnchor: [12, 12],
    });

    this.singlePointMarker = L.marker([lat, lon], { icon: customIcon });

    const popupHtml = this.buildPopupHtml({
      id: pointId || 'Single Point',
      layerTitle: t('map.singlePointLayer'),
      lat,
      lon,
      stereoN,
      stereoE,
      elev,
      altitude,
      ellipsoidalHeight,
      zeta,
      warning: outOfGrid ? t('quality.outsideGrid') : null,
      angularUnit: coords.angularUnit || this.currentAngularUnit,
    });

    this.singlePointMarker.bindPopup(popupHtml);
    this.pointLayer.addLayer(this.singlePointMarker);

    this.renderLayersDock();
  }

  public clearSinglePoint(): void {
    this.pointLayer?.clearLayers();
    this.singlePointMarker = null;
    this.renderLayersDock();
  }

  public centerOnPoint(lat: number, lon: number): boolean {
    if (!this.map) return false;
    this.map.panTo([lat, lon], { animate: true });
    return true;
  }

  public zoomToPoint(lat: number, lon: number, zoom: number = 14): boolean {
    if (!this.map) return false;
    this.map.setView([lat, lon], zoom, { animate: true });
    return true;
  }

  public centerOnBatch(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    return this.panToPointsCenter(points);
  }

  public plotBatchPoints(points: MapPointItem[]): void {
    if (!this.map || !this.textLayer) return;

    this.lastTextPoints = points;
    this.renderLayersDock();

    const count = points.length;

    if (count > APP_CONFIG.batchWarningThreshold) {
      this.callbacks.onOpenVolumeWarning(count, (choice) => {
        if (choice === 'skip') {
          return;
        } else if (choice === 'sample') {
          this.renderPointsIntoLayer(this.textLayer, points.slice(0, APP_CONFIG.batchSampleSize), this.layerColors.text, false);
        } else {
          this.renderPointsIntoLayer(this.textLayer, points, this.layerColors.text, false);
        }
      });
    } else {
      this.renderPointsIntoLayer(this.textLayer, points, this.layerColors.text, false);
    }
  }

  public clearBatchPoints(): void {
    this.textLayer?.clearLayers();
    this.lastTextPoints = [];
    this.renderLayersDock();
  }

  public zoomToTextExtents(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    return this.fitPointsBounds(points);
  }

  public plotFilePoints(points: MapPointItem[]): void {
    if (!this.map || !this.filesLayer) return;

    this.lastFilePoints = points;
    this.renderLayersDock();

    const count = points.length;

    if (count > APP_CONFIG.batchWarningThreshold) {
      this.callbacks.onOpenVolumeWarning(count, (choice) => {
        if (choice === 'skip') {
          return;
        } else if (choice === 'sample') {
          this.renderPointsIntoLayer(this.filesLayer, points.slice(0, APP_CONFIG.batchSampleSize), this.layerColors.files, true);
        } else {
          this.renderPointsIntoLayer(this.filesLayer, points, this.layerColors.files, true);
        }
      });
    } else {
      this.renderPointsIntoLayer(this.filesLayer, points, this.layerColors.files, true);
    }
  }

  public clearFilePoints(): void {
    this.filesLayer?.clearLayers();
    this.lastFilePoints = [];
    this.renderLayersDock();
  }

  public centerOnFiles(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    return this.panToPointsCenter(points);
  }

  public zoomToFilesExtents(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    return this.fitPointsBounds(points);
  }

  private renderPointsIntoLayer(
    targetLayer: any,
    points: MapPointItem[],
    fillColor: string,
    isFromFile: boolean
  ): void {
    if (!targetLayer) return;

    targetLayer.clearLayers();

    const markers: L.Layer[] = [];

    for (const pt of points) {
      if (isNaN(pt.lat) || isNaN(pt.lon)) continue;

      const isOutOfGrid = pt.warning === 'Out of grid';
      const marker = L.circleMarker([pt.lat, pt.lon], {
        radius: 6,
        fillColor: isOutOfGrid ? '#ef4444' : fillColor,
        color: '#ffffff',
        weight: 1.5,
        opacity: 1,
        fillOpacity: 0.85,
      });

      const popupHtml = this.buildPopupHtml({
        id: pt.id || 'Point',
        layerTitle: isFromFile ? t('map.filesLayer') : t('map.textLayer'),
        fileName: pt.fileName,
        lat: pt.lat,
        lon: pt.lon,
        stereoN: pt.stereoN,
        stereoE: pt.stereoE,
        elev: pt.elev,
        altitude: pt.altitude,
        ellipsoidalHeight: pt.ellipsoidalHeight,
        zeta: pt.zeta,
        warning: pt.warning,
        angularUnit: pt.angularUnit || this.currentAngularUnit,
      });

      marker.bindPopup(popupHtml);

      // Subtle CAD point label for Zoom >= 15
      marker.bindTooltip(pt.id || 'Point', {
        permanent: true,
        direction: 'top',
        className: 'point-label',
        offset: [0, -6],
      });

      markers.push(marker);
    }

    targetLayer.addLayers(markers);
  }

  // Restructured popup with unified coordinate labels, File attribute, bottom grid status, and safe inline copy
  private buildPopupHtml(data: {
    id: string;
    layerTitle: string;
    fileName?: string;
    lat: number;
    lon: number;
    stereoN?: number;
    stereoE?: number;
    elev?: number;
    altitude?: number;
    ellipsoidalHeight?: number;
    zeta?: number | null;
    warning?: string | null;
    angularUnit?: AngularUnitType;
  }): string {
    const isOutOfGrid = !!data.warning && data.warning.toLowerCase().includes('grid');
    const statusVal = isOutOfGrid ? t('quality.outsideGrid') : t('quality.inGrid');

    const unit = data.angularUnit || this.currentAngularUnit;
    let latStr: string;
    let lonStr: string;
    if (unit === 'dms') {
      latStr = formatDegToDMS(data.lat);
      lonStr = formatDegToDMS(data.lon);
    } else if (unit === 'radians') {
      latStr = degToRad(data.lat).toFixed(10);
      lonStr = degToRad(data.lon).toFixed(10);
    } else {
      latStr = data.lat.toFixed(8);
      lonStr = data.lon.toFixed(8);
    }

    const copyPairs: string[] = [
      `${t('map.pointNameNumber')}: ${data.id}`,
    ];
    if (data.fileName) {
      copyPairs.push(`${t('map.sourceFile')}: ${data.fileName}`);
    }
    if (data.stereoN !== undefined) {
      copyPairs.push(`${t('map.northing')}: ${data.stereoN.toFixed(3)}`);
    }
    if (data.stereoE !== undefined) {
      copyPairs.push(`${t('map.easting')}: ${data.stereoE.toFixed(3)}`);
    }
    if (data.altitude !== undefined) {
      copyPairs.push(`${t('map.altitude')}: ${data.altitude.toFixed(3)}`);
    }
    copyPairs.push(`${t('map.latitude')}: ${latStr}`);
    copyPairs.push(`${t('map.longitude')}: ${lonStr}`);
    if (data.ellipsoidalHeight !== undefined) {
      copyPairs.push(`${t('map.height')}: ${data.ellipsoidalHeight.toFixed(3)}`);
    }
    if (data.zeta !== undefined && data.zeta !== null) {
      copyPairs.push(`${t('map.quasigeoidHeight')}: ${data.zeta.toFixed(3)}`);
    }
    copyPairs.push(`${t('map.gridCoverage')}: ${statusVal}`);

    const copyText = copyPairs.join(', ');

    return `
      <div class="geodetic-popup">
        <!-- 1. Header: Layer Name -->
        <div class="geodetic-popup-header">
          <span class="geodetic-popup-layer-title">${data.layerTitle}</span>
        </div>

        <!-- 2. Attributes Table with Unified Names (No Units) -->
        <table class="geodetic-popup-table">
          <tr>
            <td class="label">${t('map.pointNameNumber')}</td>
            <td class="val">${data.id}</td>
          </tr>
          ${
            data.fileName
              ? `
            <tr>
              <td class="label">${t('map.sourceFile')}</td>
              <td class="val" title="${data.fileName}">${data.fileName}</td>
            </tr>
          `
              : ''
          }
          <tr>
            <td class="label">${t('map.northing')}</td>
            <td class="val">${data.stereoN !== undefined ? data.stereoN.toFixed(3) : '—'}</td>
          </tr>
          <tr>
            <td class="label">${t('map.easting')}</td>
            <td class="val">${data.stereoE !== undefined ? data.stereoE.toFixed(3) : '—'}</td>
          </tr>
          <tr>
            <td class="label">${t('map.altitude')}</td>
            <td class="val">${data.altitude !== undefined ? data.altitude.toFixed(3) : '—'}</td>
          </tr>
          <tr>
            <td class="label">${t('map.latitude')}</td>
            <td class="val">${latStr}</td>
          </tr>
          <tr>
            <td class="label">${t('map.longitude')}</td>
            <td class="val">${lonStr}</td>
          </tr>
          <tr>
            <td class="label">${t('map.height')}</td>
            <td class="val">${data.ellipsoidalHeight !== undefined ? data.ellipsoidalHeight.toFixed(3) : '—'}</td>
          </tr>
          ${
            data.zeta !== undefined && data.zeta !== null
              ? `
            <tr class="metric">
              <td class="label">${t('map.quasigeoidHeight')}</td>
              <td class="val">${data.zeta.toFixed(3)}</td>
            </tr>
          `
              : ''
          }
          <tr>
            <td class="label">${t('map.gridCoverage')}</td>
            <td class="val ${isOutOfGrid ? 'out-of-grid-val' : 'in-grid-val'}">${statusVal}</td>
          </tr>
        </table>

        <!-- 3. Actions: Safe In-Place Copy Feedback (No Inline onclick) -->
        <div class="geodetic-popup-actions">
          <button 
            class="geodetic-popup-copy-btn" 
            data-copy-text="${encodeURIComponent(copyText)}"
          >
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
              <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
              <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
            </svg>
            <span>${t('map.copyAttributes')}</span>
          </button>
        </div>
      </div>
    `;
  }

  private panToPointsCenter(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    if (!this.map || points.length === 0) return false;

    const validPts = points.filter((p) => !isNaN(p.lat) && !isNaN(p.lon));
    if (validPts.length === 0) return false;

    if (validPts.length === 1) {
      return this.centerOnPoint(validPts[0].lat, validPts[0].lon);
    }

    const bounds = L.latLngBounds(validPts.map((p) => [p.lat, p.lon] as [number, number]));
    if (bounds.isValid()) {
      this.map.panTo(bounds.getCenter(), { animate: true });
      return true;
    }
    return false;
  }

  private fitPointsBounds(points: Array<{ lat: number; lon: number; warning?: string | null }>): boolean {
    if (!this.map || points.length === 0) return false;

    const validPts = points.filter((p) => !isNaN(p.lat) && !isNaN(p.lon));
    if (validPts.length === 0) return false;

    if (validPts.length === 1) {
      return this.zoomToPoint(validPts[0].lat, validPts[0].lon, 14);
    }

    const bounds = L.latLngBounds(validPts.map((p) => [p.lat, p.lon] as [number, number]));
    if (bounds.isValid()) {
      this.map.fitBounds(bounds, { padding: [50, 50], maxZoom: 14 });
      return true;
    }
    return false;
  }
}
