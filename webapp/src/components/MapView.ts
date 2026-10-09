import L from 'leaflet';
import 'leaflet/dist/leaflet.css';
import 'leaflet.markercluster';
import 'leaflet.markercluster/dist/MarkerCluster.css';
import 'leaflet.markercluster/dist/MarkerCluster.Default.css';

import { t, onLanguageChange } from '../i18n';
import { apiClient } from '../utils/apiClient';

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
  private coordsBadge: HTMLElement;
  private callbacks: MapViewCallbacks;

  private map: L.Map | null = null;
  private singlePointMarker: L.Marker | null = null;
  private clusterGroup: any = null; // L.MarkerClusterGroup
  private boundaryLayer: L.Polygon | null = null;

  constructor(callbacks: MapViewCallbacks) {
    this.callbacks = callbacks;

    this.container = document.createElement('div');
    this.container.className = 'map-panel';

    this.mapElement = document.createElement('div');
    this.mapElement.id = 'map';
    this.container.appendChild(this.mapElement);

    this.coordsBadge = document.createElement('div');
    this.coordsBadge.className = 'map-coords-badge';
    this.coordsBadge.textContent = 'φ: 45.9997186°  λ: 24.9984476°';
    this.container.appendChild(this.coordsBadge);

    onLanguageChange(() => this.updateLayerLabels());
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public initMap(): void {
    if (this.map) return;

    // Center of Romania: 45.9432° N, 24.9668° E
    this.map = L.map(this.mapElement, {
      center: [45.9432, 24.9668],
      zoom: 7,
      minZoom: 6,
      maxZoom: 18,
    });

    // Basemaps
    const cartoDark = L.tileLayer(
      'https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
      {
        attribution: '&copy; <a href="https://carto.com/">CARTO</a>',
        subdomains: 'abcd',
        maxZoom: 19,
      }
    );

    const osm = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; OpenStreetMap contributors',
      maxZoom: 19,
    });

    const satellite = L.tileLayer(
      'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
      {
        attribution: 'Tiles &copy; Esri',
        maxZoom: 18,
      }
    );

    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    if (isDark) {
      cartoDark.addTo(this.map);
    } else {
      osm.addTo(this.map);
    }

    const baseMaps = {
      [t('map.baseDark')]: cartoDark,
      [t('map.baseOSM')]: osm,
      [t('map.baseSatellite')]: satellite,
    };

    // Placeholder layers for future raster grid shifts
    const shiftEastPlaceholder = L.layerGroup();
    const shiftNorthPlaceholder = L.layerGroup();
    const shiftVectorPlaceholder = L.layerGroup();
    const shiftHeightPlaceholder = L.layerGroup();

    const overlayMaps: Record<string, L.Layer> = {
      [`${t('map.shiftsEast')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: shiftEastPlaceholder,
      [`${t('map.shiftsNorth')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: shiftNorthPlaceholder,
      [`${t('map.shiftsVector')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: shiftVectorPlaceholder,
      [`${t('map.shiftsHeight')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: shiftHeightPlaceholder,
    };

    L.control.layers(baseMaps, overlayMaps, { position: 'topright' }).addTo(this.map);

    // Initialize Marker Cluster Group
    this.clusterGroup = (L as any).markerClusterGroup({
      showCoverageOnHover: false,
      maxClusterRadius: 50,
      spiderfyOnMaxZoom: true,
    });
    this.map.addLayer(this.clusterGroup);

    // Map Event Listeners
    this.map.on('mousemove', (e: L.LeafletMouseEvent) => {
      this.coordsBadge.textContent = `φ: ${e.latlng.lat.toFixed(6)}°  λ: ${e.latlng.lng.toFixed(6)}°`;
    });

    this.map.on('click', (e: L.LeafletMouseEvent) => {
      this.callbacks.onMapPointClicked(e.latlng.lat, e.latlng.lng);
    });

    // Load Grid Bounds from API
    this.loadGridBoundary();
  }

  public updateSinglePoint(lat: number, lon: number, outOfGrid: boolean): void {
    if (!this.map) return;

    if (this.singlePointMarker) {
      this.map.removeLayer(this.singlePointMarker);
    }

    const customIcon = L.divIcon({
      className: `geodetic-marker ${outOfGrid ? 'out-of-grid' : ''}`,
      html: `
        <div class="geodetic-pin-pulse"></div>
        <div class="geodetic-pin-dot"></div>
      `,
      iconSize: [24, 24],
      iconAnchor: [12, 12],
    });

    this.singlePointMarker = L.marker([lat, lon], { icon: customIcon }).addTo(this.map);

    // Smoothly pan map to new coordinates
    this.map.panTo([lat, lon], { animate: true });
  }

  public plotBatchPoints(
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ): void {
    if (!this.map || !this.clusterGroup) return;

    const count = points.length;

    // High volume pre-warning threshold (> 2,000 points)
    if (count > 2000) {
      this.callbacks.onOpenVolumeWarning(count, (choice) => {
        if (choice === 'skip') {
          return;
        } else if (choice === 'sample') {
          this.renderBatchPoints(points.slice(0, 1000));
        } else {
          this.renderBatchPoints(points);
        }
      });
    } else {
      this.renderBatchPoints(points);
    }
  }

  private renderBatchPoints(
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ): void {
    if (!this.clusterGroup || !this.map) return;

    this.clusterGroup.clearLayers();

    const markers: L.Layer[] = [];
    const bounds = L.latLngBounds([]);

    for (const pt of points) {
      if (isNaN(pt.lat) || isNaN(pt.lon)) continue;

      const isOutOfGrid = pt.warning === 'Out of grid';
      const marker = L.circleMarker([pt.lat, pt.lon], {
        radius: 6,
        fillColor: isOutOfGrid ? '#ef4444' : '#0284c7',
        color: '#ffffff',
        weight: 1.5,
        opacity: 1,
        fillOpacity: 0.85,
      });

      let popupContent = `<strong>${pt.id || 'Point'}</strong><br/>φ: ${pt.lat.toFixed(6)}°<br/>λ: ${pt.lon.toFixed(6)}°`;
      if (pt.elev !== undefined) {
        popupContent += `<br/>Elev: ${pt.elev.toFixed(3)} m`;
      }
      if (pt.warning) {
        popupContent += `<br/><span style="color: #ef4444; font-weight: bold;">⚠ ${pt.warning}</span>`;
      }

      marker.bindPopup(popupContent);
      markers.push(marker);
      bounds.extend([pt.lat, pt.lon]);
    }

    this.clusterGroup.addLayers(markers);

    if (markers.length > 0) {
      this.map.fitBounds(bounds, { padding: [50, 50], maxZoom: 14 });
    }
  }

  private async loadGridBoundary(): Promise<void> {
    try {
      await apiClient.getGridInfo();

      // Romania National TransDatRO Envelope:
      // Approximate geodetic perimeter bounding the national territory
      const romaniaEnvelope: [number, number][] = [
        [43.618, 28.583],
        [43.833, 25.967],
        [43.683, 24.083],
        [44.350, 22.683],
        [44.750, 21.433],
        [45.300, 20.800],
        [46.133, 20.267],
        [46.600, 21.233],
        [47.883, 22.883],
        [48.267, 23.217],
        [47.950, 25.267],
        [48.267, 26.700],
        [47.167, 28.050],
        [45.467, 28.217],
        [45.200, 29.650],
        [44.400, 28.850],
        [43.618, 28.583],
      ];

      if (this.map) {
        this.boundaryLayer = L.polygon(romaniaEnvelope, {
          color: '#38bdf8',
          weight: 2,
          dashArray: '6, 6',
          fillColor: '#0284c7',
          fillOpacity: 0.05,
        }).addTo(this.map);

        this.boundaryLayer.bindTooltip(t('map.gridBoundary'), { sticky: true });
      }
    } catch (err) {
      console.warn('Failed to fetch grid bounds from API', err);
    }
  }

  private updateLayerLabels(): void {
    // When language changes, update any dynamic labels
  }
}
