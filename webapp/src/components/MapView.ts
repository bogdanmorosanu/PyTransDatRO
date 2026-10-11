import L from 'leaflet';
import 'leaflet/dist/leaflet.css';
import 'leaflet.markercluster';
import 'leaflet.markercluster/dist/MarkerCluster.css';
import 'leaflet.markercluster/dist/MarkerCluster.Default.css';

import { t, onLanguageChange } from '../i18n';
import { apiClient } from '../utils/apiClient';

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
  private pointLayer: L.LayerGroup | null = null;
  private textLayer: any = null; // L.MarkerClusterGroup
  private filesLayer: any = null; // L.MarkerClusterGroup
  private boundaryGroup: L.LayerGroup | null = null;
  private boundaryPolygon: L.Polygon | null = null;
  private layerControl: L.Control.Layers | null = null;

  private cartoDark: L.TileLayer | null = null;
  private osm: L.TileLayer | null = null;
  private satellite: L.TileLayer | null = null;

  private shiftEastPlaceholder = L.layerGroup();
  private shiftNorthPlaceholder = L.layerGroup();
  private shiftVectorPlaceholder = L.layerGroup();
  private shiftHeightPlaceholder = L.layerGroup();

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

    onLanguageChange(() => {
      this.updateLayerControl();
      this.refreshBoundaryTooltip();
    });

    // Sync basemap with theme
    const observer = new MutationObserver((mutations) => {
      for (const m of mutations) {
        if (m.type === 'attributes' && m.attributeName === 'data-theme') {
          const isDarkTheme = document.documentElement.getAttribute('data-theme') === 'dark';
          if (this.map && this.cartoDark && this.osm) {
            if (isDarkTheme) {
              if (this.map.hasLayer(this.osm)) this.map.removeLayer(this.osm);
              if (!this.map.hasLayer(this.cartoDark) && !this.map.hasLayer(this.satellite!)) {
                this.cartoDark.addTo(this.map);
              }
            } else {
              if (this.map.hasLayer(this.cartoDark)) this.map.removeLayer(this.cartoDark);
              if (!this.map.hasLayer(this.osm) && !this.map.hasLayer(this.satellite!)) {
                this.osm.addTo(this.map);
              }
            }
          }
        }
      }
    });
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
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
    this.cartoDark = L.tileLayer(
      'https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
      {
        attribution: '&copy; <a href="https://carto.com/">CARTO</a>',
        subdomains: 'abcd',
        maxZoom: 19,
      }
    );

    this.osm = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; OpenStreetMap contributors',
      maxZoom: 19,
    });

    this.satellite = L.tileLayer(
      'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
      {
        attribution: 'Tiles &copy; Esri',
        maxZoom: 18,
      }
    );

    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    if (isDark) {
      this.cartoDark.addTo(this.map);
    } else {
      this.osm.addTo(this.map);
    }

    // Initialize 3 independent data layers + boundary layer
    this.pointLayer = L.layerGroup().addTo(this.map);
    this.textLayer = (L as any).markerClusterGroup({
      showCoverageOnHover: false,
      maxClusterRadius: 50,
      spiderfyOnMaxZoom: true,
    }).addTo(this.map);
    this.filesLayer = (L as any).markerClusterGroup({
      showCoverageOnHover: false,
      maxClusterRadius: 50,
      spiderfyOnMaxZoom: true,
    }).addTo(this.map);
    this.boundaryGroup = L.layerGroup().addTo(this.map);

    this.updateLayerControl();

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

  private updateLayerControl(): void {
    if (
      !this.map ||
      !this.cartoDark ||
      !this.osm ||
      !this.satellite ||
      !this.pointLayer ||
      !this.textLayer ||
      !this.filesLayer ||
      !this.boundaryGroup
    ) {
      return;
    }

    if (this.layerControl) {
      this.map.removeControl(this.layerControl);
    }

    const baseMaps: Record<string, L.Layer> = {
      [t('map.baseDark')]: this.cartoDark,
      [t('map.baseOSM')]: this.osm,
      [t('map.baseSatellite')]: this.satellite,
    };

    const overlayMaps: Record<string, L.Layer> = {
      [t('map.singlePointLayer')]: this.pointLayer,
      [t('map.textLayer')]: this.textLayer,
      [t('map.filesLayer')]: this.filesLayer,
      [t('map.gridBoundary')]: this.boundaryGroup,
      [`${t('map.shiftsEast')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: this.shiftEastPlaceholder,
      [`${t('map.shiftsNorth')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: this.shiftNorthPlaceholder,
      [`${t('map.shiftsVector')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: this.shiftVectorPlaceholder,
      [`${t('map.shiftsHeight')} <span class="layer-badge-soon">${t('map.comingSoon')}</span>`]: this.shiftHeightPlaceholder,
    };

    this.layerControl = L.control.layers(baseMaps, overlayMaps, { position: 'topright' }).addTo(this.map);
  }

  private refreshBoundaryTooltip(): void {
    if (this.boundaryPolygon) {
      this.boundaryPolygon.unbindTooltip();
      this.boundaryPolygon.bindTooltip(t('map.gridBoundary'), { sticky: true });
    }
  }

  public updateSinglePoint(lat: number, lon: number, outOfGrid: boolean): void {
    if (!this.map || !this.pointLayer) return;

    this.pointLayer.clearLayers();
    this.singlePointMarker = null;

    if (isNaN(lat) || isNaN(lon) || lat < -90 || lat > 90 || lon < -180 || lon > 180) return;

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
    this.pointLayer.addLayer(this.singlePointMarker);

    // Decoupled non-panning update: do NOT panTo automatically.
  }

  public clearSinglePoint(): void {
    this.pointLayer?.clearLayers();
    this.singlePointMarker = null;
  }

  public centerOnPoint(lat: number, lon: number): boolean {
    if (!this.map) return false;
    if (!isInsideRomaniaGrid(lat, lon)) {
      return false;
    }
    this.map.panTo([lat, lon], { animate: true });
    return true;
  }

  public zoomToPoint(lat: number, lon: number, zoom: number = 14): boolean {
    if (!this.map) return false;
    if (!isInsideRomaniaGrid(lat, lon)) {
      return false;
    }
    this.map.setView([lat, lon], zoom, { animate: true });
    return true;
  }

  public centerOnBatch(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    return this.panToPointsCenter(points);
  }

  public plotBatchPoints(
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ): void {
    if (!this.map || !this.textLayer) return;

    const count = points.length;

    if (count > 2000) {
      this.callbacks.onOpenVolumeWarning(count, (choice) => {
        if (choice === 'skip') {
          return;
        } else if (choice === 'sample') {
          this.renderPointsIntoLayer(this.textLayer, points.slice(0, 1000));
        } else {
          this.renderPointsIntoLayer(this.textLayer, points);
        }
      });
    } else {
      this.renderPointsIntoLayer(this.textLayer, points);
    }
  }

  public clearBatchPoints(): void {
    this.textLayer?.clearLayers();
  }

  public zoomToTextExtents(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    return this.fitPointsBounds(points);
  }

  public plotFilePoints(
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ): void {
    if (!this.map || !this.filesLayer) return;

    const count = points.length;

    if (count > 2000) {
      this.callbacks.onOpenVolumeWarning(count, (choice) => {
        if (choice === 'skip') {
          return;
        } else if (choice === 'sample') {
          this.renderPointsIntoLayer(this.filesLayer, points.slice(0, 1000));
        } else {
          this.renderPointsIntoLayer(this.filesLayer, points);
        }
      });
    } else {
      this.renderPointsIntoLayer(this.filesLayer, points);
    }
  }

  public clearFilePoints(): void {
    this.filesLayer?.clearLayers();
  }

  public centerOnFiles(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    return this.panToPointsCenter(points);
  }

  public zoomToFilesExtents(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    return this.fitPointsBounds(points);
  }

  private renderPointsIntoLayer(
    targetLayer: any,
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ): void {
    if (!targetLayer) return;

    targetLayer.clearLayers();

    const markers: L.Layer[] = [];

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
    }

    targetLayer.addLayers(markers);
    // Note: Decoupled non-panning: do NOT call fitBounds automatically.
  }

  private panToPointsCenter(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    if (!this.map || points.length === 0) return false;

    const validInGrid = points.filter(
      (p) => (!p.warning || p.warning !== 'Out of grid') && isInsideRomaniaGrid(p.lat, p.lon)
    );

    if (validInGrid.length === 0) return false;

    if (validInGrid.length === 1) {
      return this.centerOnPoint(validInGrid[0].lat, validInGrid[0].lon);
    }

    const bounds = L.latLngBounds(validInGrid.map((p) => [p.lat, p.lon] as [number, number]));
    if (bounds.isValid()) {
      this.map.panTo(bounds.getCenter(), { animate: true });
      return true;
    }
    return false;
  }

  private fitPointsBounds(
    points: Array<{ lat: number; lon: number; warning?: string | null }>
  ): boolean {
    if (!this.map || points.length === 0) return false;

    const validInGrid = points.filter(
      (p) => (!p.warning || p.warning !== 'Out of grid') && isInsideRomaniaGrid(p.lat, p.lon)
    );

    if (validInGrid.length === 0) return false;

    if (validInGrid.length === 1) {
      return this.zoomToPoint(validInGrid[0].lat, validInGrid[0].lon, 14);
    }

    const bounds = L.latLngBounds(validInGrid.map((p) => [p.lat, p.lon] as [number, number]));
    if (bounds.isValid()) {
      this.map.fitBounds(bounds, { padding: [50, 50], maxZoom: 14 });
      return true;
    }
    return false;
  }

  private async loadGridBoundary(): Promise<void> {
    try {
      await apiClient.getGridInfo();

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

      if (this.boundaryGroup) {
        this.boundaryGroup.clearLayers();
        this.boundaryPolygon = L.polygon(romaniaEnvelope, {
          color: '#38bdf8',
          weight: 2,
          dashArray: '6, 6',
          fillColor: '#0284c7',
          fillOpacity: 0.05,
        });
        this.boundaryPolygon.bindTooltip(t('map.gridBoundary'), { sticky: true });
        this.boundaryGroup.addLayer(this.boundaryPolygon);
      }
    } catch (err) {
      console.warn('Failed to fetch grid bounds from API', err);
    }
  }
}
