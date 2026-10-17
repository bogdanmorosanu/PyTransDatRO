import { t, onLanguageChange } from '../i18n';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { parseInputLines, detectDelimiter, type CoordinateOrder } from '../utils/pointParser';
import { formatDegToDMS, degToRad } from '../utils/dmsParser';
import { exportToCSV, exportToGeoJSON, exportToDXF, type ExportPointItem } from '../utils/exporters';
import { showToast } from '../utils/toast';
import { renderCopyButton } from '../ui';

export type AngularUnitType = 'degrees' | 'dms' | 'radians';

export interface BatchReportData {
  total: number;
  valid: number;
  outOfGrid: number;
  invalid: number;
  items: ExportPointItem[];
}

export interface BatchMapPoint {
  id: string | null;
  lat: number;
  lon: number;
  elev?: number;
  altitude?: number;
  ellipsoidalHeight?: number;
  zeta?: number | null;
  stereoN?: number;
  stereoE?: number;
  warning?: string | null;
  angularUnit?: AngularUnitType;
}

export interface BatchConverterCallbacks {
  onPointsPlotted: (points: BatchMapPoint[]) => void;
  onCenterOnMap?: (points: BatchMapPoint[]) => void;
  onZoomToExtents?: (points: BatchMapPoint[]) => void;
}

export class BatchConverter {
  private container: HTMLElement;
  private callbacks: BatchConverterCallbacks;

  // Settings inherited from shared config
  private op: TransformationOp = 'ETRS89ToStereo70';
  private angularUnit: AngularUnitType = 'degrees';
  private order: CoordinateOrder = 'NE';
  private delimiter: string = 'auto';

  // State
  private inputText: string = '';
  private outputText: string = '';
  private lastTransformedItems: ExportPointItem[] = [];
  private lastMapPoints: BatchMapPoint[] = [];
  private lastSummary: {
    inGridCount: number;
    outOfGridCount: number;
    syntaxErrorCount: number;
  } | null = null;
  private isProcessing: boolean = false;

  constructor(callbacks: BatchConverterCallbacks) {
    this.callbacks = callbacks;
    this.container = document.createElement('div');
    this.container.className = 'tab-content';

    // Sample default points for immediate exploration (ETRS89 Romania)
    this.inputText = [
      '101, 46.00000000, 25.00000000, 100.000',
      '102, 46.00225000, 25.00195000, 102.500',
      '103, 46.00900000, 24.99350000, 98.300',
      '104, 10.00000000, 10.00000000, 0.000', // Intentional out of grid sample
    ].join('\n');

    this.render();
    onLanguageChange(() => this.render());
  }

  private getDynamicPlaceholder(): string {
    const isEtrs = this.op === 'ETRS89ToStereo70';
    if (isEtrs) {
      if (this.angularUnit === 'dms') {
        return this.order === 'EN'
          ? '101, 25° 00\' 00.00000", 46° 00\' 00.00000", 100.000\n102, 25° 01\' 15.20000", 46° 02\' 30.10000", 102.500'
          : '101, 46° 00\' 00.00000", 25° 00\' 00.00000", 100.000\n102, 46° 02\' 30.10000", 25° 01\' 15.20000", 102.500';
      }
      return this.order === 'EN'
        ? '101, 25.00000000, 46.00000000, 100.000\n102, 25.02000000, 46.05000000, 102.500'
        : '101, 46.00000000, 25.00000000, 100.000\n102, 46.05000000, 25.02000000, 102.500';
    }
    return this.order === 'EN'
      ? '101, 500000.000, 500000.000, 100.000\n102, 500150.000, 500200.000, 102.500'
      : '101, 500000.000, 500000.000, 100.000\n102, 500200.000, 500150.000, 102.500';
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public setDirection(op: TransformationOp): void {
    if (this.op === op) return;
    this.op = op;

    const defaultEtrsSample = [
      '101, 46.00000000, 25.00000000, 100.000',
      '102, 46.00225000, 25.00195000, 102.500',
      '103, 46.00900000, 24.99350000, 98.300',
      '104, 10.00000000, 10.00000000, 0.000',
    ].join('\n');

    const defaultS70Sample = [
      '101, 500000.000, 500000.000, 100.000',
      '102, 500250.000, 500200.000, 102.500',
      '103, 501000.000, 499300.000, 98.300',
      '104, 100000.000, 100000.000, 0.000',
    ].join('\n');

    if (this.op === 'Stereo70ToETRS89' && this.inputText.trim() === defaultEtrsSample.trim()) {
      this.inputText = defaultS70Sample;
    } else if (this.op === 'ETRS89ToStereo70' && this.inputText.trim() === defaultS70Sample.trim()) {
      this.inputText = defaultEtrsSample;
    }

    this.render();
  }

  public setAngularUnit(unit: AngularUnitType): void {
    if (this.angularUnit === unit) return;
    this.angularUnit = unit;
    this.render();
  }

  public setCoordinateOrder(order: CoordinateOrder): void {
    if (this.order === order) return;
    this.order = order;
    this.render();
  }

  public getCoordinateOrder(): CoordinateOrder {
    return this.order;
  }


  private async processBatch(): Promise<void> {
    if (this.isProcessing) return;
    if (!this.inputText.trim()) {
      showToast(t('batch.warnNoRows'), 'warning');
      return;
    }

    this.isProcessing = true;
    this.render();

    try {
      const isSourceGeo = this.op === 'ETRS89ToStereo70';
      const parsedRows = parseInputLines(this.inputText, this.delimiter, this.order, isSourceGeo, this.angularUnit);

      const validRows = parsedRows.filter((r) => !r.error && !r.isComment && !r.isEmpty && r.coords.length >= 2);
      const syntaxErrorRows = parsedRows.filter((r) => !r.isComment && !r.isEmpty && !!r.error);

      if (validRows.length === 0 && syntaxErrorRows.length === 0) {
        throw new Error('No valid coordinate rows could be parsed');
      }

      let batchResults: any[] = [];
      if (validRows.length > 0) {
        const payloadPoints = validRows.map((r) => r.coords);
        const res = await apiClient.transformBatch({
          op: this.op,
          points: payloadPoints,
          unit: 'degrees',
        });
        batchResults = res.results;
      }

      // Assemble results strictly maintaining 1-to-1 input-to-output sequence
      const items: ExportPointItem[] = [];
      const outputLines: string[] = [];
      const mapPoints: BatchMapPoint[] = [];

      let validCount = 0;
      let outOfGridCount = 0;
      let validRowIndex = 0;
      let pointIndex = 0;

      for (const row of parsedRows) {
        if (row.isEmpty) {
          outputLines.push('');
          continue;
        }
        if (row.isComment) {
          outputLines.push(row.rawLine);
          continue;
        }
        if (row.error) {
          outputLines.push(`${row.rawLine} # Error: ${row.error}`);
          continue;
        }

        // Process valid coordinate row
        pointIndex++;
        const pointLabel = row.pointId || `Point ${pointIndex}`;
        const apiRes = batchResults[validRowIndex++];
        const isGeoOutput = this.op === 'Stereo70ToETRS89';
        const isOutOfGrid = apiRes.warning === 'Out of grid';

        if (isOutOfGrid) {
          outOfGridCount++;
          outputLines.push(row.pointId ? `${row.pointId}, # Out of grid` : '# Out of grid');
          items.push({
            pointId: row.pointId || pointLabel,
            coords: [],
            warning: 'Out of grid',
            isGeo: isGeoOutput,
          });
          continue;
        }

        validCount++;

        const item: ExportPointItem = {
          pointId: pointLabel,
          coords: apiRes.coos,
          warning: apiRes.warning,
          isGeo: isGeoOutput,
        };
        items.push(item);

        const parts: string[] = [];
        if (row.pointId) parts.push(row.pointId);

        if (isGeoOutput) {
          const lat = apiRes.coos[0];
          const lon = apiRes.coos[1];
          const elev = apiRes.coos[2];

          let latStr: string;
          let lonStr: string;
          if (this.angularUnit === 'dms') {
            latStr = formatDegToDMS(lat);
            lonStr = formatDegToDMS(lon);
          } else if (this.angularUnit === 'radians') {
            latStr = degToRad(lat).toFixed(10);
            lonStr = degToRad(lon).toFixed(10);
          } else {
            latStr = lat.toFixed(8);
            lonStr = lon.toFixed(8);
          }

          if (this.order === 'EN') {
            parts.push(lonStr, latStr);
          } else {
            parts.push(latStr, lonStr);
          }

          if (elev !== undefined) parts.push(elev.toFixed(3));

          const stereoN = this.order === 'EN' ? row.coords[1] : row.coords[0];
          const stereoE = this.order === 'EN' ? row.coords[0] : row.coords[1];
          const alt = row.coords[2];
          const h = elev;
          const zeta = (h !== undefined && alt !== undefined) ? (h - alt) : undefined;
          mapPoints.push({
            id: pointLabel,
            lat,
            lon,
            elev,
            altitude: alt,
            ellipsoidalHeight: h,
            zeta,
            stereoN,
            stereoE,
            warning: apiRes.warning,
            angularUnit: this.angularUnit,
          });
        } else {
          const n = apiRes.coos[0];
          const e = apiRes.coos[1];
          const z = apiRes.coos[2];

          if (this.order === 'EN') {
            parts.push(e.toFixed(3), n.toFixed(3));
          } else {
            parts.push(n.toFixed(3), e.toFixed(3));
          }

          if (z !== undefined) parts.push(z.toFixed(3));

          const origLat = row.coords[0];
          const origLon = row.coords[1];
          const h = row.coords[2];
          const alt = z;
          const zeta = (h !== undefined && alt !== undefined) ? (h - alt) : undefined;
          mapPoints.push({
            id: pointLabel,
            lat: origLat,
            lon: origLon,
            elev: z,
            altitude: alt,
            ellipsoidalHeight: h,
            zeta,
            stereoN: n,
            stereoE: e,
            warning: apiRes.warning,
            angularUnit: this.angularUnit,
          });
        }

        if (apiRes.warning) {
          parts.push(`# ${apiRes.warning}`);
        }

        const detectedDelim = this.delimiter === 'auto' ? detectDelimiter(this.inputText.split(/\r?\n/)) : this.delimiter;
        const outDelim = detectedDelim === ' ' ? ', ' : detectedDelim === '\t' ? '\t' : `${detectedDelim} `;
        outputLines.push(parts.join(outDelim));
      }

      this.outputText = outputLines.join('\n');
      this.lastTransformedItems = items;
      this.lastMapPoints = mapPoints;
      this.lastSummary = {
        inGridCount: validCount,
        outOfGridCount: outOfGridCount,
        syntaxErrorCount: syntaxErrorRows.length,
      };

      this.callbacks.onPointsPlotted(mapPoints);

      showToast(t('batch.successToast', { count: validCount }), 'success');
    } catch (err: any) {
      showToast(err.message || 'Error processing batch', 'danger');
    } finally {
      this.isProcessing = false;
      this.render();
    }
  }

  private render(): void {
    this.container.innerHTML = `
      <div class="studio-card">
        <!-- Vertical Stack: Input Textarea with Direction-Aware Dynamic Placeholder -->
        <div class="form-group" style="margin-bottom: 12px;">
          <label class="form-label">${t('batch.inputLabel')}</label>
          <textarea 
            class="form-input mono" 
            id="batch-input-area" 
            aria-label="${t('batch.inputLabel')}"
            rows="6" 
            placeholder="${this.getDynamicPlaceholder()}"
            style="resize: vertical; font-size: 12px; line-height: 1.4;"
          >${this.inputText}</textarea>
        </div>

        <!-- Transform Primary Button -->
        <button class="btn btn-primary" id="batch-transform-btn" style="width: 100%; margin-bottom: ${this.outputText ? '16px' : '0'};">
          <span>${this.isProcessing ? t('batch.calculating') : t('batch.transformBtn')}</span>
        </button>

        <!-- Output Results Area -->
        ${
          this.outputText
            ? `
          <div class="form-group" style="margin-bottom: 10px;">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 6px;">
              <label class="form-label" style="margin-bottom: 0;">${t('batch.outputLabel')}</label>
              <div id="batch-copy-results-slot"></div>
            </div>
            <textarea 
              class="form-input mono" 
              id="batch-output-area" 
              aria-label="${t('batch.outputLabel')}"
              rows="6" 
              readonly 
              style="resize: vertical; font-size: 12px; line-height: 1.4; background-color: var(--bg-surface-elevated);"
            >${this.outputText}</textarea>
          </div>

          <!-- Standardized Quality Report Widget -->
          ${
            this.lastSummary
              ? `
            <div class="quality-report-card" style="margin-bottom: 12px;">
              <div class="quality-report-row">
                <div class="quality-report-badges">
                  <span class="status-pill success">✓ ${t('quality.inGrid')}: ${this.lastSummary.inGridCount}</span>
                  ${
                    this.lastSummary.outOfGridCount > 0
                      ? `<span class="status-pill warning">⚠ ${t('quality.outsideGrid')}: ${this.lastSummary.outOfGridCount}</span>`
                      : ''
                  }
                  ${
                    this.lastSummary.syntaxErrorCount > 0
                      ? `<span class="status-pill danger">✕ ${t('quality.syntaxErrors', { count: this.lastSummary.syntaxErrorCount })}</span>`
                      : ''
                  }
                </div>
              </div>
            </div>
          `
              : ''
          }

          <!-- Relocated Export Buttons below Output -->
          <div class="download-group">
            <span class="download-group-label">${t('batch.downloadAs')}</span>
            <div style="display: flex; gap: 6px;">
              <button class="btn btn-secondary" id="batch-export-csv" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''} style="padding: 6px 12px; font-size: 12px;">
                ${t('batch.downloadCSV')}
              </button>
              <button class="btn btn-secondary" id="batch-export-geojson" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''} style="padding: 6px 12px; font-size: 12px;">
                ${t('batch.downloadGeoJSON')}
              </button>
              <button class="btn btn-secondary" id="batch-export-dxf" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''} style="padding: 6px 12px; font-size: 12px;">
                ${t('batch.downloadDXF')}
              </button>
            </div>
          </div>

          <!-- Decoupled Map Navigation Actions on its own Bottom-Right Line -->
          <div class="action-footer-right" style="display: flex; gap: 8px; justify-content: flex-end; width: 100%; margin-top: 10px;">
            <button 
              class="btn btn-outline center-map-btn" 
              id="batch-center-map-btn" 
              ${this.lastMapPoints.length === 0 ? 'disabled' : ''} 
              style="padding: 4px 10px; font-size: 11.5px;"
            >
              <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                <circle cx="12" cy="12" r="10"></circle>
                <line x1="22" y1="12" x2="18" y2="12"></line>
                <line x1="6" y1="12" x2="2" y2="12"></line>
                <line x1="12" y1="6" x2="12" y2="2"></line>
                <line x1="12" y1="22" x2="12" y2="18"></line>
              </svg>
              <span>${t('common.centerOnMap')}</span>
            </button>
            <button 
              class="btn btn-outline zoom-extents-btn" 
              id="batch-zoom-extents-btn" 
              ${this.lastMapPoints.length === 0 ? 'disabled' : ''} 
              style="padding: 4px 10px; font-size: 11.5px;"
            >
              <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                <circle cx="11" cy="11" r="8"></circle>
                <line x1="21" y1="21" x2="16.65" y2="16.65"></line>
                <line x1="11" y1="8" x2="11" y2="14"></line>
                <line x1="8" y1="11" x2="14" y2="11"></line>
              </svg>
              <span>${t('common.zoomToExtents')}</span>
            </button>
          </div>
        `
            : ''
        }
      </div>
    `;

    this.attachEventListeners();
  }

  private attachEventListeners(): void {
    // Input sync
    const inputArea = this.container.querySelector('#batch-input-area') as HTMLTextAreaElement;
    if (inputArea) {
      inputArea.addEventListener('input', () => {
        this.inputText = inputArea.value;
      });
    }

    // Top-right Copy Results action button using reusable UI factory
    const copySlot = this.container.querySelector('#batch-copy-results-slot');
    if (copySlot) {
      const copyResultsBtn = renderCopyButton({
        id: 'batch-copy-results-btn',
        title: t('batch.copyAll'),
        ariaLabel: t('batch.copyAll'),
        disabled: !this.outputText.trim(),
        feedbackText: t('point.copiedTooltip'),
        getCopyText: () => this.outputText,
        onCopied: () => {
          showToast(t('batch.copySuccess'), 'success');
        },
      });
      copySlot.replaceWith(copyResultsBtn);
    }

    // Transform button
    this.container.querySelector('#batch-transform-btn')?.addEventListener('click', () => {
      if (inputArea) this.inputText = inputArea.value;
      this.processBatch();
    });

    // Center on Map button
    this.container.querySelector('#batch-center-map-btn')?.addEventListener('click', () => {
      if (this.lastMapPoints.length > 0) {
        this.callbacks.onCenterOnMap?.(this.lastMapPoints);
      }
    });

    // Zoom to Extents button
    this.container.querySelector('#batch-zoom-extents-btn')?.addEventListener('click', () => {
      if (this.lastMapPoints.length > 0) {
        this.callbacks.onZoomToExtents?.(this.lastMapPoints);
      }
    });

    // Export buttons
    this.container.querySelector('#batch-export-csv')?.addEventListener('click', () => {
      const detectedDelim = this.delimiter === 'auto' ? detectDelimiter(this.inputText.split(/\r?\n/)) : this.delimiter;
      const csvSep = detectedDelim === ' ' ? ',' : detectedDelim;
      exportToCSV(this.lastTransformedItems, csvSep, this.order);
      showToast(t('batch.csvDownloaded'), 'success');
    });

    this.container.querySelector('#batch-export-geojson')?.addEventListener('click', () => {
      exportToGeoJSON(this.lastTransformedItems);
      showToast(t('batch.geoJsonDownloaded'), 'success');
    });

    this.container.querySelector('#batch-export-dxf')?.addEventListener('click', () => {
      exportToDXF(this.lastTransformedItems);
      showToast(t('batch.dxfDownloaded'), 'success');
    });
  }
}
