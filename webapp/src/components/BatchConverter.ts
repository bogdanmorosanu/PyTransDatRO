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
    totalCount: number;
    inGridCount: number;
    outOfGridCount: number;
    noDataCount: number;
    invalidInputCount: number;
    invalidAngleCount: number;
    count2D: number;
    count3D: number;
  } | null = null;
  private savedInputHeight?: string;
  private savedOutputHeight?: string;
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

  public clear(): void {
    this.inputText = '';
    this.outputText = '';
    this.lastTransformedItems = [];
    this.lastMapPoints = [];
    this.lastSummary = null;
    this.render();
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
      let noDataCount = 0;
      let count2D = 0;
      let count3D = 0;
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
          const commentText = row.errorType === 'invalid_angle' ? '# Invalid angle' : '# Invalid input';
          outputLines.push(row.pointId ? `${row.pointId}, ${commentText}` : commentText);
          items.push({
            pointId: row.pointId || `Point ${pointIndex + 1}`,
            coords: [],
            warning: row.errorType === 'invalid_angle' ? 'Invalid angle' : 'Invalid input',
            isGeo: this.op === 'Stereo70ToETRS89',
          });
          continue;
        }

        // Process valid coordinate row
        pointIndex++;
        const pointLabel = row.pointId || `Point ${pointIndex}`;
        const apiRes = batchResults[validRowIndex++];
        const isGeoOutput = this.op === 'Stereo70ToETRS89';
        const isUntransformable =
          apiRes.warning === 'Out of grid' || apiRes.warning === 'No data on grid';

        if (isUntransformable) {
          if (apiRes.warning === 'No data on grid') {
            noDataCount++;
          } else {
            outOfGridCount++;
          }
          const warnComment = `# ${apiRes.warning}`;
          outputLines.push(row.pointId ? `${row.pointId}, ${warnComment}` : warnComment);
          items.push({
            pointId: row.pointId || pointLabel,
            coords: [],
            warning: apiRes.warning,
            isGeo: isGeoOutput,
          });
          continue;
        }

        validCount++;
        if (row.coords.length >= 3) {
          count3D++;
        } else {
          count2D++;
        }

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
      let invalidInputCount = 0;
      let invalidAngleCount = 0;
      for (const r of syntaxErrorRows) {
        if (r.errorType === 'invalid_angle') {
          invalidAngleCount++;
        } else {
          invalidInputCount++;
        }
      }

      this.lastSummary = {
        totalCount: validRows.length + syntaxErrorRows.length,
        inGridCount: validCount,
        outOfGridCount: outOfGridCount,
        noDataCount: noDataCount,
        invalidInputCount,
        invalidAngleCount,
        count2D,
        count3D,
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

  private captureTextareaHeights(): void {
    const inputArea = this.container.querySelector('#batch-input-area') as HTMLTextAreaElement | null;
    if (inputArea && inputArea.style.height) {
      this.savedInputHeight = inputArea.style.height;
    }
    const outputArea = this.container.querySelector('#batch-output-area') as HTMLTextAreaElement | null;
    if (outputArea && outputArea.style.height) {
      this.savedOutputHeight = outputArea.style.height;
    }
  }

  private restoreTextareaHeights(): void {
    if (this.savedInputHeight) {
      const inputArea = this.container.querySelector('#batch-input-area') as HTMLTextAreaElement | null;
      if (inputArea) inputArea.style.height = this.savedInputHeight;
    }
    if (this.savedOutputHeight) {
      const outputArea = this.container.querySelector('#batch-output-area') as HTMLTextAreaElement | null;
      if (outputArea) outputArea.style.height = this.savedOutputHeight;
    }
  }

  private render(): void {
    this.captureTextareaHeights();

    this.container.innerHTML = `
      <!-- Card 1: Input Section -->
      <div class="studio-card" role="group" aria-label="${t('batch.inputLabel')}">
        <div class="section-micro-label">${t('batch.inputLabel')}</div>
        <div class="form-group margin-bottom-none">
          <textarea 
            class="form-input mono batch-textarea" 
            id="batch-input-area" 
            aria-label="${t('batch.inputLabel')}"
            rows="6" 
            placeholder="${this.getDynamicPlaceholder()}"
            ${this.isProcessing ? 'disabled readonly' : ''}
          >${this.inputText}</textarea>
        </div>
      </div>

      <!-- Active Separator Button between Input & Output -->
      <div class="transform-action-row">
        <button class="btn btn-primary" id="batch-transform-btn" ${this.isProcessing ? 'disabled' : ''}>
          <span>${this.isProcessing ? t('batch.calculating') : t('batch.transformBtn')}</span>
        </button>
      </div>

      <!-- Card 2: Output Section -->
      <div class="studio-card result-card" role="group" aria-label="${t('batch.outputLabel')}">
        <div class="flex-between-center margin-bottom-sm">
          <div class="section-micro-label">${t('batch.outputLabel')}</div>
          <div id="batch-copy-results-slot"></div>
        </div>
        <div class="form-group margin-bottom-none">
          <textarea 
            class="form-input mono batch-textarea batch-textarea-readonly" 
            id="batch-output-area" 
            aria-label="${t('batch.outputLabel')}"
            rows="6" 
            readonly 
          >${this.outputText}</textarea>
        </div>

        <!-- Export Buttons inside Output Card -->
        <div class="download-group">
          <span class="download-group-label">${t('batch.downloadAs')}</span>
          <div class="flex-row-gap-sm">
            <button class="btn btn-secondary btn-secondary-action" id="batch-export-csv" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              ${t('batch.downloadCSV')}
            </button>
            <button class="btn btn-secondary btn-secondary-action" id="batch-export-geojson" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              ${t('batch.downloadGeoJSON')}
            </button>
            <button class="btn btn-secondary btn-secondary-action" id="batch-export-dxf" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              ${t('batch.downloadDXF')}
            </button>
          </div>
        </div>
      </div>

      <!-- Card 3: Report Section -->
      <div class="studio-card report-card" role="group" aria-label="${t('quality.sectionTitle')}">
        <div class="section-micro-label">${t('quality.sectionTitle')}</div>
        <div class="quality-report-content flex-column-gap-sm">
          <!-- Top Row: Coordinates -->
          <div class="quality-report-row">
            <span class="quality-label">${t('quality.sectionCoords')}</span>
            <div class="quality-report-badges">
              <span class="status-pill status-pill-neutral">
                ${t('quality.total')} (${this.lastSummary ? this.lastSummary.totalCount : 0})
              </span>
              <span class="status-pill status-pill-neutral">
                ${t('quality.dim2D')} (${this.lastSummary ? this.lastSummary.count2D : 0})
              </span>
              <span class="status-pill status-pill-neutral">
                ${t('quality.dim3D')} (${this.lastSummary ? this.lastSummary.count3D : 0})
              </span>
            </div>
          </div>

          <!-- Discrete Dashed Separator Line -->
          <div class="quality-report-divider-dashed"></div>

          <!-- Bottom Row: Transformation -->
          <div class="quality-report-row">
            <span class="quality-label">${t('quality.sectionTransform')}</span>
            <div class="quality-report-badges">
              ${
                this.lastSummary && this.lastSummary.inGridCount > 0
                  ? `<span class="status-pill success">${t('quality.valid')} (${this.lastSummary.inGridCount})</span>`
                  : ''
              }
              ${
                this.lastSummary && this.lastSummary.outOfGridCount > 0
                  ? `<span class="status-pill warning">${t('quality.outsideGrid')} (${this.lastSummary.outOfGridCount})</span>`
                  : ''
              }
              ${
                this.lastSummary && this.lastSummary.noDataCount > 0
                  ? `<span class="status-pill warning">${t('quality.noGridData')} (${this.lastSummary.noDataCount})</span>`
                  : ''
              }
              ${
                this.lastSummary && this.lastSummary.invalidAngleCount > 0
                  ? `<span class="status-pill danger">${t('quality.invalidAngle')} (${this.lastSummary.invalidAngleCount})</span>`
                  : ''
              }
              ${
                this.lastSummary && this.lastSummary.invalidInputCount > 0
                  ? `<span class="status-pill danger">${t('quality.invalidInput')} (${this.lastSummary.invalidInputCount})</span>`
                  : ''
              }
              ${
                !this.lastSummary
                  ? `<span class="status-pill neutral">${t('quality.readyToTransform')}</span>`
                  : ''
              }
            </div>
          </div>
        </div>
      </div>

      <!-- Decoupled Map Navigation Dock -->
      <div class="workbench-actions-dock">
        <button 
          class="btn btn-outline center-map-btn btn-secondary-action" 
          id="batch-center-map-btn" 
          ${this.lastMapPoints.length === 0 ? 'disabled' : ''}
        >
          <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <circle cx="12" cy="12" r="10"></circle>
            <line x1="22" y1="12" x2="18" y2="12"></line>
            <line x1="6" y1="12" x2="2" y2="12"></line>
            <line x1="12" y1="6" x2="12" y2="2"></line>
            <line x1="12" y1="22" x2="12" y2="18"></line>
          </svg>
          <span>${t('common.centerOnMap')}</span>
        </button>
        <button 
          class="btn btn-outline zoom-extents-btn btn-secondary-action" 
          id="batch-zoom-extents-btn" 
          ${this.lastMapPoints.length === 0 ? 'disabled' : ''}
        >
          <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <circle cx="11" cy="11" r="8"></circle>
            <line x1="21" y1="21" x2="16.65" y2="16.65"></line>
            <line x1="11" y1="8" x2="11" y2="14"></line>
            <line x1="8" y1="11" x2="14" y2="11"></line>
          </svg>
          <span>${t('common.zoomToExtents')}</span>
        </button>
      </div>
    `;

    this.attachEventListeners();
    this.restoreTextareaHeights();
  }

  private attachEventListeners(): void {
    // Input sync and scrollbar hover cursor detection
    const inputArea = this.container.querySelector('#batch-input-area') as HTMLTextAreaElement | null;
    if (inputArea) {
      inputArea.addEventListener('mousemove', (e: MouseEvent) => {
        const isOverScrollbar = e.offsetX > inputArea.clientWidth || e.offsetY > inputArea.clientHeight;
        inputArea.style.cursor = isOverScrollbar ? 'default' : 'text';
      });
      inputArea.addEventListener('mouseleave', () => {
        inputArea.style.cursor = 'text';
      });
      inputArea.addEventListener('mouseup', () => {
        if (inputArea.style.height) this.savedInputHeight = inputArea.style.height;
      });
      inputArea.addEventListener('input', () => {
        this.inputText = inputArea.value;
      });
    }

    const outputArea = this.container.querySelector('#batch-output-area') as HTMLTextAreaElement | null;
    if (outputArea) {
      outputArea.addEventListener('mouseup', () => {
        if (outputArea.style.height) this.savedOutputHeight = outputArea.style.height;
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

    // Export buttons: Native download triggers silently without premature toasts
    this.container.querySelector('#batch-export-csv')?.addEventListener('click', () => {
      const detectedDelim = this.delimiter === 'auto' ? detectDelimiter(this.inputText.split(/\r?\n/)) : this.delimiter;
      const csvSep = detectedDelim === ' ' ? ',' : detectedDelim;
      exportToCSV(this.lastTransformedItems, csvSep, this.order);
    });

    this.container.querySelector('#batch-export-geojson')?.addEventListener('click', () => {
      exportToGeoJSON(this.lastTransformedItems);
    });

    this.container.querySelector('#batch-export-dxf')?.addEventListener('click', () => {
      exportToDXF(this.lastTransformedItems);
    });
  }
}
