import { t, onLanguageChange } from '../i18n';
import { formatTelemetryDecimal } from '../utils/numberFormat';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { type CoordinateOrder } from '../utils/pointParser';
import { parseAngleToDeg } from '../utils/dmsParser';
import { showToast } from '../utils/toast';
import { type MapPointItem } from './MapView';

export type AngularUnitType = 'degrees' | 'dms' | 'radians';

export interface FileUploadCallbacks {
  onPointsPlotted?: (points: MapPointItem[]) => void;
  onCenterOnMap?: (points: MapPointItem[]) => void;
  onZoomToExtents?: (points: MapPointItem[]) => void;
  onClearFiles?: () => void;
}

interface QueuedFileItem {
  id: string;
  file: File;
  status: 'pending' | 'processing' | 'completed' | 'error';
  errorMsg?: string;
  downloadUrl?: string;
  downloadName?: string;
  totalRows?: number;
  validRows?: number;
  outOfGridRows?: number;
  errorRows?: number;
  count2D?: number;
  count3D?: number;
  points?: MapPointItem[];
}

export class FileUpload {
  private container: HTMLElement;
  private callbacks?: FileUploadCallbacks;

  private op: TransformationOp = 'ETRS89ToStereo70';
  private angularUnit: AngularUnitType = 'degrees';
  private order: CoordinateOrder = 'NE';

  private fileQueue: QueuedFileItem[] = [];
  private isProcessing: boolean = false;

  constructor(callbacks?: FileUploadCallbacks) {
    this.callbacks = callbacks;
    this.container = document.createElement('div');
    this.container.className = 'tab-content';
    this.render();

    onLanguageChange(() => this.render());
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public setDirection(op: TransformationOp): void {
    if (this.op === op) return;
    this.op = op;
    this.render();
  }

  public setAngularUnit(unit: AngularUnitType): void {
    if (this.angularUnit === unit) return;
    this.angularUnit = unit;
    this.render();
  }

  public setCoordinateOrder(order: CoordinateOrder): void {
    this.order = order;
  }

  public getCoordinateOrder(): CoordinateOrder {
    return this.order;
  }

  public clearQueue(): void {
    this.fileQueue = [];
    this.callbacks?.onClearFiles?.();
    this.render();
  }

  private addFilesToQueue(files: FileList | File[]): void {
    const fileArray = Array.from(files);
    let addedCount = 0;

    for (const file of fileArray) {
      if (file.size > 20 * 1024 * 1024) {
        showToast(t('file.errSizeLimit', { name: file.name }), 'danger');
        continue;
      }

      this.fileQueue.push({
        id: `${file.name}-${file.size}-${Date.now()}-${Math.random()}`,
        file,
        status: 'pending',
      });
      addedCount++;
    }

    this.render();
  }

  private removeFileFromQueue(id: string): void {
    const idx = this.fileQueue.findIndex((f) => f.id === id);
    if (idx === -1) return;
    this.fileQueue.splice(idx, 1);
    if (this.fileQueue.length === 0) {
      this.callbacks?.onClearFiles?.();
    } else {
      const remainingPoints = this.fileQueue.flatMap((f) => f.points || []);
      if (remainingPoints.length > 0) {
        this.callbacks?.onPointsPlotted?.(remainingPoints);
      } else {
        this.callbacks?.onClearFiles?.();
      }
    }
    this.render();
  }


  private parsePairedPoints(
    inText: string,
    outText: string,
    fileName: string
  ): MapPointItem[] {
    const inLines = inText.split(/\r?\n/);
    const outLines = outText.split(/\r?\n/);
    const points: MapPointItem[] = [];
    const maxRows = Math.min(inLines.length, outLines.length);

    const parseLineTokens = (raw: string): { pointId: string | null; tokens: string[] } => {
      const clean = raw.split('#')[0].trim();
      if (!clean) return { pointId: null, tokens: [] };
      const delim = clean.includes('\t') ? '\t' : clean.includes(';') ? ';' : clean.includes(',') ? ',' : /\s+/;
      const allTokens = clean.split(delim).map((t) => t.trim()).filter(Boolean);
      if (allTokens.length < 2) return { pointId: null, tokens: [] };

      if (allTokens.length >= 4) {
        return { pointId: allTokens[0], tokens: allTokens.slice(1, 4) };
      } else if (allTokens.length === 3) {
        const t0 = allTokens[0];
        const t1 = allTokens[1];
        const t2 = allTokens[2];
        const hasAngle = (s: string) => /[°'"dms]/i.test(s);
        if (hasAngle(t1) || (!hasAngle(t0) && isNaN(Number(t0)))) {
          return { pointId: t0, tokens: [t1, t2] };
        } else if (parseFloat(t0) < 50000 && parseFloat(t1) > 50000) {
          return { pointId: t0, tokens: [t1, t2] };
        }
      }
      return { pointId: null, tokens: allTokens };
    };

    let rowIdx = 0;
    for (let i = 0; i < maxRows; i++) {
      const inRaw = inLines[i]?.trim();
      const outRaw = outLines[i]?.trim();
      if (!inRaw || inRaw.startsWith('#')) continue;
      rowIdx++;

      const inParsed = parseLineTokens(inRaw);
      const outParsed = parseLineTokens(outRaw || '');
      if (inParsed.tokens.length < 2) continue;

      const pointId = inParsed.pointId || outParsed.pointId || `Point ${rowIdx}`;
      const isOutOfGrid = (outRaw || '').toLowerCase().includes('grid');

      if (this.op === 'Stereo70ToETRS89') {
        const inC1 = parseFloat(inParsed.tokens[0]);
        const inC2 = parseFloat(inParsed.tokens[1]);
        const inElev = inParsed.tokens[2] ? parseFloat(inParsed.tokens[2]) : undefined;

        const stereoN = this.order === 'EN' ? inC2 : inC1;
        const stereoE = this.order === 'EN' ? inC1 : inC2;

        if (outParsed.tokens.length >= 2) {
          let outLat: number;
          let outLon: number;
          try {
            outLat = parseAngleToDeg(outParsed.tokens[0]);
            outLon = parseAngleToDeg(outParsed.tokens[1]);
          } catch {
            continue;
          }
          const outH = outParsed.tokens[2] ? parseFloat(outParsed.tokens[2]) : inElev;
          const alt = inElev;
          const h = outH;
          const zeta = (h !== undefined && alt !== undefined) ? (h - alt) : undefined;

          points.push({
            id: pointId,
            fileName,
            lat: outLat,
            lon: outLon,
            stereoN,
            stereoE,
            elev: outH !== undefined ? outH : inElev,
            altitude: alt,
            ellipsoidalHeight: h,
            zeta,
            warning: isOutOfGrid ? 'Out of grid' : null,
            angularUnit: this.angularUnit,
          });
        }
      } else {
        // ETRS89ToStereo70
        let inLat: number;
        let inLon: number;
        try {
          inLat = parseAngleToDeg(inParsed.tokens[this.order === 'EN' ? 1 : 0]);
          inLon = parseAngleToDeg(inParsed.tokens[this.order === 'EN' ? 0 : 1]);
        } catch {
          continue;
        }
        const inH = inParsed.tokens[2] ? parseFloat(inParsed.tokens[2]) : undefined;

        let outN: number | undefined;
        let outE: number | undefined;
        let outElev: number | undefined;

        if (outParsed.tokens.length >= 2) {
          outN = parseFloat(outParsed.tokens[0]);
          outE = parseFloat(outParsed.tokens[1]);
          outElev = outParsed.tokens[2] ? parseFloat(outParsed.tokens[2]) : inH;
        }

        const alt = outElev;
        const h = inH;
        const zeta = (h !== undefined && alt !== undefined) ? (h - alt) : undefined;

        points.push({
          id: pointId,
          fileName,
          lat: inLat,
          lon: inLon,
          stereoN: outN,
          stereoE: outE,
          elev: outElev !== undefined ? outElev : inH,
          altitude: alt,
          ellipsoidalHeight: h,
          zeta,
          warning: isOutOfGrid ? 'Out of grid' : null,
          angularUnit: this.angularUnit,
        });
      }
    }

    return points;
  }

  private async processQueue(): Promise<void> {
    if (this.isProcessing) return;
    const pendingItems = this.fileQueue.filter((item) => item.status === 'pending');
    if (pendingItems.length === 0) {
      showToast(t('file.warnNoPending'), 'warning');
      return;
    }

    this.isProcessing = true;
    this.render();

    for (const item of pendingItems) {
      item.status = 'processing';
      this.render();

      try {
        const response = await apiClient.transformFile(item.file, this.op, this.angularUnit, 'auto');
        const blob = await response.blob();
        const text = await blob.text();
        item.downloadUrl = URL.createObjectURL(blob);
        item.downloadName = `transformed_${item.file.name}`;
        item.status = 'completed';

        // Parse per-file quality metrics
        const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l.length > 0 && !l.startsWith('//'));
        let valCount = 0;
        let outGridCount = 0;
        let errCount = 0;
        for (const line of lines) {
          if (line.startsWith('#')) continue;
          if (line.includes('# Out of grid') || line.includes('Out of grid') || line.includes('No data on grid')) {
            outGridCount++;
          } else if (
            line.includes('# Wrong dimension') ||
            line.includes('# Invalid coordinate data') ||
            line.includes('# Error') ||
            line.includes('Error:') ||
            /#\s*(Wrong dimension|Invalid coordinate|Limit exceeded)/i.test(line)
          ) {
            errCount++;
          } else if (line.includes('#')) {
            const afterHash = line.split('#')[1].trim();
            if (afterHash.toLowerCase().includes('grid')) {
              outGridCount++;
            } else {
              errCount++;
            }
          } else {
            valCount++;
          }
        }
        item.totalRows = valCount + outGridCount + errCount;
        item.validRows = valCount;
        item.outOfGridRows = outGridCount;
        item.errorRows = errCount;

        try {
          const inText = await item.file.text();
          item.points = this.parsePairedPoints(inText, text, item.file.name);
          let c2D = 0;
          let c3D = 0;
          const inLines = inText.split(/\r?\n/).map((l) => l.trim()).filter((l) => l.length > 0 && !l.startsWith('#') && !l.startsWith('//'));
          for (const line of inLines) {
            const tokens = line.split(/[;, \t]+/).filter(Boolean);
            if (tokens.length >= 4) {
              c3D++;
            } else if (tokens.length === 3) {
              const t2Num = parseFloat(tokens[2]);
              if (!isNaN(t2Num)) {
                c3D++;
              } else {
                c2D++;
              }
            } else if (tokens.length === 2) {
              c2D++;
            }
          }
          item.count2D = c2D;
          item.count3D = c3D;
        } catch {
          item.points = [];
          item.count2D = 0;
          item.count3D = 0;
        }
      } catch (err: any) {
        item.status = 'error';
        item.errorMsg = err.message || 'Transformation error';
      }
      this.render();
    }

    this.isProcessing = false;
    this.render();

    const cumulativePoints = this.fileQueue.flatMap((f) => f.points || []);
    if (cumulativePoints.length > 0) {
      this.callbacks?.onPointsPlotted?.(cumulativePoints);
    }

    const successCount = this.fileQueue.filter((f) => f.status === 'completed').length;
    showToast(t('file.successBatch', { count: successCount }), 'success');
  }

  private downloadFile(item: QueuedFileItem): void {
    if (!item.downloadUrl || !item.downloadName) return;
    const a = document.createElement('a');
    a.href = item.downloadUrl;
    a.download = item.downloadName;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  private downloadAll(): void {
    const completed = this.fileQueue.filter((item) => item.status === 'completed');
    completed.forEach((item, index) => {
      setTimeout(() => this.downloadFile(item), index * 300);
    });
  }

  private render(): void {
    const hasFiles = this.fileQueue.length > 0;
    const completedCount = this.fileQueue.filter((f) => f.status === 'completed').length;
    const cumulativePoints = this.fileQueue.flatMap((f) => f.points || []);

    const totalValid = this.fileQueue.reduce((acc, f) => acc + (f.validRows || 0), 0);
    const totalOutGrid = this.fileQueue.reduce((acc, f) => acc + (f.outOfGridRows || 0), 0);
    const totalErrors = this.fileQueue.reduce((acc, f) => acc + (f.errorRows || 0), 0);
    const total2D = this.fileQueue.reduce((acc, f) => acc + (f.count2D || 0), 0);
    const total3D = this.fileQueue.reduce((acc, f) => acc + (f.count3D || 0), 0);
    const totalPoints = totalValid + totalOutGrid + totalErrors;

    this.container.innerHTML = `
      <div class="studio-card">
        ${
          hasFiles
            ? `
          <div style="display: flex; justify-content: flex-end; margin-bottom: 10px;">
            <button class="btn btn-outline" id="file-clear-queue-btn" style="padding: 3px 8px; font-size: 11px;">
              ✕ ${t('file.clearQueue')}
            </button>
          </div>
        `
            : ''
        }

        <!-- 1. Drag and Drop Dropzone -->
        <div id="file-dropzone" class="file-dropzone">
          <div style="font-weight: 600; font-size: 13.5px; margin-bottom: 4px; color: var(--text-main);">
            ${t('file.dropzoneTitle')}
          </div>
          <div style="color: var(--text-muted); font-size: 11.5px;">
            ${t('file.dropzoneHint')}
          </div>
          <input type="file" id="file-input-hidden" aria-label="${t('file.dropzoneTitle')}" multiple accept=".csv,.txt,.xyz,.dat" style="display: none;" />
        </div>

        <!-- 2. Queued Files List with File Details and Per-File Quality Metrics -->
        ${
          hasFiles
            ? `
          <div class="file-queue">
            <div style="font-size: 11.5px; font-weight: 600; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.4px;">
              ${t('file.queueTitle')} (${this.fileQueue.length})
            </div>
            ${this.fileQueue
              .map((item) => {
                const sizeKb = formatTelemetryDecimal(item.file.size / 1024, 1);
                let badgeClass = 'pending';
                let statusLabel = t('file.statusPending');
                if (item.status === 'processing') {
                  badgeClass = 'processing';
                  statusLabel = t('file.statusProcessing');
                } else if (item.status === 'completed') {
                  badgeClass = 'completed';
                  statusLabel = t('file.statusCompleted');
                } else if (item.status === 'error') {
                  badgeClass = 'error';
                  statusLabel = t('file.statusError');
                }

                let metaDetails = `${sizeKb} KB`;
                if (item.status === 'completed' && item.validRows !== undefined) {
                  metaDetails += ` • ${t('quality.validCoords')}: ${item.validRows}`;
                  if (item.outOfGridRows && item.outOfGridRows > 0) {
                    metaDetails += ` • ⚠ ${t('quality.outsideGrid')}: ${item.outOfGridRows}`;
                  }
                  if (item.errorRows && item.errorRows > 0) {
                    metaDetails += ` • ✕ ${item.errorRows} errors`;
                  }
                } else if (item.status === 'error' && item.errorMsg) {
                  metaDetails += ` • ${item.errorMsg}`;
                }

                return `
                  <div class="file-queue-item" data-id="${item.id}">
                    <div class="file-queue-info">
                      <div class="file-queue-name" title="${item.file.name}">${item.file.name}</div>
                      <div class="file-queue-meta">${metaDetails}</div>
                    </div>
                    <div class="file-queue-actions">
                      <span class="file-status-badge ${badgeClass}">${statusLabel}</span>
                      ${
                        item.status === 'completed'
                          ? `
                        <button class="btn btn-outline file-download-btn" data-id="${item.id}" style="padding: 3px 8px; font-size: 11px;">
                          ${t('file.downloadBtn')}
                        </button>
                      `
                          : ''
                      }
                      <button class="btn-icon file-queue-remove" data-id="${item.id}" title="${t('file.removeItem')}" aria-label="${t('file.removeItem')}">
                        ✕
                      </button>
                    </div>
                  </div>
                `;
              })
              .join('')}
          </div>
        `
            : ''
        }

        <!-- 3. Primary Action Button -->
        <button 
          class="btn btn-primary" 
          id="file-start-btn" 
          style="width: 100%; margin-top: 14px;" 
          ${!hasFiles || this.isProcessing ? 'disabled' : ''}
        >
          <span>${this.isProcessing ? t('file.streamingProgress') : t('file.transformBtn')}</span>
        </button>

        <!-- 4. Cumulative Quality Report & Download Actions -->
        ${
          completedCount > 0
            ? `
          <div class="quality-report-card" style="margin-top: 14px;">
            <!-- Top Section: Coordinates -->
            <div class="quality-report-section coords-section">
              <span class="quality-section-header">${t('quality.sectionCoords')}</span>
              <div class="quality-report-badges">
                <span class="status-pill" style="background-color: var(--bg-surface-elevated); color: var(--text-muted); border: 1px solid var(--border-medium);">
                  ${t('quality.total')}: ${totalPoints}
                </span>
                <span class="status-pill" style="background-color: var(--bg-surface-elevated); color: var(--text-muted); border: 1px solid var(--border-medium);">
                  ${t('quality.dim2D')}: ${total2D}
                </span>
                <span class="status-pill" style="background-color: var(--bg-surface-elevated); color: var(--text-muted); border: 1px solid var(--border-medium);">
                  ${t('quality.dim3D')}: ${total3D}
                </span>
              </div>
            </div>

            <!-- Discrete Dashed Separator Line -->
            <div class="quality-report-divider-dashed"></div>

            <!-- Bottom Section: Transformation -->
            <div class="quality-report-section transform-section">
              <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
                <span class="quality-section-header">${t('quality.sectionTransform')}</span>
                ${
                  completedCount > 1
                    ? `
                  <button class="btn btn-secondary" id="file-download-all-btn" style="padding: 3px 8px; font-size: 11px;">
                    <span>${t('file.downloadAll')}</span>
                  </button>
                `
                    : ''
                }
              </div>
              <div class="quality-report-badges">
                ${
                  totalValid > 0
                    ? `<span class="status-pill success">✓ ${t('quality.valid')} (${totalValid})</span>`
                    : ''
                }
                ${
                  totalOutGrid > 0
                    ? `<span class="status-pill warning">⚠ ${t('quality.outsideGrid')} (${totalOutGrid})</span>`
                    : ''
                }
                ${
                  totalErrors > 0
                    ? `<span class="status-pill danger">✕ ${t('quality.invalidInput')} (${totalErrors})</span>`
                    : ''
                }
              </div>
            </div>
          </div>

          <!-- 5. Bottom-Right Navigation Actions -->
          <div class="action-footer-right" style="display: flex; gap: 8px; justify-content: flex-end; width: 100%; margin-top: 10px;">
            <button 
              class="btn btn-outline center-map-btn" 
              id="file-center-map-btn" 
              ${cumulativePoints.length === 0 ? 'disabled' : ''} 
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
              id="file-zoom-extents-btn" 
              ${cumulativePoints.length === 0 ? 'disabled' : ''} 
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
    const dropzone = this.container.querySelector('#file-dropzone');
    const hiddenInput = this.container.querySelector('#file-input-hidden') as HTMLInputElement;

    if (dropzone && hiddenInput) {
      dropzone.addEventListener('click', () => hiddenInput.click());

      dropzone.addEventListener('dragover', (e: any) => {
        e.preventDefault();
        dropzone.classList.add('drag-over');
      });

      dropzone.addEventListener('dragleave', (e: any) => {
        e.preventDefault();
        dropzone.classList.remove('drag-over');
      });

      dropzone.addEventListener('drop', (e: any) => {
        e.preventDefault();
        dropzone.classList.remove('drag-over');
        if (e.dataTransfer && e.dataTransfer.files.length > 0) {
          this.addFilesToQueue(e.dataTransfer.files);
        }
      });

      hiddenInput.addEventListener('change', () => {
        if (hiddenInput.files && hiddenInput.files.length > 0) {
          this.addFilesToQueue(hiddenInput.files);
        }
      });
    }

    this.container.querySelector('#file-clear-queue-btn')?.addEventListener('click', () => {
      this.fileQueue = [];
      this.callbacks?.onClearFiles?.();
      this.render();
    });

    this.container.querySelector('#file-start-btn')?.addEventListener('click', () => {
      this.processQueue();
    });

    this.container.querySelector('#file-download-all-btn')?.addEventListener('click', () => {
      this.downloadAll();
    });

    // Center on map from completed files
    this.container.querySelector('#file-center-map-btn')?.addEventListener('click', () => {
      const cumulativePoints = this.fileQueue.flatMap((f) => f.points || []);
      if (cumulativePoints.length > 0) {
        this.callbacks?.onCenterOnMap?.(cumulativePoints);
      }
    });

    // Zoom to extents from completed files
    this.container.querySelector('#file-zoom-extents-btn')?.addEventListener('click', () => {
      const cumulativePoints = this.fileQueue.flatMap((f) => f.points || []);
      if (cumulativePoints.length > 0) {
        this.callbacks?.onZoomToExtents?.(cumulativePoints);
      }
    });

    // Individual download buttons
    this.container.querySelectorAll('.file-download-btn').forEach((btn) => {
      btn.addEventListener('click', () => {
        const id = btn.getAttribute('data-id');
        const item = this.fileQueue.find((f) => f.id === id);
        if (item) this.downloadFile(item);
      });
    });

    // Individual remove queue buttons
    this.container.querySelectorAll('.file-queue-remove').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const id = btn.getAttribute('data-id');
        if (id) this.removeFileFromQueue(id);
      });
    });
  }
}
