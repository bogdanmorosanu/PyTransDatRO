import { t, onLanguageChange } from '../i18n';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { parseInputLines, type CoordinateOrder } from '../utils/pointParser';
import { formatDegToDMS, degToRad } from '../utils/dmsParser';
import { exportToCSV, exportToGeoJSON, exportToDXF, type ExportPointItem } from '../utils/exporters';
import { showToast } from '../utils/toast';

export interface BatchReportData {
  total: number;
  valid: number;
  outOfGrid: number;
  invalid: number;
  items: ExportPointItem[];
}

export interface BatchConverterCallbacks {
  onPointsPlotted: (
    points: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }>
  ) => void;
  onOpenReport: (report: BatchReportData) => void;
}

export class BatchConverter {
  private container: HTMLElement;
  private callbacks: BatchConverterCallbacks;

  // Settings
  private op: TransformationOp = 'Stereo70ToETRS89';
  private order: CoordinateOrder = 'NE';
  private delimiter: string = 'auto';
  private angularUnit: 'degrees' | 'dms' | 'radians' = 'degrees';

  // State
  private inputText: string = '';
  private outputText: string = '';
  private lastTransformedItems: ExportPointItem[] = [];
  private isProcessing: boolean = false;
  private viewLayout: 'split' | 'stacked' = 'stacked';

  constructor(callbacks: BatchConverterCallbacks) {
    this.callbacks = callbacks;
    this.container = document.createElement('div');
    this.container.className = 'tab-content';

    // Sample default points for immediate exploration
    this.inputText = [
      '101, 500000.000, 500000.000, 100.000',
      '102, 500250.000, 500150.000, 102.500',
      '103, 501000.000, 499500.000, 98.300',
      '104, 5.000, 1.000, 0.000', // Intentional out of grid sample
    ].join('\n');

    this.render();
    onLanguageChange(() => this.render());
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  private async processBatch(): Promise<void> {
    if (this.isProcessing) return;
    if (!this.inputText.trim()) {
      showToast('Please enter coordinate rows to transform', 'warning');
      return;
    }

    this.isProcessing = true;
    this.render();

    try {
      const isSourceGeo = this.op === 'ETRS89ToStereo70';
      const parsedRows = parseInputLines(this.inputText, this.delimiter, this.order, isSourceGeo);

      const validRows = parsedRows.filter((r) => !r.error && r.coords.length >= 2);
      const invalidRows = parsedRows.filter((r) => r.error);

      if (validRows.length === 0) {
        throw new Error('No valid coordinate rows could be parsed');
      }

      // Prepare batch payload: array of [North/Lat, East/Long, (Elev)]
      const payloadPoints = validRows.map((r) => r.coords);

      const res = await apiClient.transformBatch({
        op: this.op,
        points: payloadPoints,
        unit: 'degrees',
      });

      // Assemble results
      const items: ExportPointItem[] = [];
      const outputLines: string[] = [];
      const mapPoints: Array<{ id: string | null; lat: number; lon: number; elev?: number; warning?: string | null }> = [];

      let validCount = 0;
      let outOfGridCount = 0;

      // Match responses with original valid rows
      for (let i = 0; i < validRows.length; i++) {
        const row = validRows[i];
        const apiRes = res.results[i];

        const isGeoOutput = this.op === 'Stereo70ToETRS89';
        const isOutOfGrid = apiRes.warning === 'Out of grid';

        if (isOutOfGrid) outOfGridCount++;
        else validCount++;

        const item: ExportPointItem = {
          pointId: row.pointId,
          coords: apiRes.coos,
          warning: apiRes.warning,
          isGeo: isGeoOutput,
        };
        items.push(item);

        // Format for output textarea
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
            latStr = degToRad(lat).toFixed(12);
            lonStr = degToRad(lon).toFixed(12);
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

          mapPoints.push({ id: row.pointId, lat, lon, elev, warning: apiRes.warning });
        } else {
          // Stereo 70 output
          const n = apiRes.coos[0];
          const e = apiRes.coos[1];
          const z = apiRes.coos[2];

          if (this.order === 'EN') {
            parts.push(e.toFixed(3), n.toFixed(3));
          } else {
            parts.push(n.toFixed(3), e.toFixed(3));
          }

          if (z !== undefined) parts.push(z.toFixed(3));

          // For map plotting Stereo70 outputs, we approximate lat/lon from original inputs or plot coords
          const origLat = row.coords[0];
          const origLon = row.coords[1];
          mapPoints.push({ id: row.pointId, lat: origLat, lon: origLon, elev: z, warning: apiRes.warning });
        }

        if (apiRes.warning) {
          parts.push(`# ${apiRes.warning}`);
        }

        const outDelim = this.delimiter === 'auto' || this.delimiter === ' ' ? ', ' : `${this.delimiter} `;
        outputLines.push(parts.join(outDelim));
      }

      // Add errors
      for (const inv of invalidRows) {
        outputLines.push(`${inv.rawLine} # Error: ${inv.error}`);
      }

      this.outputText = outputLines.join('\n');
      this.lastTransformedItems = items;

      // Plot markers on map
      this.callbacks.onPointsPlotted(mapPoints);

      // Trigger summary report modal
      this.callbacks.onOpenReport({
        total: parsedRows.length,
        valid: validCount,
        outOfGrid: outOfGridCount,
        invalid: invalidRows.length,
        items,
      });

      showToast(`Batch completed: ${validCount} valid coordinates`, 'success');
    } catch (err: any) {
      showToast(err.message || 'Error processing batch', 'danger');
    } finally {
      this.isProcessing = false;
      this.render();
    }
  }

  private render(): void {
    const isS70 = this.op === 'Stereo70ToETRS89';

    this.container.innerHTML = `
      <div class="studio-card">
        <div class="card-title">
          <span>${t('batch.title')}</span>
          <div style="display: flex; gap: 4px; align-items: center;">
            <button class="btn btn-outline" id="batch-layout-toggle" style="padding: 3px 8px; font-size: 11px;">
              ${this.viewLayout === 'stacked' ? '⬍ ' + (t('batch.layoutStacked') || 'Suprapus') : '⬄ ' + (t('batch.layoutSplit') || 'Alăturat')}
            </button>
          </div>
        </div>

        <!-- Direction & Order Grid -->
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 12px;">
          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('point.opLabel')}</label>
            <div class="segmented-control">
              <button class="segmented-btn ${isS70 ? 'active' : ''}" id="batch-op-s70">${t('point.opS70ToETRS')}</button>
              <button class="segmented-btn ${!isS70 ? 'active' : ''}" id="batch-op-etrs">${t('point.opETRSToS70')}</button>
            </div>
          </div>

          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('batch.orderLabel')}</label>
            <div class="segmented-control">
              <button class="segmented-btn ${this.order === 'NE' ? 'active' : ''}" id="batch-order-ne">${t('batch.orderNE')}</button>
              <button class="segmented-btn ${this.order === 'EN' ? 'active' : ''}" id="batch-order-en">${t('batch.orderEN')}</button>
            </div>
          </div>
        </div>

        <!-- Delimiter & Units Grid -->
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 14px;">
          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('batch.delimiterLabel')}</label>
            <select class="form-input" id="batch-delim-select" style="padding: 7px 10px;">
              <option value="auto" ${this.delimiter === 'auto' ? 'selected' : ''}>${t('batch.delimAuto')}</option>
              <option value="," ${this.delimiter === ',' ? 'selected' : ''}>${t('batch.delimComma')}</option>
              <option value=";" ${this.delimiter === ';' ? 'selected' : ''}>${t('batch.delimSemi')}</option>
              <option value=" " ${this.delimiter === ' ' ? 'selected' : ''}>${t('batch.delimSpace')}</option>
              <option value="\t" ${this.delimiter === '\t' ? 'selected' : ''}>${t('batch.delimTab')}</option>
            </select>
          </div>

          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('point.unitLabel')}</label>
            <div class="segmented-control">
              <button class="segmented-btn ${this.angularUnit === 'degrees' ? 'active' : ''}" id="batch-unit-deg">${t('point.unitDeg')}</button>
              <button class="segmented-btn ${this.angularUnit === 'dms' ? 'active' : ''}" id="batch-unit-dms">${t('point.unitDMS')}</button>
            </div>
          </div>
        </div>

        <!-- Textareas (Stacked or Split) -->
        <div style="display: grid; grid-template-columns: ${this.viewLayout === 'split' ? '1fr 1fr' : '1fr'}; gap: 12px;">
          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('batch.inputLabel')}</label>
            <textarea class="form-input mono" id="batch-input-area" rows="${this.viewLayout === 'split' ? '12' : '7'}" style="resize: vertical; font-size: 12px;">${this.inputText}</textarea>
          </div>

          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('batch.outputLabel')}</label>
            <textarea class="form-input mono" id="batch-output-area" rows="${this.viewLayout === 'split' ? '12' : '7'}" readonly style="resize: vertical; font-size: 12px; background-color: var(--bg-surface-elevated);">${this.outputText}</textarea>
          </div>
        </div>

        <!-- Action Row -->
        <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 14px;">
          <button class="btn btn-primary" id="batch-transform-btn" style="min-width: 160px;">
            <span>${this.isProcessing ? '⏳ ' + t('batch.calculating') : '⚡ ' + t('batch.transformBtn')}</span>
          </button>

          <!-- Export Suite -->
          <div style="display: flex; gap: 6px;">
            <button class="btn btn-secondary" id="batch-export-csv" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              📥 CSV
            </button>
            <button class="btn btn-secondary" id="batch-export-geojson" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              🗺️ GeoJSON
            </button>
            <button class="btn btn-secondary" id="batch-export-dxf" ${this.lastTransformedItems.length === 0 ? 'disabled' : ''}>
              📐 DXF
            </button>
          </div>
        </div>
      </div>
    `;

    this.attachEventListeners();
  }

  private attachEventListeners(): void {
    // Layout toggle
    this.container.querySelector('#batch-layout-toggle')?.addEventListener('click', () => {
      this.viewLayout = this.viewLayout === 'stacked' ? 'split' : 'stacked';
      this.render();
    });

    // Op switches
    this.container.querySelector('#batch-op-s70')?.addEventListener('click', () => {
      this.op = 'Stereo70ToETRS89';
      this.render();
    });
    this.container.querySelector('#batch-op-etrs')?.addEventListener('click', () => {
      this.op = 'ETRS89ToStereo70';
      this.render();
    });

    // Order switches
    this.container.querySelector('#batch-order-ne')?.addEventListener('click', () => {
      this.order = 'NE';
      this.render();
    });
    this.container.querySelector('#batch-order-en')?.addEventListener('click', () => {
      this.order = 'EN';
      this.render();
    });

    // Delimiter select
    this.container.querySelector('#batch-delim-select')?.addEventListener('change', (e) => {
      this.delimiter = (e.target as HTMLSelectElement).value;
    });

    // Units
    this.container.querySelector('#batch-unit-deg')?.addEventListener('click', () => {
      this.angularUnit = 'degrees';
      this.render();
    });
    this.container.querySelector('#batch-unit-dms')?.addEventListener('click', () => {
      this.angularUnit = 'dms';
      this.render();
    });

    // Textarea sync
    const inputArea = this.container.querySelector('#batch-input-area') as HTMLTextAreaElement;
    if (inputArea) {
      inputArea.addEventListener('input', () => {
        this.inputText = inputArea.value;
      });
    }

    // Transform
    this.container.querySelector('#batch-transform-btn')?.addEventListener('click', () => {
      if (inputArea) this.inputText = inputArea.value;
      this.processBatch();
    });

    // Exports
    this.container.querySelector('#batch-export-csv')?.addEventListener('click', () => {
      exportToCSV(this.lastTransformedItems, ',', this.order);
      showToast('CSV downloaded', 'success');
    });

    this.container.querySelector('#batch-export-geojson')?.addEventListener('click', () => {
      exportToGeoJSON(this.lastTransformedItems);
      showToast('GeoJSON downloaded', 'success');
    });

    this.container.querySelector('#batch-export-dxf')?.addEventListener('click', () => {
      exportToDXF(this.lastTransformedItems);
      showToast('AutoCAD DXF downloaded', 'success');
    });
  }
}
