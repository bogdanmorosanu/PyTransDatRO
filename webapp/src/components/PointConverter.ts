import { t, onLanguageChange } from '../i18n';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { parseAngleToDeg, formatDegToDMS, degToRad } from '../utils/dmsParser';
import { showToast } from '../utils/toast';

export type AngularUnitType = 'degrees' | 'dms' | 'radians';

export interface PointConverterCallbacks {
  onPointTransformed: (coords: { lat: number; lon: number; outOfGrid: boolean }) => void;
}

export class PointConverter {
  private container: HTMLElement;
  private callbacks: PointConverterCallbacks;

  private op: TransformationOp = 'Stereo70ToETRS89';
  private is3D: boolean = true;
  private angularUnit: AngularUnitType = 'degrees';

  // State
  private inputCoo1: string = '500000.000';
  private inputCoo2: string = '500000.000';
  private inputCoo3: string = '100.000';

  private outputCoords: number[] | null = null;
  private outputWarning: string | null = null;
  private quasigeoidZeta: number | null = null;
  private isCalculating: boolean = false;

  constructor(callbacks: PointConverterCallbacks) {
    this.callbacks = callbacks;
    this.container = document.createElement('div');
    this.container.className = 'tab-content';

    this.readQueryParams();
    this.render();

    onLanguageChange(() => this.render());

    // Auto calculate initial coordinates if present in URL or default
    setTimeout(() => this.calculate(), 100);
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  public setCoordinatesFromMap(lat: number, lon: number): void {
    if (this.op === 'Stereo70ToETRS89') {
      // User clicked on map, but current mode is Stereo70 -> ETRS89.
      // Auto-switch to ETRS89 -> Stereo70 to make sense of the click!
      this.op = 'ETRS89ToStereo70';
    }

    if (this.angularUnit === 'dms') {
      this.inputCoo1 = formatDegToDMS(lat);
      this.inputCoo2 = formatDegToDMS(lon);
    } else if (this.angularUnit === 'radians') {
      this.inputCoo1 = degToRad(lat).toFixed(12);
      this.inputCoo2 = degToRad(lon).toFixed(12);
    } else {
      this.inputCoo1 = lat.toFixed(7);
      this.inputCoo2 = lon.toFixed(7);
    }

    this.render();
    this.calculate();
  }

  private readQueryParams(): void {
    const params = new URLSearchParams(window.location.search);
    if (params.has('op')) {
      const opParam = params.get('op');
      if (opParam === 'Stereo70ToETRS89' || opParam === 'ETRS89ToStereo70') {
        this.op = opParam;
      }
    }
    if (params.has('c1')) this.inputCoo1 = params.get('c1') || this.inputCoo1;
    if (params.has('c2')) this.inputCoo2 = params.get('c2') || this.inputCoo2;
    if (params.has('c3')) {
      this.inputCoo3 = params.get('c3') || this.inputCoo3;
      this.is3D = true;
    }
  }

  private updateQueryParams(): void {
    const params = new URLSearchParams();
    params.set('op', this.op);
    params.set('c1', this.inputCoo1);
    params.set('c2', this.inputCoo2);
    if (this.is3D) params.set('c3', this.inputCoo3);
    const newUrl = `${window.location.pathname}?${params.toString()}`;
    window.history.replaceState({}, '', newUrl);
  }

  private async calculate(): Promise<void> {
    if (this.isCalculating) return;
    this.isCalculating = true;

    try {
      let raw1: number;
      let raw2: number;
      let raw3: number | undefined;

      if (this.op === 'Stereo70ToETRS89') {
        raw1 = parseFloat(this.inputCoo1);
        raw2 = parseFloat(this.inputCoo2);
      } else {
        // ETRS89 is source: parse angle into decimal degrees
        raw1 = parseAngleToDeg(this.inputCoo1);
        raw2 = parseAngleToDeg(this.inputCoo2);
      }

      if (isNaN(raw1) || isNaN(raw2)) {
        throw new Error('Invalid coordinate inputs');
      }

      if (this.is3D) {
        raw3 = parseFloat(this.inputCoo3);
        if (isNaN(raw3)) raw3 = 0.0;
      }

      const coos = raw3 !== undefined ? [raw1, raw2, raw3] : [raw1, raw2];

      const res = await apiClient.transformPoint({
        op: this.op,
        coos,
        unit: 'degrees',
      });

      this.outputCoords = res.coos;
      this.outputWarning = res.warning;

      // Quasigeoid anomaly calculation in 3D:
      // When Stereo70 -> ETRS89: res.coos[2] is h, input is H => zeta = h - H
      // When ETRS89 -> Stereo70: res.coos[2] is H, input is h => zeta = h - H
      if (this.is3D && this.outputCoords.length >= 3 && raw3 !== undefined) {
        if (this.op === 'Stereo70ToETRS89') {
          this.quasigeoidZeta = this.outputCoords[2] - raw3;
        } else {
          this.quasigeoidZeta = raw3 - this.outputCoords[2];
        }
      } else {
        this.quasigeoidZeta = null;
      }

      this.updateQueryParams();

      // Emit coordinates to map
      let mapLat: number;
      let mapLon: number;
      if (this.op === 'Stereo70ToETRS89') {
        mapLat = this.outputCoords[0];
        mapLon = this.outputCoords[1];
      } else {
        mapLat = raw1;
        mapLon = raw2;
      }

      const isOutOfGrid = res.warning === 'Out of grid';
      this.callbacks.onPointTransformed({ lat: mapLat, lon: mapLon, outOfGrid: isOutOfGrid });

      if (isOutOfGrid) {
        showToast(t('point.outOfGridWarn'), 'warning');
      }
    } catch (err: any) {
      this.outputCoords = null;
      this.outputWarning = err.message || 'Error';
      showToast(err.message || 'Error calculating coordinates', 'danger');
    } finally {
      this.isCalculating = false;
      this.render();
    }
  }

  private render(): void {
    const isS70 = this.op === 'Stereo70ToETRS89';
    const c1Label = isS70 ? t('point.northingLabel') : t('point.latitudeLabel');
    const c2Label = isS70 ? t('point.eastingLabel') : t('point.longitudeLabel');
    const c3Label = isS70 ? t('point.elevationNormal') : t('point.elevationEllipsoidal');

    this.container.innerHTML = `
      <div class="studio-card">
        <div class="card-title">
          <span>${t('point.title')}</span>
          <button class="btn btn-outline" id="point-swap-op-btn" style="padding: 4px 10px; font-size: 11px;">
            ⇄ ${t('point.swapDirection')}
          </button>
        </div>

        <!-- Direction & Mode Toggles -->
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 14px;">
          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('point.opLabel')}</label>
            <div class="segmented-control">
              <button class="segmented-btn ${isS70 ? 'active' : ''}" id="op-s70-btn">
                ${t('point.opS70ToETRS')}
              </button>
              <button class="segmented-btn ${!isS70 ? 'active' : ''}" id="op-etrs-btn">
                ${t('point.opETRSToS70')}
              </button>
            </div>
          </div>

          <div class="form-group" style="margin-bottom: 0;">
            <label class="form-label">${t('point.dimensionLabel')}</label>
            <div class="segmented-control">
              <button class="segmented-btn ${!this.is3D ? 'active' : ''}" id="dim-2d-btn">${t('point.dim2D')}</button>
              <button class="segmented-btn ${this.is3D ? 'active' : ''}" id="dim-3d-btn">${t('point.dim3D')}</button>
            </div>
          </div>
        </div>

        <!-- Angular Units (for ETRS89) -->
        <div class="form-group" style="margin-bottom: 14px;">
          <label class="form-label">${t('point.unitLabel')}</label>
          <div class="segmented-control">
            <button class="segmented-btn ${this.angularUnit === 'degrees' ? 'active' : ''}" id="unit-deg-btn">${t('point.unitDeg')}</button>
            <button class="segmented-btn ${this.angularUnit === 'dms' ? 'active' : ''}" id="unit-dms-btn">${t('point.unitDMS')}</button>
            <button class="segmented-btn ${this.angularUnit === 'radians' ? 'active' : ''}" id="unit-rad-btn">${t('point.unitRad')}</button>
          </div>
        </div>

        <!-- Inputs Grid -->
        <div style="display: grid; grid-template-columns: ${this.is3D ? '1fr 1fr 1fr' : '1fr 1fr'}; gap: 10px;">
          <div class="form-group">
            <label class="form-label">${c1Label}</label>
            <input type="text" class="form-input mono" id="pt-input-1" value="${this.inputCoo1}" />
          </div>
          <div class="form-group">
            <label class="form-label">${c2Label}</label>
            <input type="text" class="form-input mono" id="pt-input-2" value="${this.inputCoo2}" />
          </div>
          ${
            this.is3D
              ? `
            <div class="form-group">
              <label class="form-label">${c3Label}</label>
              <input type="text" class="form-input mono" id="pt-input-3" value="${this.inputCoo3}" />
            </div>
          `
              : ''
          }
        </div>

        <button class="btn btn-primary" id="pt-calculate-btn" style="width: 100%; margin-top: 4px;">
          <span>${this.isCalculating ? '⏳...' : '⚡ ' + t('point.transformBtn')}</span>
        </button>
      </div>

      <!-- Output Result Card -->
      ${this.renderResultCard()}
    `;

    this.attachEventListeners();
  }

  private renderResultCard(): string {
    if (!this.outputCoords) return '';

    const isTargetGeo = this.op === 'Stereo70ToETRS89';
    let val1Str: string;
    let val2Str: string;

    if (isTargetGeo) {
      if (this.angularUnit === 'dms') {
        val1Str = formatDegToDMS(this.outputCoords[0]);
        val2Str = formatDegToDMS(this.outputCoords[1]);
      } else if (this.angularUnit === 'radians') {
        val1Str = degToRad(this.outputCoords[0]).toFixed(12) + ' rad';
        val2Str = degToRad(this.outputCoords[1]).toFixed(12) + ' rad';
      } else {
        val1Str = this.outputCoords[0].toFixed(8) + '°';
        val2Str = this.outputCoords[1].toFixed(8) + '°';
      }
    } else {
      val1Str = this.outputCoords[0].toFixed(3) + ' m';
      val2Str = this.outputCoords[1].toFixed(3) + ' m';
    }

    const val3Str = this.outputCoords[2] !== undefined ? this.outputCoords[2].toFixed(3) + ' m' : null;

    const label1 = isTargetGeo ? t('point.latitudeLabel') : t('point.northingLabel');
    const label2 = isTargetGeo ? t('point.longitudeLabel') : t('point.eastingLabel');
    const label3 = isTargetGeo ? t('point.elevationEllipsoidal') : t('point.elevationNormal');

    return `
      <div class="result-card">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px;">
          <span style="font-weight: 600; font-size: 13px; color: var(--text-main);">${t('point.outputCoords')}</span>
          <button class="btn btn-outline" id="pt-copy-btn" style="padding: 3px 8px; font-size: 11px;">
            📋 ${t('point.copyTuple')}
          </button>
        </div>

        <div class="result-row">
          <span style="color: var(--text-muted); font-size: 12px;">${label1}:</span>
          <span class="result-val">${val1Str}</span>
        </div>
        <div class="result-row">
          <span style="color: var(--text-muted); font-size: 12px;">${label2}:</span>
          <span class="result-val">${val2Str}</span>
        </div>
        ${
          val3Str
            ? `
          <div class="result-row">
            <span style="color: var(--text-muted); font-size: 12px;">${label3}:</span>
            <span class="result-val">${val3Str}</span>
          </div>
        `
            : ''
        }

        ${
          this.quasigeoidZeta !== null
            ? `
          <div class="result-row" style="margin-top: 10px; padding-top: 8px; border-top: 1px dashed var(--border-medium);">
            <span style="color: var(--color-primary-light); font-size: 12px; font-weight: 500;">
              ${t('point.quasigeoidLabel')}:
            </span>
            <span class="result-val" style="color: var(--color-success);">
              ${(this.quasigeoidZeta >= 0 ? '+' : '') + this.quasigeoidZeta.toFixed(3)} m
            </span>
          </div>
        `
            : ''
        }

        ${
          this.outputWarning
            ? `
          <div style="margin-top: 8px; padding: 6px 10px; background-color: var(--color-warning-bg); border: 1px solid var(--color-warning); border-radius: var(--radius-sm); font-size: 12px; color: var(--color-warning);">
            ⚠ ${this.outputWarning}
          </div>
        `
            : ''
        }
      </div>
    `;
  }

  private attachEventListeners(): void {
    // Op toggles
    this.container.querySelector('#op-s70-btn')?.addEventListener('click', () => {
      this.op = 'Stereo70ToETRS89';
      this.render();
      this.calculate();
    });
    this.container.querySelector('#op-etrs-btn')?.addEventListener('click', () => {
      this.op = 'ETRS89ToStereo70';
      this.render();
      this.calculate();
    });
    this.container.querySelector('#point-swap-op-btn')?.addEventListener('click', () => {
      this.op = this.op === 'Stereo70ToETRS89' ? 'ETRS89ToStereo70' : 'Stereo70ToETRS89';
      // Invert input values with previous outputs if available
      if (this.outputCoords) {
        this.inputCoo1 = this.outputCoords[0].toFixed(3);
        this.inputCoo2 = this.outputCoords[1].toFixed(3);
        if (this.outputCoords[2] !== undefined) {
          this.inputCoo3 = this.outputCoords[2].toFixed(3);
        }
      }
      this.render();
      this.calculate();
    });

    // 2D / 3D
    this.container.querySelector('#dim-2d-btn')?.addEventListener('click', () => {
      this.is3D = false;
      this.render();
      this.calculate();
    });
    this.container.querySelector('#dim-3d-btn')?.addEventListener('click', () => {
      this.is3D = true;
      this.render();
      this.calculate();
    });

    // Angular Units
    this.container.querySelector('#unit-deg-btn')?.addEventListener('click', () => {
      this.angularUnit = 'degrees';
      this.render();
    });
    this.container.querySelector('#unit-dms-btn')?.addEventListener('click', () => {
      this.angularUnit = 'dms';
      this.render();
    });
    this.container.querySelector('#unit-rad-btn')?.addEventListener('click', () => {
      this.angularUnit = 'radians';
      this.render();
    });

    // Inputs
    const i1 = this.container.querySelector('#pt-input-1') as HTMLInputElement;
    const i2 = this.container.querySelector('#pt-input-2') as HTMLInputElement;
    const i3 = this.container.querySelector('#pt-input-3') as HTMLInputElement;

    const onInputChange = () => {
      if (i1) this.inputCoo1 = i1.value.trim();
      if (i2) this.inputCoo2 = i2.value.trim();
      if (i3) this.inputCoo3 = i3.value.trim();
    };

    [i1, i2, i3].forEach((input) => {
      if (input) {
        input.addEventListener('input', onInputChange);
        input.addEventListener('keydown', (e) => {
          if (e.key === 'Enter') {
            onInputChange();
            this.calculate();
          }
        });
      }
    });

    this.container.querySelector('#pt-calculate-btn')?.addEventListener('click', () => {
      onInputChange();
      this.calculate();
    });

    // Copy tuple
    this.container.querySelector('#pt-copy-btn')?.addEventListener('click', () => {
      if (!this.outputCoords) return;
      const copyStr = this.outputCoords.join(', ');
      navigator.clipboard.writeText(copyStr);
      showToast(t('point.copiedToast'), 'success');
    });
  }
}
