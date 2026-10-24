import { t, onLanguageChange } from '../i18n';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { parseAngleToDeg, formatDegToDMS, degToRad, radToDeg } from '../utils/dmsParser';
import { isValidDecimal, type CoordinateOrder } from '../utils/pointParser';
import { isWithinStereo70Grid, isWithinETRS89Grid } from '../utils/gridBounds';
import { showToast } from '../utils/toast';
import { isInsideRomaniaGrid } from './MapView';
import { renderCopyButton } from '../ui';

function escapeHtmlAttr(str: string): string {
  return str.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/'/g, '&#39;');
}

export type AngularUnitType = 'degrees' | 'dms' | 'radians';

export interface PointTransformedEvent {
  lat: number;
  lon: number;
  outOfGrid: boolean;
  pointId?: string;
  stereoN?: number;
  stereoE?: number;
  elev?: number;
  altitude?: number;
  ellipsoidalHeight?: number;
  zeta?: number | null;
  angularUnit?: AngularUnitType;
}

export interface PointConverterCallbacks {
  onPointTransformed: (coords: PointTransformedEvent) => void;
  onCenterOnMap?: (coords: { lat: number; lon: number; outOfGrid: boolean }) => void;
  onZoomToExtents?: (coords: { lat: number; lon: number; outOfGrid: boolean }) => void;
  onDirectionChanged?: (op: TransformationOp) => void;
}

export class PointConverter {
  private container: HTMLElement;
  private callbacks: PointConverterCallbacks;

  private op: TransformationOp = 'ETRS89ToStereo70';
  private angularUnit: AngularUnitType = 'degrees';
  private order: CoordinateOrder = 'NE';

  // Coordinate state - Default Stereo 70 projection center origin (46° N, 25° E, h=100.0 m)
  private inputCoo1: string = '46.00000000';
  private inputCoo2: string = '25.00000000';
  private inputCoo3: string = '100.000';

  private outputCoords: number[] | null = null;
  private outputWarning: string | null = null;
  private quasigeoidZeta: number | null = null;
  private isCalculating: boolean = false;
  private pendingCalculation: boolean = false;
  private qualityStatus: 'ready' | 'valid' | 'invalid_input' | 'invalid_angle' | 'out_of_grid' | 'no_data' = 'ready';

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

  public setDirection(op: TransformationOp): void {
    if (this.op === op) return;
    this.op = op;

    const hasValidOutput = this.outputCoords !== null && !this.outputWarning;

    if (op === 'Stereo70ToETRS89') {
      if (hasValidOutput && this.outputCoords && this.outputCoords.length >= 2) {
        // Output from ETRS89ToStereo70: [N, E, (H)]
        const nStr = this.outputCoords[0].toFixed(3);
        const eStr = this.outputCoords[1].toFixed(3);
        this.inputCoo1 = this.order === 'EN' ? eStr : nStr;
        this.inputCoo2 = this.order === 'EN' ? nStr : eStr;
        if (this.outputCoords[2] !== undefined) {
          this.inputCoo3 = this.outputCoords[2].toFixed(3);
        }
      } else {
        // Fallback default: Stereo 70 projection center origin (500000, 500000, 100)
        this.inputCoo1 = '500000.000';
        this.inputCoo2 = '500000.000';
        if (!this.inputCoo3 || this.inputCoo3.trim() === '') {
          this.inputCoo3 = '100.000';
        }
      }
    } else {
      // Switched to ETRS89ToStereo70
      if (hasValidOutput && this.outputCoords && this.outputCoords.length >= 2) {
        // Output from Stereo70ToETRS89: [Lat, Lon, (h)] (in degrees)
        const lat = this.outputCoords[0];
        const lon = this.outputCoords[1];
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
        this.inputCoo1 = this.order === 'EN' ? lonStr : latStr;
        this.inputCoo2 = this.order === 'EN' ? latStr : lonStr;
        if (this.outputCoords[2] !== undefined) {
          this.inputCoo3 = this.outputCoords[2].toFixed(3);
        }
      } else {
        // Fallback default: projection center origin (46° N, 25° E, 100 m)
        let latStr = '46.00000000';
        let lonStr = '25.00000000';
        if (this.angularUnit === 'dms') {
          latStr = '46°00\'00.00000"';
          lonStr = '25°00\'00.00000"';
        } else if (this.angularUnit === 'radians') {
          latStr = degToRad(46.0).toFixed(10);
          lonStr = degToRad(25.0).toFixed(10);
        }
        this.inputCoo1 = this.order === 'EN' ? lonStr : latStr;
        this.inputCoo2 = this.order === 'EN' ? latStr : lonStr;
        if (!this.inputCoo3 || this.inputCoo3.trim() === '') {
          this.inputCoo3 = '100.000';
        }
      }
    }
    this.render();
    this.calculate();
  }

  public getDirection(): TransformationOp {
    return this.op;
  }

  public setAngularUnit(unit: AngularUnitType): void {
    if (this.angularUnit === unit) return;
    const oldUnit = this.angularUnit;
    this.angularUnit = unit;
    if (this.op === 'ETRS89ToStereo70') {
      try {
        const c1Val = this.order === 'NE' ? this.inputCoo1 : this.inputCoo2;
        const c2Val = this.order === 'NE' ? this.inputCoo2 : this.inputCoo1;

        const parseOldToDeg = (val: string): number => {
          if (oldUnit === 'radians') {
            if (val.includes('°') || val.includes("'") || val.includes('"')) {
              return parseAngleToDeg(val);
            }
            const num = Number(val.replace(/rad/gi, '').trim());
            return !isNaN(num) ? radToDeg(num) : NaN;
          }
          return parseAngleToDeg(val);
        };

        const degLat = parseOldToDeg(c1Val);
        const degLon = parseOldToDeg(c2Val);
        if (!isNaN(degLat) && !isNaN(degLon)) {
          const latStr = unit === 'dms' ? formatDegToDMS(degLat) : unit === 'radians' ? degToRad(degLat).toFixed(10) : degLat.toFixed(8);
          const lonStr = unit === 'dms' ? formatDegToDMS(degLon) : unit === 'radians' ? degToRad(degLon).toFixed(10) : degLon.toFixed(8);
          this.inputCoo1 = this.order === 'NE' ? latStr : lonStr;
          this.inputCoo2 = this.order === 'NE' ? lonStr : latStr;
        }
      } catch {
        // keep existing if parse fails
      }
    }
    this.render();
  }

  public setCoordinateOrder(order: CoordinateOrder): void {
    if (this.order === order) return;
    this.order = order;
    // Swap coordinate values so inputs match the newly selected column positions
    const tmp = this.inputCoo1;
    this.inputCoo1 = this.inputCoo2;
    this.inputCoo2 = tmp;
    this.render();
    this.calculate();
  }

  public getCoordinateOrder(): CoordinateOrder {
    return this.order;
  }

  public clearCoordinates(): void {
    this.inputCoo1 = '';
    this.inputCoo2 = '';
    this.inputCoo3 = '';
    this.outputCoords = null;
    this.outputWarning = null;
    this.quasigeoidZeta = null;
    this.qualityStatus = 'ready';
    this.updateQueryParams();
    this.callbacks?.onPointTransformed?.({
      lat: NaN,
      lon: NaN,
      outOfGrid: false,
    });
    this.render();
  }

  public setQualityStatus(status: 'ready' | 'valid' | 'invalid_input' | 'invalid_angle' | 'out_of_grid' | 'no_data'): void {
    this.qualityStatus = status;
    this.render();
  }

  public setCoordinatesFromMap(lat: number, lon: number): void {
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
      this.inputCoo1 = lonStr;
      this.inputCoo2 = latStr;
    } else {
      this.inputCoo1 = latStr;
      this.inputCoo2 = lonStr;
    }

    this.render();
    this.calculate();
  }

  public setStereoCoordinatesFromMap(northing: number, easting: number): void {
    const nStr = northing.toFixed(3);
    const eStr = easting.toFixed(3);
    if (this.order === 'EN') {
      this.inputCoo1 = eStr;
      this.inputCoo2 = nStr;
    } else {
      this.inputCoo1 = nStr;
      this.inputCoo2 = eStr;
    }
    this.render();
    this.calculate();
  }

  private cleanCoordinateString(val: string): string {
    let s = val.trim();
    // Strip wrapping quotes
    if ((s.startsWith('"') && s.endsWith('"')) || (s.startsWith("'") && s.endsWith("'"))) {
      s = s.slice(1, -1).trim();
    }
    // Remove trailing unit words/letters like 'm', 'meters', 'metri'
    s = s.replace(/(\d+(?:\.\d+)?)\s*(?:m|metri|meters)\b/gi, '$1');
    // Strip compass letters (N, S, E, W) unless part of DMS string with minute/sec symbols
    if (!s.includes('°') && !s.includes("'") && !s.includes('"')) {
      s = s.replace(/([0-9.-]+)\s*[NSEWnsew]\b/g, '$1');
      s = s.replace(/^[NSEWnsew]\s*([0-9.-]+)/g, '$1');
    } else {
      // Normalize DMS: remove spaces after degree and minute symbols
      s = s.replace(/°\s+/g, '°').replace(/'\s+/g, "'");
    }
    return s.trim();
  }

  private handlePasteIntoInput(e: ClipboardEvent, targetIndex: 1 | 2 | 3): void {
    const text = e.clipboardData?.getData('text');
    if (!text) return;

    const trimmed = text.trim();
    // If the entire pasted text is a single DMS string, do not treat whitespace inside it as a multi-coordinate delimiter
    const isSingleDmsString = /^[-+]?\d{1,3}°\s*\d{1,2}'\s*[\d.]+"?\s*[NSEW]?$/i.test(trimmed);

    const delimiterMatch = isSingleDmsString
      ? null
      : trimmed.includes('\t')
      ? '\t'
      : trimmed.includes(';')
      ? ';'
      : trimmed.includes(',')
      ? ','
      : /\s{2,}|\s+/.test(trimmed)
      ? /\s+/
      : null;

    if (delimiterMatch && targetIndex === 1) {
      const parts = trimmed.split(delimiterMatch).map((p) => p.trim()).filter(Boolean);
      if (parts.length >= 2) {
        e.preventDefault();
        this.inputCoo1 = this.cleanCoordinateString(parts[0]);
        this.inputCoo2 = this.cleanCoordinateString(parts[1]);
        if (parts.length >= 3) {
          this.inputCoo3 = this.cleanCoordinateString(parts[2]);
        }
        this.render();
        this.calculate();
        return;
      }
    }

    const cleaned = this.cleanCoordinateString(trimmed);
    if (cleaned !== text) {
      e.preventDefault();
      document.execCommand('insertText', false, cleaned);
    }
  }

  private readQueryParams(): void {
    const params = new URLSearchParams(window.location.search);
    if (params.has('op')) {
      const opParam = params.get('op');
      if (opParam === 'Stereo70ToETRS89' || opParam === 'ETRS89ToStereo70') {
        this.op = opParam;
        if (opParam === 'Stereo70ToETRS89' && !params.has('c1') && !params.has('c2')) {
          this.inputCoo1 = '500000.000';
          this.inputCoo2 = '500000.000';
        }
      }
    }
    if (params.has('unit')) {
      const unitParam = params.get('unit') as AngularUnitType;
      if (unitParam === 'degrees' || unitParam === 'dms' || unitParam === 'radians') {
        this.angularUnit = unitParam;
      }
    }
    if (params.has('order')) {
      const orderParam = params.get('order') as CoordinateOrder;
      if (orderParam === 'NE' || orderParam === 'EN') {
        this.order = orderParam;
      }
    }
    if (params.has('c1')) this.inputCoo1 = params.get('c1') || this.inputCoo1;
    if (params.has('c2')) this.inputCoo2 = params.get('c2') || this.inputCoo2;
    if (params.has('c3')) this.inputCoo3 = params.get('c3') || this.inputCoo3;
  }

  private updateQueryParams(): void {
    const params = new URLSearchParams();
    params.set('op', this.op);
    params.set('unit', this.angularUnit);
    params.set('order', this.order);
    params.set('c1', this.inputCoo1);
    params.set('c2', this.inputCoo2);
    if (this.inputCoo3.trim() !== '') {
      params.set('c3', this.inputCoo3);
    }
    const newUrl = `${window.location.pathname}?${params.toString()}`;
    window.history.replaceState({}, '', newUrl);
  }

  private async calculate(): Promise<void> {
    if (this.isCalculating) {
      this.pendingCalculation = true;
      return;
    }
    this.isCalculating = true;
    this.pendingCalculation = false;

    try {
      let raw1: number;
      let raw2: number;
      let raw3: number | undefined;

      const val1 = this.inputCoo1.trim();
      const val2 = this.inputCoo2.trim();

      if (!val1 && !val2) {
        this.outputCoords = null;
        this.outputWarning = null;
        this.quasigeoidZeta = null;
        this.qualityStatus = 'ready';
        return;
      }

      // Determine northing/lat vs easting/lon according to coordinate order
      const northingOrLat = this.order === 'NE' ? val1 : val2;
      const eastingOrLon = this.order === 'NE' ? val2 : val1;

      if (this.op === 'Stereo70ToETRS89') {
        const nonNumericAngleRegex = /[°'"’“”a-zA-Z]/;
        if (
          nonNumericAngleRegex.test(northingOrLat) ||
          nonNumericAngleRegex.test(eastingOrLon) ||
          !isValidDecimal(northingOrLat) ||
          !isValidDecimal(eastingOrLon)
        ) {
          this.outputCoords = null;
          this.outputWarning = null;
          this.quasigeoidZeta = null;
          this.qualityStatus = 'invalid_input';
          return;
        }

        raw1 = Number(northingOrLat);
        raw2 = Number(eastingOrLon);
      } else {
        // ETRS89 is source: parse angle into decimal degrees
        const parseInputAngle = (str: string): number => {
          if (this.angularUnit === 'radians') {
            if (str.includes('°') || str.includes("'") || str.includes('"')) {
              return parseAngleToDeg(str);
            }
            const clean = str.replace(/rad/gi, '').trim();
            const num = Number(clean);
            return !isNaN(num) ? radToDeg(num) : NaN;
          }
          return parseAngleToDeg(str);
        };

        try {
          raw1 = parseInputAngle(northingOrLat);
        } catch {
          raw1 = NaN;
        }
        try {
          raw2 = parseInputAngle(eastingOrLon);
        } catch {
          raw2 = NaN;
        }

        if (isNaN(raw1) || isNaN(raw2)) {
          this.outputCoords = null;
          this.outputWarning = null;
          this.quasigeoidZeta = null;
          this.qualityStatus = 'invalid_input';
          return;
        }

        // Dedicated check for physically impossible angles: Latitude [-90, 90], Longitude [-180, 180]
        if (raw1 < -90 || raw1 > 90 || raw2 < -180 || raw2 > 180) {
          this.outputCoords = null;
          this.outputWarning = null;
          this.quasigeoidZeta = null;
          this.qualityStatus = 'invalid_angle';
          return;
        }
      }

      // Auto-detect 2D vs 3D based on 3rd coordinate input
      const trimmed3 = this.inputCoo3.trim();
      const is3D = trimmed3 !== '';

      if (is3D) {
        if (!isValidDecimal(trimmed3)) {
          this.outputCoords = null;
          this.outputWarning = null;
          this.quasigeoidZeta = null;
          this.qualityStatus = 'invalid_input';
          return;
        }
        raw3 = Number(trimmed3);
      }

      const coos = is3D && raw3 !== undefined ? [raw1, raw2, raw3] : [raw1, raw2];

      // Coarse pre-flight envelope check
      const isCoarseInGrid =
        this.op === 'Stereo70ToETRS89'
          ? isWithinStereo70Grid(raw1, raw2)
          : isWithinETRS89Grid(raw1, raw2);

      let res: { coos: number[]; warning?: string | null };
      if (!isCoarseInGrid) {
        res = { coos, warning: 'Out of grid' };
      } else {
        res = await apiClient.transformPoint({
          op: this.op,
          coos,
          unit: 'degrees',
        });
      }

      this.outputCoords = res.coos;
      this.outputWarning = res.warning || null;
      const isOutOfGrid = res.warning === 'Out of grid';
      const isNoData = res.warning === 'No data on grid';
      const isUntransformable = isOutOfGrid || isNoData;

      if (isOutOfGrid) {
        this.qualityStatus = 'out_of_grid';
      } else if (isNoData) {
        this.qualityStatus = 'no_data';
      } else {
        this.qualityStatus = 'valid';
      }

      // Quasigeoid height calculation in 3D:
      // When Stereo70 -> ETRS89: res.coos[2] is h, input is H => zeta = h - H
      // When ETRS89 -> Stereo70: res.coos[2] is H, input is h => zeta = h - H
      if (!isUntransformable && is3D && this.outputCoords.length >= 3 && raw3 !== undefined) {
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
      let mapLat = NaN;
      let mapLon = NaN;
      let stereoN: number | undefined;
      let stereoE: number | undefined;
      let elev: number | undefined;
      let altitude: number | undefined;
      let ellipsoidalHeight: number | undefined;

      if (this.op === 'Stereo70ToETRS89') {
        stereoN = this.order === 'EN' ? raw2 : raw1;
        stereoE = this.order === 'EN' ? raw1 : raw2;
        altitude = raw3;
        if (!isUntransformable && this.outputCoords) {
          mapLat = this.outputCoords[0];
          mapLon = this.outputCoords[1];
          elev = this.outputCoords[2];
          ellipsoidalHeight = this.outputCoords[2];
        }
      } else {
        mapLat = raw1;
        mapLon = raw2;
        ellipsoidalHeight = raw3;
        if (!isUntransformable && this.outputCoords) {
          stereoN = this.outputCoords[0];
          stereoE = this.outputCoords[1];
          elev = this.outputCoords[2];
          altitude = this.outputCoords[2];
        }
      }

      this.callbacks.onPointTransformed({
        lat: mapLat,
        lon: mapLon,
        outOfGrid: isUntransformable,
        pointId: 'Single Point',
        stereoN,
        stereoE,
        elev,
        altitude,
        ellipsoidalHeight,
        zeta: this.quasigeoidZeta,
        angularUnit: this.angularUnit,
      });

      if (isOutOfGrid) {
        showToast(t('point.outOfGridWarn'), 'warning');
      }
    } catch (err: any) {
      this.outputCoords = null;
      this.outputWarning = err.message || 'Error';
      this.qualityStatus = 'invalid_input';
      showToast(err.message || 'Error calculating coordinates', 'danger');
    } finally {
      this.isCalculating = false;
      this.render();
      if (this.pendingCalculation) {
        this.pendingCalculation = false;
        this.calculate();
      }
    }
  }

  private render(): void {
    const isS70 = this.op === 'Stereo70ToETRS89';
    let c1Label: string;
    let c2Label: string;

    if (isS70) {
      c1Label = this.order === 'EN' ? t('point.eastingLabel') : t('point.northingLabel');
      c2Label = this.order === 'EN' ? t('point.northingLabel') : t('point.eastingLabel');
    } else {
      c1Label = this.order === 'EN' ? t('point.longitudeLabel') : t('point.latitudeLabel');
      c2Label = this.order === 'EN' ? t('point.latitudeLabel') : t('point.longitudeLabel');
    }
    const c3Label = isS70 ? t('point.altitudeLabel') : t('point.heightLabel');

    const isTargetGeo = this.op === 'Stereo70ToETRS89';
    let outCol1Val = '';
    let outCol2Val = '';
    let val3Str = '';
    let outCol1Raw = '';
    let outCol2Raw = '';
    let raw3Str = '';

    const isOutOfGrid = this.outputWarning === 'Out of grid';
    const isNoData = this.outputWarning === 'No data on grid';
    const isUntransformable = isOutOfGrid || isNoData;
    const is3D = this.outputCoords !== null && this.outputCoords.length >= 3;

    if (this.outputCoords && !isUntransformable) {
      if (isTargetGeo) {
        const lat = this.outputCoords[0];
        const lon = this.outputCoords[1];

        if (this.angularUnit === 'dms') {
          outCol1Raw = this.order === 'EN' ? formatDegToDMS(lon) : formatDegToDMS(lat);
          outCol2Raw = this.order === 'EN' ? formatDegToDMS(lat) : formatDegToDMS(lon);
        } else if (this.angularUnit === 'radians') {
          outCol1Raw = this.order === 'EN' ? degToRad(lon).toFixed(10) : degToRad(lat).toFixed(10);
          outCol2Raw = this.order === 'EN' ? degToRad(lat).toFixed(10) : degToRad(lon).toFixed(10);
        } else {
          outCol1Raw = this.order === 'EN' ? lon.toFixed(8) : lat.toFixed(8);
          outCol2Raw = this.order === 'EN' ? lat.toFixed(8) : lon.toFixed(8);
        }
        outCol1Val = outCol1Raw;
        outCol2Val = outCol2Raw;
      } else {
        const n = this.outputCoords[0];
        const e = this.outputCoords[1];
        outCol1Raw = this.order === 'EN' ? e.toFixed(3) : n.toFixed(3);
        outCol2Raw = this.order === 'EN' ? n.toFixed(3) : e.toFixed(3);
        outCol1Val = outCol1Raw;
        outCol2Val = outCol2Raw;
      }

      if (this.outputCoords[2] !== undefined) {
        raw3Str = this.outputCoords[2].toFixed(3);
        val3Str = raw3Str;
      }
    }

    let outCol1Label: string;
    let outCol2Label: string;
    if (isTargetGeo) {
      outCol1Label = this.order === 'EN' ? t('point.longitudeLabel') : t('point.latitudeLabel');
      outCol2Label = this.order === 'EN' ? t('point.latitudeLabel') : t('point.longitudeLabel');
    } else {
      outCol1Label = this.order === 'EN' ? t('point.eastingLabel') : t('point.northingLabel');
      outCol2Label = this.order === 'EN' ? t('point.northingLabel') : t('point.eastingLabel');
    }
    const label3 = isTargetGeo ? t('point.heightLabel') : t('point.altitudeLabel');

    let zetaDisplay = '';
    if (this.quasigeoidZeta !== null) {
      zetaDisplay = (this.quasigeoidZeta >= 0 ? '+' : '') + this.quasigeoidZeta.toFixed(3);
    }

    // Map point coordinates
    let mapLat = NaN;
    let mapLon = NaN;
    if (this.outputCoords) {
      if (this.op === 'Stereo70ToETRS89') {
        if (!isUntransformable) {
          mapLat = this.outputCoords[0];
          mapLon = this.outputCoords[1];
        }
      } else {
        const val1 = this.inputCoo1.trim();
        const val2 = this.inputCoo2.trim();
        const northingOrLat = this.order === 'NE' ? val1 : val2;
        const eastingOrLon = this.order === 'NE' ? val2 : val1;
        const parseInputAngle = (str: string): number => {
          if (this.angularUnit === 'radians') {
            if (str.includes('°') || str.includes("'") || str.includes('"')) {
              return parseAngleToDeg(str);
            }
            const clean = str.replace(/rad/gi, '').trim();
            const num = Number(clean);
            return !isNaN(num) ? radToDeg(num) : NaN;
          }
          return parseAngleToDeg(str);
        };
        try {
          mapLat = parseInputAngle(northingOrLat);
          mapLon = parseInputAngle(eastingOrLon);
        } catch {
          mapLat = NaN;
          mapLon = NaN;
        }
      }
    }

    const canCenterOnMap = !isUntransformable && !isNaN(mapLat) && !isNaN(mapLon) && isInsideRomaniaGrid(mapLat, mapLon);

    let transformQualityBadge = '';
    switch (this.qualityStatus) {
      case 'valid':
        transformQualityBadge = `<span class="status-pill success">${t('quality.valid')}</span>`;
        break;
      case 'invalid_input':
        transformQualityBadge = `<span class="status-pill danger">${t('quality.invalidInput')}</span>`;
        break;
      case 'invalid_angle':
        transformQualityBadge = `<span class="status-pill danger">${t('quality.invalidAngle')}</span>`;
        break;
      case 'out_of_grid':
        transformQualityBadge = `<span class="status-pill warning">${t('quality.outsideGrid')}</span>`;
        break;
      case 'no_data':
        transformQualityBadge = `<span class="status-pill warning">${t('quality.noGridData')}</span>`;
        break;
      case 'ready':
      default:
        transformQualityBadge = `<span class="status-pill neutral">${t('quality.readyToTransform')}</span>`;
        break;
    }

    this.container.innerHTML = `
      <!-- Card 1: Input Section -->
      <div class="studio-card input-card" role="group" aria-label="${t('point.inputCoordsHint')}">
        <div class="section-micro-label">${t('point.inputCoordsHint')}</div>

        <div class="coord-grid-3">
          <div class="form-group margin-bottom-none">
            <label class="form-label">${c1Label}</label>
            <input 
              type="text" 
              class="form-input mono" 
              id="pt-input-1" 
              aria-label="${c1Label}" 
              value="${escapeHtmlAttr(this.inputCoo1)}" 
              autocomplete="off" 
              autocorrect="off" 
              autocapitalize="none" 
              spellcheck="false" 
            />
          </div>
          <div class="form-group margin-bottom-none">
            <label class="form-label">${c2Label}</label>
            <input 
              type="text" 
              class="form-input mono" 
              id="pt-input-2" 
              aria-label="${c2Label}" 
              value="${escapeHtmlAttr(this.inputCoo2)}" 
              autocomplete="off" 
              autocorrect="off" 
              autocapitalize="none" 
              spellcheck="false" 
            />
          </div>
          <div class="form-group margin-bottom-none">
            <label class="form-label">${c3Label}</label>
            <input 
              type="text" 
              class="form-input mono" 
              id="pt-input-3" 
              aria-label="${c3Label}"
              value="${escapeHtmlAttr(this.inputCoo3)}" 
              placeholder="${t('point.optionalPlaceholder')}" 
              autocomplete="off" 
              autocorrect="off" 
              autocapitalize="none" 
              spellcheck="false" 
            />
          </div>
        </div>
      </div>

      <!-- Active Separator Button between Input & Output -->
      <div class="transform-action-row">
        <button class="btn btn-primary" id="pt-calculate-btn">
          <span>${this.isCalculating ? t('point.calculating') : t('point.transformBtn')}</span>
        </button>
      </div>

      <!-- Card 2: Output Section -->
      <div class="studio-card result-card" role="group" aria-label="${t('point.outputSectionLabel')}">
        <div class="flex-between-center margin-bottom-sm">
          <div class="section-micro-label">${t('point.outputSectionLabel')}</div>
          <div id="pt-copy-btn-slot" data-c1="${escapeHtmlAttr(outCol1Raw)}" data-c2="${escapeHtmlAttr(outCol2Raw)}" data-c3="${escapeHtmlAttr(raw3Str)}"></div>
        </div>

        <div class="output-grid-3">
          <div class="form-group margin-bottom-none">
            <label class="form-label">${outCol1Label}</label>
            <input type="text" readonly class="form-input mono" aria-label="${outCol1Label}" value="${escapeHtmlAttr(outCol1Val)}" />
          </div>
          <div class="form-group margin-bottom-none">
            <label class="form-label">${outCol2Label}</label>
            <input type="text" readonly class="form-input mono" aria-label="${outCol2Label}" value="${escapeHtmlAttr(outCol2Val)}" />
          </div>
          <div class="form-group margin-bottom-none">
            <label class="form-label">${label3}</label>
            <input type="text" readonly class="form-input mono" aria-label="${label3}" value="${escapeHtmlAttr(val3Str)}" />
          </div>
        </div>

        <!-- Quasi-geoid Height Info Archetype -->
        <div class="output-quasigeoid-row">
          <span class="info-label">${t('point.quasigeoidLabel')}</span>
          <span class="info-field-badge">${escapeHtmlAttr(zetaDisplay)}</span>
        </div>
      </div>

      <!-- Card 3: Report Section -->
      <div class="studio-card report-card" role="group" aria-label="${t('quality.sectionTitle')}">
        <div class="section-micro-label">${t('quality.sectionTitle')}</div>
        <div class="quality-report-content flex-column-gap-sm">
          <div class="quality-report-row">
            <span class="quality-label">${t('quality.sectionCoords')}</span>
            <div class="quality-report-badges">
              <span class="status-pill status-pill-neutral">
                ${is3D ? t('quality.dim3D') : t('quality.dim2D')}
              </span>
            </div>
          </div>

          <div class="quality-report-divider-dashed"></div>

          <div class="quality-report-row">
            <span class="quality-label">${t('quality.sectionTransform')}</span>
            <div class="quality-report-badges">
              ${transformQualityBadge}
            </div>
          </div>
        </div>
      </div>

      <!-- Decoupled Map Navigation Dock -->
      <div class="workbench-actions-dock">
        <button 
          class="btn btn-outline center-map-btn btn-secondary-action" 
          id="pt-center-map-btn" 
          ${canCenterOnMap ? '' : 'disabled'}
          data-lat="${mapLat}"
          data-lon="${mapLon}"
          data-oog="${isOutOfGrid ? '1' : '0'}"
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
          id="pt-zoom-extents-btn" 
          ${canCenterOnMap ? '' : 'disabled'}
          data-lat="${mapLat}"
          data-lon="${mapLon}"
          data-oog="${isOutOfGrid ? '1' : '0'}"
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
  }

  private attachEventListeners(): void {
    const i1 = this.container.querySelector('#pt-input-1') as HTMLInputElement;
    const i2 = this.container.querySelector('#pt-input-2') as HTMLInputElement;
    const i3 = this.container.querySelector('#pt-input-3') as HTMLInputElement;

    const onInputChange = () => {
      if (i1) this.inputCoo1 = i1.value.trim();
      if (i2) this.inputCoo2 = i2.value.trim();
      if (i3) this.inputCoo3 = i3.value.trim();
    };

    [i1, i2, i3].forEach((input, idx) => {
      if (input) {
        input.addEventListener('input', onInputChange);
        input.addEventListener('paste', (e) => this.handlePasteIntoInput(e, (idx + 1) as 1 | 2 | 3));
        input.addEventListener('keydown', (e) => {
          if (e.key === 'Enter') {
            onInputChange();
            this.calculate();
          }
        });
      }
    });

    // Auto-select all text on focus on all inputs
    this.container.querySelectorAll('input.form-input').forEach((el) => {
      const input = el as HTMLInputElement;
      input.addEventListener('focus', () => input.select());
    });

    this.container.querySelector('#pt-calculate-btn')?.addEventListener('click', () => {
      onInputChange();
      this.calculate();
    });

    // Top-right code-block copy button using reusable UI factory
    const copySlot = this.container.querySelector('#pt-copy-btn-slot');
    if (copySlot) {
      const c1 = copySlot.getAttribute('data-c1') || '';
      const c2 = copySlot.getAttribute('data-c2') || '';
      const c3 = copySlot.getAttribute('data-c3') || '';
      const isOutOfGrid = this.outputWarning === 'Out of grid';

      const copyBtn = renderCopyButton({
        id: 'pt-copy-code-btn',
        title: t('point.copyTuple'),
        ariaLabel: t('point.copyTuple'),
        disabled: !this.outputCoords || isOutOfGrid,
        feedbackText: t('point.copiedTooltip'),
        getCopyText: () => {
          const parts = [c1, c2].filter(Boolean);
          if (parts.length === 0) return '';
          if (c3) parts.push(c3);
          return parts.join(', ');
        },
        onCopied: () => {
          showToast(t('point.copiedToast'), 'success');
        },
      });
      copySlot.replaceWith(copyBtn);
    }

    // Center on Map action button
    const centerBtn = this.container.querySelector('#pt-center-map-btn') as HTMLButtonElement;
    if (centerBtn) {
      centerBtn.addEventListener('click', () => {
        const lat = parseFloat(centerBtn.getAttribute('data-lat') || '0');
        const lon = parseFloat(centerBtn.getAttribute('data-lon') || '0');
        const outOfGrid = centerBtn.getAttribute('data-oog') === '1';
        this.callbacks.onCenterOnMap?.({ lat, lon, outOfGrid });
      });
    }

    // Zoom to Extents action button
    const zoomExtentsBtn = this.container.querySelector('#pt-zoom-extents-btn') as HTMLButtonElement;
    if (zoomExtentsBtn) {
      zoomExtentsBtn.addEventListener('click', () => {
        const lat = parseFloat(zoomExtentsBtn.getAttribute('data-lat') || '0');
        const lon = parseFloat(zoomExtentsBtn.getAttribute('data-lon') || '0');
        const outOfGrid = zoomExtentsBtn.getAttribute('data-oog') === '1';
        this.callbacks.onZoomToExtents?.({ lat, lon, outOfGrid });
      });
    }
  }
}
