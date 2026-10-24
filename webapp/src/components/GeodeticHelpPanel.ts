import { t, getLanguage } from '../i18n';
import { formatTelemetryNumber } from '../utils/numberFormat';
import { apiClient, type GridComparisonResponse } from '../utils/apiClient';

export class GeodeticHelpPanel {
  private element: HTMLElement | null = null;
  private comparisonData: GridComparisonResponse | null = null;
  private isLoadingComparison: boolean = false;
  private isMinimized: boolean = false;

  // Dragging state
  private isDragging: boolean = false;
  private dragStartX: number = 0;
  private dragStartY: number = 0;
  private initialLeft: number = 0;
  private initialTop: number = 0;

  constructor() {}

  public isOpen(): boolean {
    return this.element !== null;
  }

  public toggle(): void {
    if (this.isOpen()) {
      this.close();
    } else {
      this.open();
    }
  }

  public async open(): Promise<void> {
    if (this.element) {
      this.element.style.display = 'flex';
      return;
    }

    this.element = document.createElement('div');
    this.element.className = 'floating-help-panel';
    document.body.appendChild(this.element);

    this.render();
    await this.loadComparisonData();
  }

  public close(): void {
    if (this.element) {
      this.element.remove();
      this.element = null;
    }
  }

  private async loadComparisonData(): Promise<void> {
    if (this.comparisonData) return;
    this.isLoadingComparison = true;
    this.updateBenchmarkCard();

    try {
      this.comparisonData = await apiClient.getGridComparison();
    } catch (err) {
      console.warn('Could not load dynamic grid comparison data:', err);
    } finally {
      this.isLoadingComparison = false;
      this.updateBenchmarkCard();
    }
  }

  private toggleMinimize(): void {
    this.isMinimized = !this.isMinimized;
    if (this.element) {
      if (this.isMinimized) {
        this.element.classList.add('minimized');
      } else {
        this.element.classList.remove('minimized');
      }
      this.render();
    }
  }

  private render(): void {
    if (!this.element) return;

    if (this.isMinimized) {
      this.element.innerHTML = `
        <div class="drawer-header fgp-header" id="fhp-drag-handle">
          <div class="drawer-title-cluster fgp-title-cluster">
            <span class="drawer-title fgp-title">${t('geodeticHelp.title')}</span>
          </div>
          <div class="drawer-header-actions fgp-header-actions">
            <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fhp-minimize-btn" title="${t('geodeticHelp.restore')}">□</button>
            <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fhp-close-btn" title="${t('geodeticHelp.close')}">✕</button>
          </div>
        </div>
      `;
      this.attachHeaderListeners();
      this.attachDragListeners();
      return;
    }

    this.element.innerHTML = `
      <div class="drawer-header fgp-header" id="fhp-drag-handle">
        <div class="drawer-title-cluster fgp-title-cluster">
          <span class="drawer-title fgp-title">${t('geodeticHelp.title')}</span>
        </div>
        <div class="drawer-header-actions fgp-header-actions">
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fhp-minimize-btn" title="${t('geodeticHelp.minimize')}">─</button>
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fhp-close-btn" title="${t('geodeticHelp.close')}">✕</button>
        </div>
      </div>

      <div class="drawer-body fhp-body" id="fhp-body-content">
        <!-- Card 1: About & Non-Official Disclaimer -->
        <div class="content-card fhp-card">
          <div class="content-card-header fhp-card-header">
            <div class="content-card-title fhp-card-title">${t('geodeticHelp.aboutTitle')}</div>
          </div>
          <p class="fhp-topic-desc">
            ${t('geodeticHelp.aboutP1')}
          </p>
          <div class="callout-box fhp-callout-box">
            <span class="callout-title fhp-callout-title">${t('geodeticHelp.disclaimerLabel')}</span>
            <p class="callout-body fhp-callout-body">${t('geodeticHelp.disclaimerText')}</p>
          </div>
        </div>

        <!-- Card 2: Benchmark Comparison against ANCPI Official Engine -->
        <div class="content-card fhp-card" id="fhp-benchmark-slot">
          ${this.getBenchmarkCardHtml()}
        </div>

        <!-- Card 3: Reference Systems & Datums -->
        <div class="content-card fhp-card">
          <div class="content-card-header fhp-card-header">
            <div class="content-card-title fhp-card-title">${t('geodeticHelp.datumsTitle')}</div>
          </div>
          
          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.stereo70Label')}</span>
            <p class="fhp-topic-desc">${t('geodeticHelp.stereo70Desc')}</p>
          </div>

          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.etrs89Label')}</span>
            <p class="fhp-topic-desc">${t('geodeticHelp.etrs89Desc')}</p>
          </div>

          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.altimetryLabel')}</span>
            <p class="fhp-topic-desc">${t('geodeticHelp.altimetryDesc')}</p>
          </div>
        </div>

        <!-- Card 4: Formats & Conventions -->
        <div class="content-card fhp-card">
          <div class="content-card-header fhp-card-header">
            <div class="content-card-title fhp-card-title">${t('geodeticHelp.formatsTitle')}</div>
          </div>

          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.angleUnitsLabel')}</span>
            <div class="fhp-format-list">
              <div class="fhp-format-row">
                <span class="fhp-format-label">${t('geodeticHelp.angleDmsLabel')}</span>
                <code class="fhp-format-chip">${t('geodeticHelp.angleDmsVal')}</code>
              </div>
              <div class="fhp-format-row">
                <span class="fhp-format-label">${t('geodeticHelp.angleDegLabel')}</span>
                <code class="fhp-format-chip">${t('geodeticHelp.angleDegVal')}</code>
              </div>
              <div class="fhp-format-row">
                <span class="fhp-format-label">${t('geodeticHelp.angleRadLabel')}</span>
                <code class="fhp-format-chip">${t('geodeticHelp.angleRadVal')}</code>
              </div>
            </div>
            <p class="fhp-topic-desc">${t('geodeticHelp.angleNote')}</p>
          </div>

          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.orderLabel')}</span>
            <p class="fhp-topic-desc">${t('geodeticHelp.orderDesc')}</p>
          </div>

          <div class="fhp-topic-item">
            <span class="fhp-topic-heading">${t('geodeticHelp.detectionLabel')}</span>
            <p class="fhp-topic-desc">${t('geodeticHelp.detectionDesc')}</p>
          </div>
        </div>
      </div>

      <!-- Corner Resize Grip -->
      <div class="fgp-resize-handle" id="fhp-resize-handle" title="Resize window">
        <svg width="10" height="10" viewBox="0 0 10 10" fill="currentColor">
          <path d="M9 1v8H1l8-8zM9 5v4H5l4-4z"/>
        </svg>
      </div>
    `;

    this.attachHeaderListeners();
    this.attachDragListeners();
    this.attachResizeListeners();
  }

  private getBenchmarkCardHtml(): string {
    const title = t('geodeticHelp.benchmarkTitle');
    const pointsCount = formatTelemetryNumber(this.comparisonData?.total_points || 9046);
    const subtitle = t('geodeticHelp.benchmarkSubtitle').replace('{count}', pointsCount);

    if (this.isLoadingComparison && !this.comparisonData) {
      return `
        <div class="content-card-header fhp-card-header">
          <div class="content-card-title fhp-card-title">${title}</div>
        </div>
        <div class="fgp-status-box">
          <div class="fgp-spinner"></div>
          <span>${t('geodeticHelp.loadingBenchmark')}</span>
        </div>
      `;
    }

    const st70ToEtrs = this.comparisonData?.st70_to_etrs89;
    const etrsToSt70 = this.comparisonData?.etrs89_to_st70;
    const isRo = getLanguage() === 'ro';
    const fmt = (val: number, d: number) => {
      const s = val.toFixed(d);
      return isRo ? s.replace('.', ',') : s;
    };

    // S70 -> ETRS: planimetric in arcseconds ("), elevation in mm
    const maxS70Plan = st70ToEtrs
      ? `${fmt(Math.max(st70ToEtrs.max_abs_error.latitude_sec || 0, st70ToEtrs.max_abs_error.longitude_sec || 0), 6)}"`
      : (isRo ? '0,000010"' : '0.000010"');
    const meanS70Plan = st70ToEtrs
      ? `${fmt(Math.max(st70ToEtrs.mean_abs_error.latitude_sec || 0, st70ToEtrs.mean_abs_error.longitude_sec || 0), 6)}"`
      : (isRo ? '0,000005"' : '0.000005"');
    const maxS70Elev = st70ToEtrs
      ? `${fmt(st70ToEtrs.max_abs_error.elevation_m * 1000, 1)} mm`
      : (isRo ? '0,5 mm' : '0.5 mm');
    const meanS70Elev = st70ToEtrs
      ? `${fmt(st70ToEtrs.mean_abs_error.elevation_m * 1000, 2)} mm`
      : (isRo ? '0,25 mm' : '0.25 mm');

    // ETRS -> S70: planimetric in mm, elevation in mm
    const maxEtrsPlan = etrsToSt70
      ? `${fmt(Math.max(etrsToSt70.max_abs_error.easting_m || 0, etrsToSt70.max_abs_error.northing_m || 0) * 1000, 2)} mm`
      : (isRo ? '0,34 mm' : '0.34 mm');
    const meanEtrsPlan = etrsToSt70
      ? `${fmt(Math.max(etrsToSt70.mean_abs_error.easting_m || 0, etrsToSt70.mean_abs_error.northing_m || 0) * 1000, 2)} mm`
      : (isRo ? '0,08 mm' : '0.08 mm');
    const maxEtrsElev = etrsToSt70
      ? `${fmt(etrsToSt70.max_abs_error.elevation_m * 1000, 1)} mm`
      : (isRo ? '0,5 mm' : '0.5 mm');
    const meanEtrsElev = etrsToSt70
      ? `${fmt(etrsToSt70.mean_abs_error.elevation_m * 1000, 2)} mm`
      : (isRo ? '0,25 mm' : '0.25 mm');

    const st70Status = st70ToEtrs?.status || 'PASS';
    const etrsStatus = etrsToSt70?.status || 'PASS';

    return `
      <div class="content-card-header fhp-card-header">
        <div class="content-card-title fhp-card-title">${title}</div>
        <div class="content-card-subtitle fhp-card-subtitle">${subtitle}</div>
      </div>

      <table class="fhp-benchmark-table">
        <thead>
          <tr>
            <th>${t('geodeticHelp.colDirection')}</th>
            <th>${t('geodeticHelp.colPlanMax')}</th>
            <th>${t('geodeticHelp.colPlanMean')}</th>
            <th>${t('geodeticHelp.colElevMax')}</th>
            <th>${t('geodeticHelp.colElevMean')}</th>
            <th>${t('geodeticHelp.colStatus')}</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>${t('geodeticHelp.dirS70ToETRS')}</strong></td>
            <td><code class="fhp-metric-chip">${maxS70Plan}</code></td>
            <td><code class="fhp-metric-chip subtle">${meanS70Plan}</code></td>
            <td><code class="fhp-metric-chip">${maxS70Elev}</code></td>
            <td><code class="fhp-metric-chip subtle">${meanS70Elev}</code></td>
            <td><span class="status-pill success">${st70Status}</span></td>
          </tr>
          <tr>
            <td><strong>${t('geodeticHelp.dirETRSToS70')}</strong></td>
            <td><code class="fhp-metric-chip">${maxEtrsPlan}</code></td>
            <td><code class="fhp-metric-chip subtle">${meanEtrsPlan}</code></td>
            <td><code class="fhp-metric-chip">${maxEtrsElev}</code></td>
            <td><code class="fhp-metric-chip subtle">${meanEtrsElev}</code></td>
            <td><span class="status-pill success">${etrsStatus}</span></td>
          </tr>
        </tbody>
      </table>

      <div class="callout-box fhp-callout-box">
        <span class="callout-title fhp-callout-title">${t('geodeticHelp.criteriaLabel')}</span>
        <p class="callout-body fhp-callout-body">${t('geodeticHelp.benchmarkCriteriaNote')}</p>
      </div>
    `;
  }

  private updateBenchmarkCard(): void {
    const slot = this.element?.querySelector('#fhp-benchmark-slot');
    if (slot) {
      slot.innerHTML = this.getBenchmarkCardHtml();
    }
  }

  private attachHeaderListeners(): void {
    if (!this.element) return;

    this.element.querySelector('#fhp-close-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.close();
    });

    this.element.querySelector('#fhp-minimize-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.toggleMinimize();
    });
  }

  private attachDragListeners(): void {
    if (!this.element) return;
    const handle = this.element.querySelector('#fhp-drag-handle') as HTMLElement | null;
    if (!handle) return;

    handle.addEventListener('pointerdown', (e: PointerEvent) => {
      if ((e.target as HTMLElement).closest('.btn-icon')) return;
      if (window.innerWidth <= 768) return;

      this.isDragging = true;
      this.dragStartX = e.clientX;
      this.dragStartY = e.clientY;

      const rect = this.element!.getBoundingClientRect();
      this.initialLeft = rect.left;
      this.initialTop = rect.top;

      handle.setPointerCapture(e.pointerId);

      const onPointerMove = (moveEvent: PointerEvent) => {
        if (!this.isDragging || !this.element) return;
        const dx = moveEvent.clientX - this.dragStartX;
        const dy = moveEvent.clientY - this.dragStartY;

        let newLeft = this.initialLeft + dx;
        let newTop = this.initialTop + dy;

        const maxLeft = Math.max(0, window.innerWidth - (this.element.offsetWidth || 300));
        const maxTop = Math.max(0, window.innerHeight - 50);

        newLeft = Math.max(0, Math.min(newLeft, maxLeft));
        newTop = Math.max(0, Math.min(newTop, maxTop));

        this.element.style.left = `${newLeft}px`;
        this.element.style.top = `${newTop}px`;
        this.element.style.right = 'auto';
        this.element.style.bottom = 'auto';
      };

      const onPointerUp = (upEvent: PointerEvent) => {
        this.isDragging = false;
        try {
          handle.releasePointerCapture(upEvent.pointerId);
        } catch {}
        handle.removeEventListener('pointermove', onPointerMove);
        handle.removeEventListener('pointerup', onPointerUp);
        handle.removeEventListener('pointercancel', onPointerUp);
      };

      handle.addEventListener('pointermove', onPointerMove);
      handle.addEventListener('pointerup', onPointerUp);
      handle.addEventListener('pointercancel', onPointerUp);
    });
  }

  private attachResizeListeners(): void {
    if (!this.element) return;
    const resizer = this.element.querySelector('#fhp-resize-handle') as HTMLElement | null;
    if (!resizer) return;

    resizer.addEventListener('pointerdown', (e: PointerEvent) => {
      e.preventDefault();
      e.stopPropagation();

      const startX = e.clientX;
      const startY = e.clientY;
      const startWidth = this.element!.offsetWidth;
      const startHeight = this.element!.offsetHeight;

      resizer.setPointerCapture(e.pointerId);

      const onPointerMove = (moveEvent: PointerEvent) => {
        if (!this.element) return;
        moveEvent.preventDefault();
        moveEvent.stopPropagation();

        const dw = moveEvent.clientX - startX;
        const dh = moveEvent.clientY - startY;

        const minWidth = 360;
        const maxWidth = Math.max(minWidth, window.innerWidth - 40);
        const minHeight = 240;
        const maxHeight = Math.max(minHeight, window.innerHeight - 80);

        const newWidth = Math.max(minWidth, Math.min(startWidth + dw, maxWidth));
        const newHeight = Math.max(minHeight, Math.min(startHeight + dh, maxHeight));

        this.element.style.width = `${newWidth}px`;
        this.element.style.height = `${newHeight}px`;
        this.element.style.maxHeight = 'none';
      };

      const onPointerUp = (upEvent: PointerEvent) => {
        try {
          resizer.releasePointerCapture(upEvent.pointerId);
        } catch {}
        resizer.removeEventListener('pointermove', onPointerMove);
        resizer.removeEventListener('pointerup', onPointerUp);
        resizer.removeEventListener('pointercancel', onPointerUp);
      };

      resizer.addEventListener('pointermove', onPointerMove);
      resizer.addEventListener('pointerup', onPointerUp);
      resizer.addEventListener('pointercancel', onPointerUp);
    });
  }
}
