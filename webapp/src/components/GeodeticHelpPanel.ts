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
        <div class="fgp-header" id="fhp-drag-handle">
          <div class="fgp-title-cluster">
            <svg class="fgp-header-icon" width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">
              <path d="M5.5 5.5a2.5 2.5 0 0 1 4.5 1.5c0 1.5-2 2-2 3v.5M8 13.5h.01"/>
            </svg>
            <span class="fgp-title">${t('geodeticHelp.title')}</span>
          </div>
          <div class="fgp-header-actions">
            <button class="btn-icon fgp-icon-btn" id="fhp-minimize-btn" title="Restore">□</button>
            <button class="btn-icon fgp-icon-btn" id="fhp-close-btn" title="Close">✕</button>
          </div>
        </div>
      `;
      this.attachHeaderListeners();
      this.attachDragListeners();
      return;
    }

    this.element.innerHTML = `
      <div class="fgp-header" id="fhp-drag-handle">
        <div class="fgp-title-cluster">
          <svg class="fgp-header-icon" width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">
            <path d="M5.5 5.5a2.5 2.5 0 0 1 4.5 1.5c0 1.5-2 2-2 3v.5M8 13.5h.01"/>
          </svg>
          <span class="fgp-title">${t('geodeticHelp.title')}</span>
        </div>
        <div class="fgp-header-actions">
          <button class="btn-icon fgp-icon-btn" id="fhp-minimize-btn" title="${t('geodeticHelp.minimize')}">─</button>
          <button class="btn-icon fgp-icon-btn" id="fhp-close-btn" title="${t('geodeticHelp.close')}">✕</button>
        </div>
      </div>

      <div class="fhp-body" id="fhp-body-content">
        <!-- Card 1: About & Non-Official Disclaimer -->
        <div class="fhp-card">
          <div class="fhp-card-title">${t('geodeticHelp.aboutTitle')}</div>
          <p style="margin: 0; color: var(--text-main); font-size: 12px;">
            ${t('geodeticHelp.aboutP1')}
          </p>
          <div class="fhp-disclaimer-box">
            <strong style="color: var(--text-main);">${t('geodeticHelp.disclaimerLabel')}</strong>
            ${t('geodeticHelp.disclaimerText')}
          </div>
        </div>

        <!-- Card 2: Benchmark Comparison against ANCPI Official Engine -->
        <div class="fhp-card" id="fhp-benchmark-slot">
          ${this.getBenchmarkCardHtml()}
        </div>

        <!-- Card 3: Reference Systems & Datums -->
        <div class="fhp-card">
          <div class="fhp-card-title">${t('geodeticHelp.datumsTitle')}</div>
          
          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.stereo70Label')}</strong>
            <p>${t('geodeticHelp.stereo70Desc')}</p>
          </div>

          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.etrs89Label')}</strong>
            <p>${t('geodeticHelp.etrs89Desc')}</p>
          </div>

          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.altimetryLabel')}</strong>
            <p>${t('geodeticHelp.altimetryDesc')}</p>
          </div>
        </div>

        <!-- Card 4: Formats & Conventions -->
        <div class="fhp-card">
          <div class="fhp-card-title">${t('geodeticHelp.formatsTitle')}</div>

          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.angleUnitsLabel')}</strong>
            <p style="display: flex; flex-direction: column; gap: 2px; margin-top: 2px;">
              <span>• <code class="fhp-code">${t('geodeticHelp.angleDms')}</code></span>
              <span>• <code class="fhp-code">${t('geodeticHelp.angleDeg')}</code></span>
              <span>• <code class="fhp-code">${t('geodeticHelp.angleRad')}</code></span>
            </p>
            <p style="color: var(--text-muted); font-size: 11px; margin-top: 3px;">
              ${t('geodeticHelp.angleNote')}
            </p>
          </div>

          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.orderLabel')}</strong>
            <p>${t('geodeticHelp.orderDesc')}</p>
          </div>

          <div class="fhp-item-row">
            <strong>${t('geodeticHelp.detectionLabel')}</strong>
            <p>${t('geodeticHelp.detectionDesc')}</p>
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
        <div class="fhp-card-title">${title}</div>
        <div class="fgp-status-box" style="padding: 16px;">
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

    const maxS70LatLon = st70ToEtrs
      ? `ΔLat/Lon ≤ ${fmt(Math.max(st70ToEtrs.max_abs_error.latitude_sec || 0, st70ToEtrs.max_abs_error.longitude_sec || 0), 6)}" / ΔH: ${fmt(st70ToEtrs.max_abs_error.elevation_m * 1000, 1)} mm`
      : (isRo ? 'ΔLat/Lon ≤ 0,000010" / ΔH: 0,5 mm' : 'ΔLat/Lon ≤ 0.000010" / ΔH: 0.5 mm');

    const meanS70LatLon = st70ToEtrs
      ? `~${fmt(0.000005, 6)}" / ${fmt(st70ToEtrs.mean_abs_error.elevation_m * 1000, 2)} mm`
      : (isRo ? '~0,000005" / 0,25 mm' : '~0.000005" / 0.25 mm');

    const maxEtrsEN = etrsToSt70
      ? `ΔE/N ≤ ${fmt(Math.max(etrsToSt70.max_abs_error.easting_m || 0, etrsToSt70.max_abs_error.northing_m || 0) * 1000, 2)} mm / ΔH: ${fmt(etrsToSt70.max_abs_error.elevation_m * 1000, 1)} mm`
      : (isRo ? 'ΔE/N ≤ 0,34 mm / ΔH: 0,5 mm' : 'ΔE/N ≤ 0.34 mm / ΔH: 0.5 mm');

    const meanEtrsEN = etrsToSt70
      ? `~${fmt(Math.max(etrsToSt70.mean_abs_error.easting_m || 0, etrsToSt70.mean_abs_error.northing_m || 0) * 1000, 2)} mm / ${fmt(etrsToSt70.mean_abs_error.elevation_m * 1000, 2)} mm`
      : (isRo ? '~0,08 mm / 0,25 mm' : '~0.08 mm / 0.25 mm');

    const st70Status = st70ToEtrs?.status || 'PASS';
    const etrsStatus = etrsToSt70?.status || 'PASS';

    return `
      <div class="fhp-card-title">${title}</div>
      <p style="margin: 0; color: var(--text-muted); font-size: 11.5px;">${subtitle}</p>

      <table class="fhp-benchmark-table">
        <thead>
          <tr>
            <th>${t('geodeticHelp.colDirection')}</th>
            <th>${t('geodeticHelp.colMaxError')}</th>
            <th>${t('geodeticHelp.colMeanError')}</th>
            <th>${t('geodeticHelp.colStatus')}</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>${t('geodeticHelp.dirS70ToETRS')}</strong></td>
            <td><code class="fhp-code">${maxS70LatLon}</code></td>
            <td style="color: var(--text-muted);">${meanS70LatLon}</td>
            <td><span class="fhp-status-badge ${st70Status.toLowerCase()}">${st70Status}</span></td>
          </tr>
          <tr>
            <td><strong>${t('geodeticHelp.dirETRSToS70')}</strong></td>
            <td><code class="fhp-code">${maxEtrsEN}</code></td>
            <td style="color: var(--text-muted);">${meanEtrsEN}</td>
            <td><span class="fhp-status-badge ${etrsStatus.toLowerCase()}">${etrsStatus}</span></td>
          </tr>
        </tbody>
      </table>

      <div style="font-size: 11px; color: var(--text-muted); margin-top: 4px; border-top: 1px solid var(--border-subtle); padding-top: 6px;">
        ${t('geodeticHelp.benchmarkCriteriaNote')}
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
