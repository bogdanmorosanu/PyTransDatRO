import { t } from '../i18n';
import { apiClient, type GridMetadataResponse } from '../utils/apiClient';
import { showToast } from '../utils/toast';

export class GridMetadataPanel {
  private element: HTMLElement | null = null;
  private data: GridMetadataResponse | null = null;
  private isLoading: boolean = false;
  private searchQuery: string = '';
  private expandedPaths: Set<string> = new Set();
  private isMinimized: boolean = false;

  // Dragging state
  private isDragging: boolean = false;
  private dragStartX: number = 0;
  private dragStartY: number = 0;
  private initialLeft: number = 0;
  private initialTop: number = 0;

  constructor() {
    // Initial auto-expanded top-level roots
    this.expandedPaths.add('file_info');
    this.expandedPaths.add('metadata');
    this.expandedPaths.add('params');
    this.expandedPaths.add('grids');
    this.expandedPaths.add('grids.geodetic_shifts');
    this.expandedPaths.add('grids.geoid_heights');
  }

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
    this.element.className = 'floating-grid-panel';
    document.body.appendChild(this.element);

    this.render();
    await this.loadData();
  }

  public close(): void {
    if (this.element) {
      this.element.remove();
      this.element = null;
    }
  }

  private async loadData(): Promise<void> {
    this.isLoading = true;
    this.renderBody();

    try {
      this.data = await apiClient.getGridMetadata();
    } catch (err) {
      console.error('Failed to load grid metadata:', err);
    } finally {
      this.isLoading = false;
      this.render();
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

  private expandAll(): void {
    if (!this.data) return;
    const collectPaths = (obj: any, currentPath: string = '') => {
      if (typeof obj === 'object' && obj !== null) {
        if (currentPath) this.expandedPaths.add(currentPath);
        for (const k of Object.keys(obj)) {
          const next = currentPath ? `${currentPath}.${k}` : k;
          collectPaths(obj[k], next);
        }
      }
    };
    collectPaths(this.data);
    this.renderBody();
  }

  private collapseAll(): void {
    this.expandedPaths.clear();
    this.renderBody();
  }

  private copyRawJson(): void {
    if (!this.data) return;
    const jsonStr = JSON.stringify(this.data, null, 2);
    navigator.clipboard.writeText(jsonStr).then(() => {
      showToast(t('gridMetadata.copiedJson') || 'Copied JSON to clipboard');
    });
  }

  private render(): void {
    if (!this.element) return;

    const sizeStr = this.data?.file_info?.size_formatted || '';
    const shaStr = this.data?.file_info?.sha256 || '';

    this.element.innerHTML = `
      <!-- Draggable Header -->
      <div class="drawer-header fgp-header" id="fgp-drag-handle">
        <div class="drawer-title-cluster fgp-title-cluster">
          <span class="drawer-title fgp-title">${t('gridMetadata.title')}</span>
        </div>
        <div class="drawer-header-actions fgp-header-actions">
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fgp-minimize-btn" title="${this.isMinimized ? t('common.restore') : t('common.minimize')}">
            ${this.isMinimized ? '□' : '—'}
          </button>
          <button class="btn-icon drawer-icon-btn fgp-icon-btn" id="fgp-close-btn" title="${t('common.close')}">✕</button>
        </div>
      </div>

      <!-- Main Body Container -->
      <div class="fgp-content" style="${this.isMinimized ? 'display: none;' : ''}">
        <!-- Authoritative Origin Notice -->
        <div class="callout-box fhp-callout-box fgp-callout-box">
          <span class="callout-title fhp-callout-title">${t('gridMetadata.originTitle')}</span>
          <p class="callout-body fhp-callout-body">${t('gridMetadata.sourceNotice')}</p>
          ${
            shaStr
              ? `
            <div class="fgp-meta-chips-row">
              <div class="fgp-meta-chip-item">
                <span class="fgp-meta-label">${t('gridMetadata.sizeLabel')}</span>
                <code class="fhp-metric-chip mono">${sizeStr}</code>
              </div>
              <div class="fgp-meta-chip-item">
                <span class="fgp-meta-label">${t('gridMetadata.sha256Label')}</span>
                <code class="fhp-metric-chip mono" title="${shaStr}">${shaStr.slice(0, 16)}...</code>
              </div>
            </div>
          `
              : ''
          }
        </div>

        <!-- Controls Toolbar -->
        <div class="fgp-toolbar">
          <div class="fgp-search-box">
            <svg width="12" height="12" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5">
              <circle cx="7" cy="7" r="5"/>
              <path d="M11 11l3.5 3.5"/>
            </svg>
            <input 
              type="text" 
              class="fgp-search-input" 
              id="fgp-search" 
              placeholder="${t('gridMetadata.searchPlaceholder') || 'Filter keys or values...'}"
              value="${this.searchQuery}"
            />
            ${this.searchQuery ? `<button class="btn-icon fgp-clear-search" id="fgp-clear-search">✕</button>` : ''}
          </div>
          <div class="fgp-btn-row">
            <button class="btn btn-outline fgp-tool-btn" id="fgp-expand-all-btn">${t('gridMetadata.expandAll')}</button>
            <button class="btn btn-outline fgp-tool-btn" id="fgp-collapse-all-btn">${t('gridMetadata.collapseAll')}</button>
            <button class="btn btn-outline fgp-tool-btn" id="fgp-copy-json-btn">${t('gridMetadata.copyJson')}</button>
          </div>
        </div>

        <!-- Tree View Container -->
        <div class="fgp-tree-container" id="fgp-tree-slot">
          <!-- Dynamic Tree Rendered Here -->
        </div>

        <!-- Corner Resize Grip -->
        <div class="fgp-resize-handle" id="fgp-resize-handle" title="Resize window">
          <svg width="10" height="10" viewBox="0 0 10 10" fill="currentColor">
            <path d="M9 1v8H1l8-8zM9 5v4H5l4-4z"/>
          </svg>
        </div>
      </div>
    `;

    this.attachHeaderListeners();
    this.attachDragListeners();
    this.attachResizeListeners();
    this.renderBody();
  }

  private attachHeaderListeners(): void {
    if (!this.element) return;

    this.element.querySelector('#fgp-close-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.close();
    });

    this.element.querySelector('#fgp-minimize-btn')?.addEventListener('click', (e) => {
      e.stopPropagation();
      this.toggleMinimize();
    });

    this.element.querySelector('#fgp-expand-all-btn')?.addEventListener('click', () => this.expandAll());
    this.element.querySelector('#fgp-collapse-all-btn')?.addEventListener('click', () => this.collapseAll());
    this.element.querySelector('#fgp-copy-json-btn')?.addEventListener('click', () => this.copyRawJson());

    const searchInput = this.element.querySelector('#fgp-search') as HTMLInputElement | null;
    searchInput?.addEventListener('input', (e) => {
      this.searchQuery = (e.target as HTMLInputElement).value;
      this.renderBody();
    });

    this.element.querySelector('#fgp-clear-search')?.addEventListener('click', () => {
      this.searchQuery = '';
      if (searchInput) searchInput.value = '';
      this.renderBody();
    });
  }

  private renderBody(): void {
    const slot = this.element?.querySelector('#fgp-tree-slot');
    if (!slot) return;

    if (this.isLoading) {
      slot.innerHTML = `
        <div class="fgp-status-box">
          <div class="fgp-spinner"></div>
          <span>${t('gridMetadata.loading') || 'Loading grid metadata...'}</span>
        </div>
      `;
      return;
    }

    if (!this.data) {
      slot.innerHTML = `
        <div class="fgp-status-box error">
          <span>${t('gridMetadata.errorLoading') || 'Failed to load grid metadata.'}</span>
        </div>
      `;
      return;
    }

    const treeHtml = this.buildTreeHtml(this.data, '');
    if (!treeHtml.trim()) {
      slot.innerHTML = `
        <div class="fgp-status-box">
          <span>${t('gridMetadata.noMatches') || 'No metadata items match your filter.'}</span>
        </div>
      `;
      return;
    }

    slot.innerHTML = `<div class="fgp-tree">${treeHtml}</div>`;
    this.attachTreeToggleListeners(slot);
  }

  private buildTreeHtml(obj: any, currentPath: string, depth: number = 0): string {
    if (typeof obj !== 'object' || obj === null) {
      return '';
    }

    let html = '';
    const query = this.searchQuery.trim().toLowerCase();
    const entries = Object.entries(obj);

    for (const [key, val] of entries) {
      const fullPath = currentPath ? `${currentPath}.${key}` : key;
      const isObject = typeof val === 'object' && val !== null;
      const isArray = Array.isArray(val);

      // Filtering logic
      if (query && !this.pathOrChildMatches(key, val, query)) {
        continue;
      }

      const isExpanded = this.expandedPaths.has(fullPath) || Boolean(query);
      const indentPx = depth * 18;

      if (isObject) {
        const count = isArray ? val.length : Object.keys(val).length;
        const countTag = isArray ? `[${count}]` : `{${count}}`;

        html += `
          <div class="fgp-tree-branch" data-path="${fullPath}">
            <div class="fgp-tree-row branch-row" style="padding-left: ${indentPx}px;">
              <span class="fgp-expander ${isExpanded ? 'open' : ''}">▸</span>
              <span class="fgp-key branch-key">${key}</span>
              <span class="fgp-type-badge">${countTag}</span>
            </div>
            <div class="fgp-tree-children ${isExpanded ? 'open' : ''}" style="${isExpanded ? '' : 'display: none;'}">
              ${this.buildTreeHtml(val, fullPath, depth + 1)}
            </div>
          </div>
        `;
      } else {
        const valStr = String(val);
        const valType = typeof val;
        let formattedVal = valStr;

        if (valType === 'string') {
          // If HTML (like attribution), escape or render cleanly
          if (valStr.includes('<') && valStr.includes('>')) {
            formattedVal = `<span class="fgp-val-html" title="${valStr.replace(/"/g, '&quot;')}">${valStr}</span>`;
          } else {
            formattedVal = `<span class="fgp-val-string">"${valStr}"</span>`;
          }
        } else if (valType === 'number') {
          formattedVal = `<span class="fgp-val-number mono">${valStr}</span>`;
        } else if (valType === 'boolean') {
          formattedVal = `<span class="fgp-val-bool ${val ? 'true' : 'false'}">${valStr}</span>`;
        } else if (val === null) {
          formattedVal = `<span class="fgp-val-null">null</span>`;
        }

        html += `
          <div class="fgp-tree-row leaf-row" style="padding-left: ${indentPx + 14}px;">
            <span class="fgp-key leaf-key">${key}:</span>
            <span class="fgp-val">${formattedVal}</span>
            <button class="btn-icon fgp-copy-val-btn" data-copy="${valStr.replace(/"/g, '&quot;')}" title="Copy value">
              <svg width="10" height="10" viewBox="0 0 16 16" fill="currentColor">
                <path d="M4 2a2 2 0 00-2 2v8a2 2 0 002 2h6a2 2 0 002-2V4a2 2 0 00-2-2H4zm0 1h6a1 1 0 011 1v8a1 1 0 01-1 1H4a1 1 0 01-1-1V4a1 1 0 011-1z"/>
              </svg>
            </button>
          </div>
        `;
      }
    }

    return html;
  }

  private pathOrChildMatches(key: string, val: any, query: string): boolean {
    if (key.toLowerCase().includes(query)) return true;
    if (val === null || val === undefined) return false;
    if (typeof val !== 'object') {
      return String(val).toLowerCase().includes(query);
    }
    for (const [k, v] of Object.entries(val)) {
      if (this.pathOrChildMatches(k, v, query)) return true;
    }
    return false;
  }

  private attachTreeToggleListeners(container: Element): void {
    container.querySelectorAll('.branch-row').forEach((row) => {
      row.addEventListener('click', (e) => {
        e.stopPropagation();
        const branch = row.closest('.fgp-tree-branch') as HTMLElement | null;
        if (!branch) return;
        const path = branch.getAttribute('data-path');
        if (!path) return;

        const children = branch.querySelector(':scope > .fgp-tree-children') as HTMLElement | null;
        const expander = branch.querySelector(':scope > .branch-row .fgp-expander');

        if (this.expandedPaths.has(path)) {
          this.expandedPaths.delete(path);
          if (children) children.style.display = 'none';
          expander?.classList.remove('open');
        } else {
          this.expandedPaths.add(path);
          if (children) children.style.display = 'block';
          expander?.classList.add('open');
        }
      });
    });

    container.querySelectorAll('.fgp-copy-val-btn').forEach((btn) => {
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const text = btn.getAttribute('data-copy');
        if (text) {
          navigator.clipboard.writeText(text).then(() => {
            showToast('✓ Copied');
          });
        }
      });
    });
  }

  private attachDragListeners(): void {
    if (!this.element) return;
    const handle = this.element.querySelector('#fgp-drag-handle') as HTMLElement | null;
    if (!handle) return;

    // Prevent any click/mousedown inside the floating window from bubbling up to Leaflet or other handlers
    this.element.addEventListener('mousedown', (e) => {
      e.stopPropagation();
    });
    this.element.addEventListener('pointerdown', (e) => {
      e.stopPropagation();
    });

    handle.addEventListener('pointerdown', (e: PointerEvent) => {
      if ((e.target as HTMLElement).closest('.fgp-header-actions')) return;
      if (window.innerWidth <= 768) return;
      e.preventDefault();
      e.stopPropagation();

      this.isDragging = true;
      this.dragStartX = e.clientX;
      this.dragStartY = e.clientY;

      const rect = this.element!.getBoundingClientRect();
      this.initialLeft = rect.left;
      this.initialTop = rect.top;

      // Capture pointer so dragging is smooth and cannot be intercepted by Leaflet map tiles
      handle.setPointerCapture(e.pointerId);

      const onPointerMove = (moveEvent: PointerEvent) => {
        if (!this.isDragging || !this.element) return;
        moveEvent.preventDefault();
        moveEvent.stopPropagation();

        const dx = moveEvent.clientX - this.dragStartX;
        const dy = moveEvent.clientY - this.dragStartY;

        // Keep inside viewport bounds
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
    const resizer = this.element.querySelector('#fgp-resize-handle') as HTMLElement | null;
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

        const minWidth = 340;
        const maxWidth = Math.max(minWidth, window.innerWidth - 40);
        const minHeight = 220;
        const maxHeight = Math.max(minHeight, window.innerHeight - 80);

        const newWidth = Math.max(minWidth, Math.min(startWidth + dw, maxWidth));
        const newHeight = Math.max(minHeight, Math.min(startHeight + dh, maxHeight));

        this.element.style.width = `${newWidth}px`;
        this.element.style.height = `${newHeight}px`;
        this.element.style.maxHeight = 'none'; // allow manual sizing beyond default max-height
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
