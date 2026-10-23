# Task 39: Usage Telemetry Tabs Harmonization, Modal Refinements & Brand Geometry

## Objectives
- [x] **Usage Navigation Header**: Remove subtitle and vertical separator so tabs are flush left on the Computation, Spatial, and Legacy views.
- [x] **Archetype 8.1 KPI Cards**: Standardize all KPI cards with flex-centered alignment and clean sentence-case labels across Computation (3 cards), Spatial (2 cards), and Legacy (5 cards).
- [x] **Archetype 4.1 Date Range Limits**: Replace wide pill banners with paired Info field archetypes (`.info-label` + `.info-field-readonly`) displaying full calendar day resolution (`YYYY-MM-DD`).
- [x] **Archetype 2.3 Enclosed Multi-Select Trays**: Group multi-select options (Sources, Operations, Metrics) in unified enclosed pill containers (`.telemetry-btn-toggle-group.multi-select`) with outer borders, borderless buttons, subtle hover fill, and illuminated active dot indicators (`::before`).
- [x] **Archetype 6.6 Stationary Total Toggle**: Reposition the Total line toggle to top-left (`.legend-total-cluster`) before swatches, stripping redundant icon `✓` and label prefix, ensuring 0px layout shift during toggle interactions.
- [x] **Spatial Tab Ergonomics**: Split the merged points/cells KPI into two distinct cards, retain mutually-exclusive segmented tray styling for Points/Calls metric selector, and apply clean sentence-case titles with cardinal abbreviations for quadrants.
- [x] **Volume Warning Modal Polish**: Remove all hardcoded emojis and icons (`⚠`, `⚡`, `✕`) from title and actions; restyle the close button to a compact 24×24px square icon button matching `.fgp-icon-btn`.
- [x] **Header Brand Geometry**: Increase vertical gap between `PyTrans.ro` title and version badge by ~3.5px while preserving header vertical centering, left-aligned version badge, and right-aligned online health dot.
- [x] **E2E & Living Documentation**: Author `computation.spec.ts` and `spatial.spec.ts` E2E suites; synchronize `actionable_inventory.md`, `archetypes.md`, `design_system.md`, `architecture.md`, `testing.md`, and `webapp/README.md`.

---

## Implementation Details

### 1. Telemetry Dashboard & Chart Visual Harmonization (`webapp/src/components/`, `webapp/src/style/telemetry.css`)
- **Navigation Shell**: Removed `.telemetry-sub-title` and `.telemetry-tab-sep` across all telemetry views, aligning navigation tabs flush left with the studio layout.
- **KPI Cards (Archetype 8.1)**:
  - Centered value readouts and labels using `display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center;`.
  - Replaced uppercase typography with clean sentence case (`Total calls`, `Total points`, `REST API`, `Web app`, `Total volume`, `Active cells`, `Stereo 70 points`, `Stereo 30 points`).
  - Added `.telemetry-kpi-banner.two-cards` and `.telemetry-kpi-banner.five-cards` grid templates with responsive column collapse.
- **Date Range Limits (Archetype 4.1)**:
  - Replaced `.telemetry-range-limits-pill` with `.telemetry-date-limits-cluster`.
  - Structured each endpoint with `<label class="info-label">` and `<span class="info-field-readonly">`.
  - Guaranteed full calendar day resolution (`YYYY-MM-DD`) by falling back to `-01` when ISO day strings are truncated.
  - Added localized strings for `archiveStart` and `archiveEnd` in `en.ts` and `ro.ts`.
- **Enclosed Multi-Select Trays (Archetype 2.3)**:
  - Grouped multi-select controls into `.telemetry-btn-toggle-group.multi-select`.
  - Applied outer bounding container border (`1px solid var(--border-medium); border-radius: var(--radius-sm); overflow: hidden;`).
  - Stripped individual button borders (`border: none; border-radius: 0; background: transparent;`).
  - Implemented dot indicator on `.telemetry-toggle-btn.btn-active::before` with `--color-accent-teal` illumination.
  - Kept radio segmented pill tray styling for mutually-exclusive metrics in Spatial view.
- **Stationary Total Toggle (Archetype 6.6)**:
  - Restructured timeline legend to anchor `.legend-total-cluster` at top-left.
  - Removed icon prefix and redundant "Total:" label.
  - Dynamic legend swatches expand to the right of the button, eliminating cursor layout shift.
- **SVG Chart Axis Typography**:
  - Reduced vertical axis titles (`Calls`, `Points`) to 9.5px, 500 font-weight, sentence-case.

### 2. Volume Warning Modal & Header Brand Geometry (`VolumeWarningModal.ts`, `layout.css`)
- **Volume Warning Modal**:
  - Removed icons from title and buttons (`Render sample`, `Render all`, `Skip map rendering`).
  - Restyled `#vw-close-icon` with `.modal-close-btn` (`width: 24px; height: 24px; padding: 0; border: none; font-size: 12px; border-radius: var(--radius-sm)`).
- **Header Brand Typography**:
  - Adjusted `.brand-text-column` to `gap: 3.5px; justify-content: center;`.
  - Set `.brand-title-text` to `text-align: center; width: 100%;`.
  - Maintained `.brand-meta-subrow` with `width: 100%; justify-content: space-between;` ensuring square proportions in the 52px header bar.

### 3. Living Documentation & Tests (`docs/webapp/`, `webapp/e2e/`)
- Created Playwright browser tests:
  - `webapp/e2e/computation.spec.ts`: Validates centered KPI cards, Archetype 4.1 date limits, multi-select trays, and stationary total toggle.
  - `webapp/e2e/spatial.spec.ts`: Validates 2 KPI cards, quadrant matrix, and segmented radio metric trays.
- Updated `docs/webapp/actionable_inventory.md`, `docs/webapp/archetypes.md`, `docs/webapp/design_system.md`, `docs/webapp/architecture.md`, `docs/testing.md`, and `webapp/README.md`.

---

## Verification & Quality Gates
- **Vitest Unit Suite**: 80/80 passed across all 6 test files.
- **TypeScript & Vite Build**: Passed cleanly with zero compilation errors (252ms).
- **Core Geodetic Suite** (`tests/test_trans_ro.py`): 54/54 passed with zero numerical regression ($< 10^{-9}$ tolerance).
- **API & Configuration Suites**: 46/46 passed.
- **Python Code Style & Hygiene**: Ruff check passed cleanly with zero linting or formatting defects.
