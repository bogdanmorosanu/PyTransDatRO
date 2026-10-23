# Task 36: Workbench Architecture, Header Redesign & Multi-Tab Harmonization

## Objectives
- [x] Streamline Global Header: GitHub UnderlineNav primary tabs, neutral gray indicator palette, status-only health indicator, text-only language selector, and 2-line attribution.
- [x] Unify Secondary Navigation Tabs: Eliminate redundant "TRANSFORM:" label prefix, standardize sleek secondary underline tabs across Point, Text, Files, and Usage sub-views.
- [x] Refine Global Configuration Card: Elevated `.shared-config-bar` matching `.studio-card`, discrete `SETTINGS` micro-label, 3-row layout with 28px sliding gliders (`Operation Direction`, `Coordinate Order`), compact vertically-centered Angular Unit dropdown, and format Example info chip.
- [x] Establish Three-Card Workbench Architecture: Symmetrically applied across `PointConverter.ts`, `BatchConverter.ts`, and `FileUpload.ts`:
  - **Card 1: `INPUT`**: Source input controls, text-only actions.
  - **Active Separator**: Full-width primary `Transform` action button (34px `.btn-primary`).
  - **Card 2: `OUTPUT`**: Read-only outputs, code-block copy trigger, downloads (`CSV`, `GeoJSON`, `DXF`, file downloads, `Download All`, `Clear List`).
  - **Card 3: `REPORT`**: Two-row dataset quality summary (`Coordinates` $\to$ neutral pills `Total (N)`, `2D (N)`, `3D (N)`; `Transformation` $\to$ semantic status pills without icons `Valid (N)`, `Outside Grid (N)`, `Invalid Input (N)`).
  - **Bottom Dock**: Decoupled `.workbench-actions-dock` containing 26px `.btn-secondary-action` buttons (`Center on Map`, `Zoom to Extents`).
- [x] Refine Point Workbench: Fix DMS seconds symbol (`"`) escaping via `escapeHtmlAttr()`, non-editable Quasigeoid height info chip (`.info-field-badge`), empty outputs when uncalculated.
- [x] Refine Text Workbench: Scope `height: 30px` to `input.form-input` so textareas preserve 6+ visible lines (`min-height: 130px; line-height: 1.45`), standardize 26px secondary download buttons.
- [x] Refine Files Workbench: Integrate "Ignore header / first line" toggle switch (Archetype 2.4), add `Clear List` button on left of `Download All`, remove icons from `Clear Queue`, rename button to `Transform`, strip per-file coordinate count clutter while preserving `✕` remove button.
- [x] Color & Contrast Tokens: Calibrate `--bg-surface-readonly` in dark mode to `#181d24` (distinct from `#0b0d10` input black) and ensure neutral pills and readonly surfaces are unmistakably non-editable in both themes.
- [x] Living Documentation: Synchronize `docs/webapp/archetypes.md`, `docs/webapp/design_system.md`, and `docs/webapp/actionable_inventory.md`.

---

## Implementation Details

### 1. Global Header Chrome & Identity (`Header.ts`, `badges.css`, `layout.css`)
- **Primary UnderlineNav Tabs**: Translucent slate hover background (`--tab-hover-bg`) with a 10px vertical gap above the header line; neutral gray indicator underline (`--tab-indicator-color`: `#d0d7de` dark, `#1f2328` light).
- **Status-Only Health Indicator**: Replaced clickable discrepancy modal button with `<span class="health-dot" id="header-health-dot">` providing native hover tooltips without modal triggers.
- **Language Selector**: Removed decorative SVG globe icon, providing a clean `RO ▾` / `EN ▾` disclosure trigger.
- **Attribution**: 2-line attribution block (`powered by` / `PyTransDatRO`) pointing to GitHub without oblique external link arrows.
- **Dynamic Title**: Synchronized static `<title>` in `index.html` and dynamic `document.title` with `APP_CONFIG.appName`.

### 2. Global Configuration Card (`PointConverter.ts`, `forms.css`, `layout.css`)
- Framed `.shared-config-bar` as an elevated card matching `.studio-card` (`border-radius: var(--radius-lg); background: var(--bg-surface-elevated); border: 1px solid var(--border-medium);`).
- Added discrete `.section-micro-label`: `SETTINGS` / `SETĂRI`.
- Compressed into 3 horizontal rows (`grid-template-columns: 125px 1fr`):
  - Row 1: `Operation Direction` $\to$ 28px sliding glider.
  - Row 2: `Coordinate Order` $\to$ 28px sliding glider.
  - Row 3: `Angular Unit` $\to$ 28px vertically centered select (`.config-unit-select`) + `Example` info badge.

### 3. Three-Card Workbench Architecture across Converters
- **`PointConverter.ts`**:
  - Organized into Card 1 (`INPUT`), Active Separator (`Transform`), Card 2 (`OUTPUT`), Card 3 (`REPORT`), and `.workbench-actions-dock`.
  - Added `escapeHtmlAttr()` to escape `"` in DMS formatted strings (`46°00'00.00000"`).
  - Quasigeoid height converted to non-editable info badge (`.info-field-badge`).
  - Cleared / uncalculated outputs render empty string `""` instead of persistent dashes.
- **`BatchConverter.ts`**:
  - Organized into Card 1 (`INPUT`), Active Separator (`Transform`), Card 2 (`OUTPUT`), Card 3 (`REPORT`), and `.workbench-actions-dock`.
  - Fixed textarea height: scoped `height: 30px` to `input.form-input` and added `textarea.form-input { min-height: 130px; line-height: 1.45; }`.
  - Standardized `.btn-secondary-action` (26px height) across `Download CSV`, `GeoJSON`, `DXF`, `Center on Map`, and `Zoom to Extents`.
  - Removed decorative icons (`✓`, `⚠`, `✕`) from status pills; unified parentheses count notation `Total (4)`, `Valid (3)`.
- **`FileUpload.ts`**:
  - Organized into Card 1 (`INPUT`), Active Separator (`Transform`), Card 2 (`OUTPUT`), Card 3 (`REPORT`), and `.workbench-actions-dock`.
  - Added "Ignore header / first line" toggle switch (`#file-skip-header-toggle`, Archetype 2.4).
  - Renamed primary action button to `Transform` across all locales.
  - Removed icon from `Clear Queue` button in Card 1.
  - Added `Clear List` button (`#file-clear-output-btn`, 26px) adjacent to `Download All` in Card 2.
  - Stripped per-file coordinate count clutter from output cards, showing only file name, status pill, download button, and `✕` remove button.

---

## Verification & Testing

- **TypeScript & Bundle Build**:
  - `npm --prefix webapp run build` (`tsc && vite build`) passed in 399ms (0 errors).
- **Vitest Unit Suite**:
  - `npm --prefix webapp test` passed all 80 tests across 6 test suites (`pointParser`, `dmsParser`, `exporters`, `gridBounds`, `ui`, `numberFormat`).
- **Playwright End-to-End Suite**:
  - `npm --prefix webapp run test:e2e -- e2e/file.spec.ts` passed in 6.7s.
- **Living Documentation**:
  - Synchronized `docs/webapp/archetypes.md` (Archetype 2.4).
  - Synchronized `docs/webapp/design_system.md` (Three-Card Workbench Architecture).
  - Synchronized `docs/webapp/actionable_inventory.md` (Table controls).
  - Logged session outcomes in `walkthrough.md`.
