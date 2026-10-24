# PyTransDatRO Web Application — Design System & Visual Manual

This document is the authoritative visual design manual for the **PyTransDatRO** web application (`webapp/`). It defines the core design philosophy, the **GitHub Primer neutral Dark & Light theme**, the **Subtle Surface + Border (`#8`) primary action button**, typography, structural layout invariants, header navigation rules, Web GIS map widget styling, and the **Map Projection Graticule dual-theme brand logo**.

> [!NOTE]
> Component grouping and multi-instance styling classifications live in [`archetypes.md`](./archetypes.md). Purely algorithmic and data-processing specifications live in [`architecture.md`](./architecture.md).

---

## 1. Core Design Philosophy: The Trinity & Singularity of Purpose

### 1.1 Singularity of Purpose
PyTransDatRO has a dedicated, singular mission: **accurate, reliable coordinate transformations** between Romania's national projection (**Stereo 70 / Black Sea 1975**) and European/global GIS standards (**ETRS89 / GRS80**).
- Avoid decorative clutter or visual distractions that detract from this core spatial workflow.
- Prioritize modern GIS and cartographic workflows: speed, legibility, uncluttered hierarchy, and frictionless copy-paste interoperability with GIS and CAD environments.

### 1.2 The Trinity Balance

```
                 [ Simplicity ]
                 /            \
                / PyTransDatRO \
               /  Design Center \
        [ Optimal ] ---------- [ Modern ]
```

1. **Simplicity (Clarity over Decoration)**:
   - **Zero superfluous emojis or decorative icons**: Tab labels, primary action buttons, dialog choices, and segmented options never contain informal emojis (no ⚡, ⬍, ⬄, 🧭). Clean typography and precise vector geometry communicate authority.
   - **Subtle micro-headers**: Label functional groups with quiet, structured micro-labels (`Input coordinates (or click on map)`, `Computed coordinates`, `Download As:`).
   - **Visual quietness**: Secondary configuration toggles remain neutral while exhibiting clear ON/OFF state contrast, and the primary transformation button integrates naturally across the `480px` panel without harsh color glare.

2. **Optimal (Spatial Ergonomics & Layout Invariants)**:
   - **Fixed `480px` Left Workbench Panel (`.utility-panel`)**: Specifically proportioned to house multi-column coordinate tuples comfortably without encroaching on the interactive Web GIS map viewport.
   - **Asymmetrical 3-Column Coordinate Grid (`1.35fr 1.35fr 0.8fr`)**: Northing/Latitude (`minmax(0, 1.35fr)`), Easting/Longitude (`minmax(0, 1.35fr)`), and Altitude/Height (`minmax(0, 0.8fr)`) are sized directly to the `480px` panel so long unspaced DMS strings (`46°00'00.00000"`) never wrap or clip.
   - **Top-right code-block copy buttons**: Icon-based copy triggers with instant inline feedback (`✓ Copied`) provide zero-friction extraction.

3. **Modern (Precision GIS & Cartographic Engineering)**:
   - **GitHub Primer Neutral Dark & Daylight Surfaces (Zero Dark Navy-Blue)**: Dark mode uses an authentic carbon-slate `#0d1117` canvas and `#161b22` card elevation matching the dark Web GIS map widgets—never saturated dark navy or indigo (`#090d16`, `#0f172a`, `#1e293b`). Light mode uses `#f6f8fa` canvas and `#ffffff` surfaces.
   - **Subtle Surface + Border Primary Action Button**: Because the `Transform` button spans the full width of the `480px` left workbench panel, it uses an elevated neutral surface fill with a crisp structural border (`--btn-primary-bg: #1f2630`, `--btn-primary-border: #48586a`, `--btn-primary-text: #e6edf3` in dark mode; `--btn-primary-bg: #eaeef2`, `--btn-primary-border: #8c959f`, `--btn-primary-text: #1f2328` in light mode) rather than a loud saturated blue or green block.
   - **Strict Semantic & Metric Color Isolation**: Physical spatial measurements—most notably **Quasigeoid Height** ($\zeta = h - H$)—use the dedicated cyan/teal metric palette (`--color-metric-text`), strictly isolated from system status colors (green, amber, red).

---

## 2. GitHub Primer Dual-Theme Palette & Primary Button Tokens

All visual styling is governed by semantic CSS custom properties in `webapp/src/style/variables.css` across Dark (`data-theme="dark"`) and Light (`data-theme="light"`) modes.

### 2.1 The Canonical Thematic Color Palette Matrix

All visual styling is governed by semantic CSS custom properties in `webapp/src/style/variables.css` across Dark (`data-theme="dark"`) and Light (`data-theme="light"`) modes. To uphold our **simplicity directive** without compromising contrast or usability, all visual elements draw strictly from this defined palette:

| Token Category | Semantic Token | Dark Theme (`data-theme="dark"`) | Light Theme (`data-theme="light"`) | Usability / Architectural Role |
| :--- | :--- | :--- | :--- | :--- |
| **Surfaces (Tiers 0–4)** | `--bg-app` | `#0e1013` (Carbon) | `#f6f8fa` (Daylight) | Application root canvas background. |
| | `--bg-surface` | `#14171c` (Graphite) | `#ffffff` (Pure White) | Primary card surfaces, header bar, workbench panel. |
| | `--bg-surface-elevated`| `#1c2026` (Slate) | `#f3f4f6` (Light Gray) | Elevated content cards, drawer cards, secondary buttons. |
| | `--bg-surface-input` | `#0b0d10` (Dark Inset) | `#ffffff` (Inset White) | Inset coordinate inputs, batch textareas, notice boxes. |
| | `--bg-surface-readonly`| `#181d24` (Subtle Inset)| `#eaeef2` (Muted Daylight)| Non-editable output boxes, copy badges. |
| | `--bg-surface-hover` | `#252a32` | `#eaeef2` | Interactive button & menu item hover state. |
| **Borders (Tiers 1–4)** | `--border-subtle` | `#22272e` | `#d8dee4` | Card outlines, solid section dividers, inactive boundaries. |
| | `--border-medium` | `#30363d` | `#d0d7de` | Input boundaries, elevated control outlines, drawer borders. |
| | `--border-focus` | `#7c93ab` | `#57606a` | Keyboard navigation `:focus-visible` ring. |
| | `--border-dashed` | `#30363d` | `#afb8c1` | Topic separation dashed lines, total curve indicators. |
| **Typography (Tiers 1–3)**| `--text-main` | `#e6edf3` | `#1f2328` | Primary readable text, coordinate values, active tab text. |
| | `--text-muted` | `#8b949e` | `#57606a` | Standard field labels (`.form-label`), inactive options. |
| | `--text-dim` | `#6e7681` | `#6e7781` | Subordinate labels (`.info-label`), micro-headers. |
| **Primary CTA Button** | `--btn-primary-bg` | `#1f2630` | `#eaeef2` | Option #8 Subtle Surface primary action fill. |
| | `--btn-primary-border`| `#48586a` | `#8c959f` | Option #8 Subtle Surface structural border. |
| | `--btn-primary-text` | `#e6edf3` | `#1f2328` | Option #8 High-contrast readable button label. |
| | `--btn-primary-hover`| `#28313e` | `#dde3ea` | Option #8 Interactive hover surface. |
| **Isolated Metric** | `--color-metric-text`| `#56b6c2` (Cyan/Teal) | `#0e7490` (Deep Teal) | **Quasigeoid Height ($\zeta$)** physical spatial measurement. |
| | `--color-metric-border`| `rgba(86, 182, 194, 0.32)` | `rgba(14, 116, 144, 0.28)` | Quasigeoid height readout border. |
| | `--color-metric-bg` | `rgba(86, 182, 194, 0.08)` | `rgba(14, 116, 144, 0.06)` | Quasigeoid height readout background. |
| **System Status** | `--color-success` | `#10b981` (Green) | `#10b981` (Green) | In-grid validity, API online health. |
| | `--color-warning` | `#f59e0b` (Amber) | `#f59e0b` (Amber) | Out-of-grid notice, large point volume warning. |
| | `--color-danger` | `#ef4444` (Red) | `#ef4444` (Red) | Input syntax errors, network communication failures. |
| **Spatial Heatmap** | `--heatmap-step-1..4`| `#1e3a5f` .. `#a5f3fc` | `#0969da` .. `#80ccff` | Spatial choropleth matrix heat scale (isolated to map/chart). |

### 2.2 Primary Action Button Tokens (`Subtle Surface + Border`)

```css
/* Dark Theme (:root, :root[data-theme="dark"]) */
--btn-primary-bg: #1f2630;
--btn-primary-hover: #28313e;
--btn-primary-text: #e6edf3;
--btn-primary-border: #48586a;

/* Light Theme (:root[data-theme="light"]) */
--btn-primary-bg: #eaeef2;
--btn-primary-hover: #dde3ea;
--btn-primary-text: #1f2328;
--btn-primary-border: #8c959f;
```

### 2.3 Semantic & Metric Colors (Strict Isolation)
- **Primary Technical Metric (`--color-metric-text`)**:
  - Dark Theme: `#56d4dd` | Light Theme: `#0550ae`
  - **Purpose**: Strictly reserved for spatial and cartographic scalar measurements, notably the **Quasigeoid Height** ($\zeta = h - H$).
- **Success (`#10b981`)**: Strictly for verified API health and valid in-grid coordinates (`In Grid` / `Valid`).
- **Warning (`#f59e0b`)**: Strictly for cautionary notices (`Outside Grid`, `No Grid Data`, or high-volume point warnings > 2,000 points).
- **Danger / Error (`#ef4444`)**: Strictly for input syntax failures (`Invalid Values`, `Invalid Angle`), network failures, or unparseable input.

> [!IMPORTANT]
> **Data Values Rule**: Quasigeoid height $\zeta$, coordinate shifts $\Delta E, \Delta N$, and coordinate tuples must **never** be rendered in green, red, or yellow. They are physical scalar measurements and must use either neutral monospace typography (`--text-main`) or the dedicated technical metric accent palette (`--color-metric-text`).

> [!NOTE]
> **Map Marker Symbology & Layer Color Preservation**: Single-point map markers do not turn red when clicked or positioned outside the Romanian grid boundary. The marker strictly honors its configured layer color (e.g., `#3b82f6` for Stereo 70, `#10b981` for ETRS89) and features continuous soft pulsing animation (`.marker-pulsing`). Grid validity is communicated exclusively via the Two-Section Quality Report card and coordinate readouts.

---

## 3. Typography, Geometry & Formatting Conventions

### 3.1 Typography
- **UI Sans Font (`--font-sans`)**: `'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif`
- **Coordinate & Metric Monospace Font (`--font-mono`)**: `'JetBrains Mono', 'Fira Code', Menlo, Monaco, Consolas, monospace`
  - Mandatory for all coordinate inputs, read-only coordinate outputs, Quasigeoid height values, angle format badges, and map coordinate readouts.

### 3.2 Structural Layout Invariants (`480px` Workbench & 3-Column Coordinate Grid)
1. **Fixed `480px` Left Workbench Panel (`.utility-panel`)**:
   - `width: 480px; min-width: 480px; max-width: 480px;` on desktop viewports.
   - Calibrated directly to the 3-column coordinate input and output rows in the **Point** tab.
2. **Asymmetrical 3-Column Coordinate Grid (`.coord-grid-3`)**:
   - `grid-template-columns: minmax(0, 1.35fr) minmax(0, 1.35fr) minmax(0, 0.8fr);`
   - Columns 1 & 2 (`1.35fr` each) hold planimetric/geographic coordinates (Northing/Latitude and Easting/Longitude), guaranteeing that 16-character unspaced DMS strings (`45°50'28.44420"`) fit comfortably on a single line without horizontal scrolling or clipping.
   - Column 3 (`0.8fr`) holds Altitude/Height ($H$ or $h$).
3. **Corner Radii Hierarchy**:
   - `--radius-xs: 4px` (badges, status pills)
   - `--radius-sm: 6px` (segmented buttons, dropdown items)
   - `--radius-md: 6px` (coordinate inputs, outputs, action buttons)
   - `--radius-lg: 8px` (studio cards, floating panels, map docks)
4. **Global Configuration Section Card (`.shared-config-bar`)**:
   - Framed elevated card (`margin: 12px 16px 10px; padding: 10px 14px 10px; border-radius: var(--radius-lg); background-color: var(--bg-surface-elevated); border: 1px solid var(--border-medium);`) containing the discrete `SETTINGS` landmark, optimized into 3 compact horizontal rows (`grid-template-columns: 125px 1fr`):
     - **Row 1**: `Operation Direction` label $\to$ 28px sliding toggle (`ETRS89 to Stereo70` vs `Stereo70 to ETRS89`).
     - **Row 2**: `Coordinate Order` label $\to$ 28px sliding toggle (`NE(H) / LatLon(h)` vs `EN(H) / LonLat(h)`).
     - **Row 3**: `Angular Unit` label $\to$ compact dropdown (`.config-unit-select`, 95px) $\to$ adjacent `Example` label and copyable coordinate badge.
5. **Three-Card Workbench Architecture (`PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)**:
   - Organizes the converters into 3 dedicated section cards (`.studio-card`), a central active separator button, and a decoupled bottom action dock:
     - **Card 1: Input (`INPUT`)**: Point input grid, multi-line batch text area, or drag-and-drop file dropzone with "Ignore header / first line" toggle (white background in light mode, pitch black in dark mode).
     - **Active Separator**: Full-width primary `Transform` action button (`.transform-action-row`, 34px height) bridging input data entry with calculated results.
     - **Card 2: Output (`OUTPUT`)**: Readonly outputs (grey background in both light and dark themes) with top-right copy action, Quasigeoid Height info chip (`.info-field-badge`), export buttons (`CSV`, `GeoJSON`, `DXF`), or transformed file downloads (`.file-download-btn`, `Download All`).
     - **Card 3: Report (`REPORT`)**: Validation status card featuring sentence-case discrete labels (`Coordinates`, `Transformation`), grey non-editable 2D/3D status pills (`.status-pill-neutral`), and semantic quality indicators without icon symbols.
     - **Bottom Actions Dock (`.workbench-actions-dock`)**: `Center on Map` and `Zoom to Extents` secondary outlined buttons (26px height, `.btn-secondary-action`) anchored at the bottom of the workbench pane.


### 3.3 Angle Formatting Display Rules
- **Unspaced DMS**: `DD°MM'SS.SSSSS"` (e.g., `46°00'00.00000"`), with zero whitespace after `°` or `'`.
- **Degrees**: `DD.DDDDDDDD` (8 decimal places).
- **Radians**: `R.RRRRRRRRRR` (10 decimal places).
- **Pure Numeric Readouts**: Input and output boxes display pure numbers only (never append `m`, `°`, or `rad` inside the text field). Cleared or untransformable output boxes display empty fields (`""`) to prevent displaying misleading numbers.

### 3.4 Responsive Mobile Breakpoints & Dockable Workbench Invariants
1. **Desktop Left Docking (`@media (min-width: 769px)`)**:
   - The 480px left workbench can be collapsed to the left (`.utility-panel.is-collapsed { margin-left: -480px; }`) via the vertical divider toggle button (`#workbench-dock-toggle`) or keyboard shortcut (`Ctrl+\` / `Alt+W`).
   - The Leaflet map dynamically occupies 100% of the browser window without requiring OS fullscreen mode.
   - The collapse state is persisted in `localStorage['pytransdat_workbench_collapsed']`.
2. **Mobile Bottom Sheet (`@media (max-width: 768px)`)**:
   - The Leaflet map occupies the full screen (`position: absolute; inset: 0;`).
   - The workbench converts into a Google Maps-style bottom sheet (`.utility-panel`) with a grab handle bar (`.bottom-sheet-handle-bar`).
   - Supports 3 ergonomic states: `sheet-peek` (~72px peek bar at bottom), `sheet-half` (45% height, default), and `sheet-full` (85vh height).
   - Tapping the handle cycles states; touch swipe gestures track smoothly. Focusing coordinate inputs automatically elevates the sheet to `sheet-full` to prevent mobile virtual keyboards from obscuring data fields.
3. **Compact Phone 2-Row Coordinate Grid (`@media (max-width: 540px)`)**:
   - `.coord-grid-3` switches from 3 columns to an ergonomic 2-row layout: Row 1 = Coordinates 1 & 2 (`1fr 1fr`), Row 2 = Coordinate 3 (Height / Z) spanning full width underneath.
   - Guarantees 17-character DMS coordinates (`45°50'28.44420"`) fit with zero horizontal clipping.
   - Form inputs and select dropdowns enforce `font-size: 16px !important;` on mobile ($\le 768\text{px}$) to prevent iOS Safari auto-zoom on focus.
4. **Mobile Bottom Modal Dialogs**:
   - Geodetic Help, Feedback, and Grid Metadata floating windows render as responsive full-width bottom modal sheets on mobile viewports ($\le 768\text{px}$).

---

## 4. Header Architecture & Navigation Invariants (`Header.ts`)

The top navigation bar (`52px` fixed height) serves as the unified application chrome:

1. **Single-Line Compactness (`height: 52px`)**:
   - Preserves a strict `52px` single-line profile across desktop viewports, organized into Left (Brand + Repository Attribution + Primary View Tabs) and Right (Help Menu + Language Selector Menu + Theme Toggle) clusters.
2. **Brand Identity & Centralized Naming**:
   - The application title is centrally governed by `APP_CONFIG.appName` (`VITE_APP_NAME` fallback to `PyTransDatRO`) and synchronized with the browser tab `document.title`.
   - The brand logo (`36px × 36px` inline SVG inside `.brand-logo-wrap`) sits next to a stacked text column where the lower metadata subrow (`v1.0.0-beta` badge on the left, `7.5px` semantic `.health-dot` on the right) spans `width: 100%; justify-content: space-between;` to match the exact width of the title.
   - The `.health-dot` is a **status-only** indicator with native `title` hover tooltips (green for healthy, red for offline) with no click handler; live API outage reporting is decoupled to the global warning banner.
3. **2-Line Repository Attribution (`powered by` / `PyTransDatRO`)**:
   - Located in the **left header zone** immediately to the right of the brand block, separated by a `1px solid var(--border-subtle)` left divider, pointing to GitHub:
     - **Line 1 (`.repo-link-prefix`)**: `powered by` (lowercase) in smaller, muted typography (`9.5px`, `color: var(--text-dim)`).
     - **Line 2 (`.repo-link-name`)**: `PyTransDatRO` in `12px` semi-bold (`font-weight: 600`, `color: var(--text-muted)`).
4. **Primary View Tabs (`Transform` & `Usage`)**:
   - Top-level application navigation (`.header-tab-btn`) styled at `13.5px`, height `32px`, `font-weight: 500` (inactive) / `font-weight: 600` (active).
   - On hover, displays a discrete neutral gray pill background (`--tab-hover-bg`: `rgba(177, 186, 196, 0.12)` in dark mode, `rgba(208, 215, 222, 0.38)` in light mode) with a clean 10px vertical gap above the bottom header line, mirroring GitHub's UnderlineNav.
   - On active, renders a crisp neutral gray `2.5px` bottom indicator line (`--tab-indicator-color`: `#d0d7de` in dark mode, `#1f2328` in light mode) flush with the header boundary.
5. **External Link Indicator Rule (Menu Items Only)**:
   - The oblique external link arrow (`↗` SVG icon) applies **strictly to external items inside dropdown menus** (such as ROMPOS portal, RomGEO desktop, etc.).
   - Top-level header brand and attribution links (`powered by PyTransDatRO`) **never** display the oblique `↗` arrow, keeping the primary header chrome clean and uncluttered.
6. **Help Menu Structure & `RomGEO` Capitalization**:
   - Contains PyTransDatRO guides/diagnostics and official ANCPI/ROMPOS resource links.
   - All references to the official ANCPI transformation suite and API use capital **`GEO`** (**`RomGEO`**).
7. **Action Cluster & Language Selector Menu (`#header-lang-details`)**:
   - The right action cluster flows naturally: `[Help ▾] [RO ▾] [Theme ⊙]` without artificial vertical dividers separating the menus.
   - The language selector uses a clean text-only trigger (`RO ▾` / `EN ▾`), eliminating decorative globe icons or national flags, maintaining visual quietness.
8. **Brand Logo Reset vs. Brand Title Navigation**:
   - **Brand Title (`.brand-title-text`)**: Standard SPA route activation preserving active coordinate inputs.
   - **Brand Logo (`#header-brand-logo-wrap`)**: Explicit first-load reset action purging all inputs, outputs, batch text, file queues, map markers, and resetting to Romania extents.

---

## 5. Dual-Theme Brand Logo (`Map Projection Graticule`)

The official application brand logo is the **Map Projection Graticule** (`cartography-2`), defined in `webapp/src/assets/logos.ts` and `webapp/public/favicon.svg`.
- **Cartographic Visual Concept**: Curved latitude/longitude graticule arcs intersecting a rounded stereographic tangent frame with a high-contrast central origin node.
- **Automatic Dual-Theme Reactivity**:
  - In the application header, `getLogoInlineSvg(36)` renders an inline `<svg class="brand-logo-svg">` parameterized by CSS variables derived via `color-mix()` in `variables.css`. Switching between Dark and Light themes re-themes the logo instantaneously with zero contrast inversion bugs.
  - Simultaneously, `syncFaviconWithTheme(theme)` updates the browser tab `<link rel="icon">` favicon dynamically.

---

## 6. Web GIS Map Widgets & Telemetry Visual Synchronization

### 6.1 Map Overlay Widgets (`webapp/src/style/map.css`)
All floating widgets inside the Leaflet map canvas (`.map-controls-dock`, `.map-legend-dock`, `.map-scale-dock`, `.map-layers-dock`, and `.geodetic-popup`) are bound to the active theme tokens:
- **Zero Hardcoded Slate/Navy**: Map widgets use `var(--bg-glass)`, `var(--bg-surface)`, `var(--border-medium)`, and `var(--radius-md)` / `var(--radius-lg)` so their background blackness and corner curvature match the left workbench panel in both Dark and Light modes.
- **Segmented GIS Graphic Scale Bar (`.map-scale-dock`)**: High-contrast alternating monochrome segments (`--scale-seg-dark` / `--scale-seg-light`) positioned beneath the persistent top-right Legend dock.

### 6.2 Telemetry Charts & Spatial Matrix (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`)
1. **Neutral Axes & Labels**: Both X and Y axes, tick marks, and axis title labels use neutral grey typography (`var(--text-muted)` / `var(--border-subtle)`), never saturated series colors.
2. **Visual Symbology Brightness Hierarchy**:
   - **Calls (Invocations)**: Lighter, luminous/pastel tints of the source hue family.
   - **Points (Volume)**: Darker, deeper tones of the same hue family.
   - **Total Curve**: High-visibility dashed curve (`stroke-dasharray="7 5"`, `2.6px`) in neutral theme text (`var(--text-main)`).
3. **2-Line Grouped Symbology Legend**: Line 1 groups active sources (`[Source]: [Calls] Calls, [Points] Points`); Line 2 isolates `Total`.
4. **Theme-Reactive Spatial Matrix Canvas**: `SpatialMatrixChart.ts` observes `data-theme` on `document.documentElement` via `MutationObserver`, switching ghost cells (`rgba(255,255,255,0.04)` in dark mode vs `rgba(0,0,0,0.05)` in light mode), hover outlines, and the 5-tier color ramp automatically.
