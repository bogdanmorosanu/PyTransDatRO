# PyTransDatRO Web Application & Map — Design System Guide

This document is the authoritative visual design manual for the **PyTransDatRO** web application (`webapp/`). It defines the core design philosophy, the **GitHub Primer neutral Dark & Light theme**, the **Subtle Surface + Border (`#8`) primary action button**, typography, structural layout invariants, the **Component Catalogue (`C1`–`C6`)**, header navigation rules, Web GIS map widget styling, and the **Map Projection Graticule dual-theme brand logo**.

> [!NOTE]
> Purely functional, algorithmic, and data-processing specifications (e.g., coordinate regex validation, byte-stream character decoding, and REST API payload contracts) live in [`docs/web_app.md`](../docs/web_app.md). This guide focuses strictly on **visual architecture, UI ergonomics, and styling invariants**.

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
   - **Zero superfluous emojis or decorative icons**: Tab labels, primary action buttons, and segmented options never contain informal emojis (no ⚡, ⬍, ⬄, 🧭). Clean typography and precise vector geometry communicate authority.
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

All visual styling is governed by semantic CSS custom properties in [`webapp/src/style/variables.css`](src/style/variables.css) across Dark (`data-theme="dark"`) and Light (`data-theme="light"`) modes.

### 2.1 Core Surface, Border & Typography Tokens

| Token | Dark Theme (`data-theme="dark"`) | Light Theme (`data-theme="light"`) | Architectural Role |
| :--- | :--- | :--- | :--- |
| `--bg-app` | `#0d1117` | `#f6f8fa` | Application canvas background. |
| `--bg-surface` | `#161b22` | `#ffffff` | Primary cards, header bar, and left workbench panel. |
| `--bg-surface-elevated` | `#21262d` | `#f3f4f6` | Elevated controls, segmented trays, and secondary buttons. |
| `--bg-surface-input` | `#090c10` | `#f6f8fa` | Recessed coordinate inputs, batch textareas, and inactive toggle trays. |
| `--bg-surface-hover` | `#30363d` | `#eaeef2` | Hover states on interactive buttons and menu items. |
| `--border-subtle` | `#21262d` | `#d0d7de` | Card borders, section dividers, and inactive toggle borders. |
| `--border-medium` | `#30363d` | `#afb8c1` | Input borders, elevated control outlines, and active toggle borders. |
| `--border-focus` | `#58a6ff` | `#0969da` | Keyboard and input focus ring accent. |
| `--text-main` | `#e6edf3` | `#1f2328` | Primary typography, active toggle labels, and coordinate values. |
| `--text-muted` | `#8b949e` | `#59636e` | Secondary labels and inactive segmented options. |
| `--text-dim` | `#6e7681` | `#6e7781` | Micro-section headers and disabled/inactive filter text. |
| `--tab-indicator-color` | `#f78166` | `#fd8c73` | GitHub Primer coral active tab underline indicator. |

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
   - This width is an architectural design invariant calibrated directly to the 3-column coordinate input and output rows in the **Point** tab.
2. **Asymmetrical 3-Column Coordinate Grid (`.coord-grid-3`)**:
   - `grid-template-columns: minmax(0, 1.35fr) minmax(0, 1.35fr) minmax(0, 0.8fr);`
   - Columns 1 & 2 (`1.35fr` each) hold planimetric/geographic coordinates (Northing/Latitude and Easting/Longitude), guaranteeing that 16-character unspaced DMS strings (`45°50'28.44420"`) fit comfortably on a single line without horizontal scrolling or clipping.
   - Column 3 (`0.8fr`) holds Altitude/Height ($H$ or $h$).
3. **Corner Radii Hierarchy**:
   - `--radius-xs: 4px` (badges, status pills)
   - `--radius-sm: 6px` (segmented buttons, dropdown items)
   - `--radius-md: 6px` (coordinate inputs, outputs, action buttons)
   - `--radius-lg: 8px` (studio cards, floating panels, map docks)

### 3.3 Angle Formatting Display Rules
- **Unspaced DMS**: `DD°MM'SS.SSSSS"` (e.g., `46°00'00.00000"`), with zero whitespace after `°` or `'`.
- **Degrees**: `DD.DDDDDDDD` (8 decimal places).
- **Radians**: `R.RRRRRRRRRR` (10 decimal places).
- **Pure Numeric Readouts**: Input and output boxes display pure numbers only (never append `m`, `°`, or `rad` inside the text field). Cleared or uncomputed output boxes display an em-dash (`—`) to prevent vertical layout shift.

---

## 4. Full Component Catalogue (`C1`–`C6`)

Every interactive and informational element in the web application belongs to one of 6 canonical categories (**24 component types total**):

### Category C1: Primary Toggles & Mode Switchers (Clear ON/OFF State Contrast)
- **`C1.1` Sliding Direction Toggle (`.sliding-toggle`, `.sliding-toggle-glider`, `.sliding-toggle-btn`)**:
  - Top-level transformation direction switcher (`ETRS89 ➔ Stereo 70` vs `Stereo 70 ➔ ETRS89`) in Tier 1 of the Shared Config Bar.
  - Uses a recessed `--bg-surface-input` track and a smooth CSS-animated glider (`width: calc(50% - 4px)`, `--btn-primary-bg` fill, `1px solid var(--btn-primary-border)` outline) so the active direction is immediately distinguishable from the dimmed (`opacity: 0.65`) inactive option.
- **`C1.2` Segmented Control (`renderSegmentedControl`, `.segmented-control`, `.segmented-btn`)**:
  - Compact two-option radio switcher used for Coordinate Order (`NE(H)` vs `EN(H)`) in Tier 2 of the Shared Config Bar, featuring an elevated `--btn-primary-bg` active pill with a `1px solid var(--btn-primary-border)` border.
- **`C1.3` Mutually Exclusive Telemetry Segmented Tray (`.telemetry-btn-toggle-group`, `.telemetry-toggle-btn`)**:
  - Enclosed single-choice pill tray used in the Usage dashboard for Time Presets (`All Time`, `12M`, `YTD`, `30D`, `7D`), Granularity (`Daily`, `Weekly`, `Monthly`), and Spatial Scaling (`Logarithmic`, `Linear`).
- **`C1.4` Standalone Multi-Choice Filter Toggles (`.telemetry-individual-toggle-group`, `.telemetry-toggle-btn.individual`, `#comp-toggle-total`)**:
  - Physically separated toggle buttons (`gap: 8px`, independent borders, no outer tray) used when multiple series can be enabled simultaneously (`SOURCES`: Web App, REST API, Desktop / CLI; `METRIC`: Points, Calls; and computation timeline `#comp-toggle-total`).
  - Features an embedded `6px × 6px` circular ON/OFF status dot (`::before`): **ON** renders a solid accent dot (`--color-primary`), `--btn-primary-bg` surface fill, and `--btn-primary-border` border; **OFF** renders a hollow ring dot, recessed `--bg-surface-input` background, and `opacity: 0.62`.

### Category C2: Navigation & Section Headers
- **`C2.1` Global Header View Switcher Tabs (`.header-nav-tabs`, `.header-tab-btn`, `.header-tab-indicator`)**:
  - Top-level application view tabs (`Transform` vs `Usage`) inside the `52px` header, featuring a `2px` bottom active bar (`--tab-indicator-color`).
- **`C2.2` Workbench Folder Bookmark Tabs (`.tabs-bar`, `.tabs-prefix`, `.tab-btn`)**:
  - Mode switcher (`Transform: Point | Text | Files`) that physically merges with the `.studio-card` surface below (`margin-bottom: -1px`, matching surface background and top border, zero gap).
- **`C2.3` Micro-Section Headers (`.coord-section-label`, `.dropdown-group-title`, `.filter-group-label`)**:
  - Quiet uppercase structural labels (`font-size: 10.5px–11px`, `letter-spacing: 0.4px–0.5px`, `color: var(--text-dim)`).
- **`C2.4` Disclosure Menu Triggers (`renderDropdownMenu`, `.nav-menu-trigger`, `.header-dropdown-menu`)**:
  - Native HTML5 `<details><summary>` flyout menus (used for both the `Help` menu and the status-displaying `Language` selector menu in the main header) with rotating chevron indicator, optional active checkmark (`✓`, `role="menuitemradio"`), and keyboard navigation.

### Category C3: Informational Fields, Badges & Status Indicators
- **`C3.1` Interactive Copyable Format Badge (`renderCopyBadge`, `.format-example-row`, `.format-example-badge`)**:
  - Tier 3 Shared Config Bar readout showing a live clickable angle format sample (`45.84123450` / `45°50'28.44420"` / `0.7997813210`).
- **`C3.2` Dedicated Quasigeoid Metric Field (`.quasigeoid-field`, `.quasigeoid-label`, `.quasigeoid-value`)**:
  - Full-width row below single-point coordinate outputs displaying $\zeta = h - H$ in the technical accent palette (`--color-metric-text`, `--color-metric-bg`, `--color-metric-border`).
- **`C3.3` Two-Section Golden Rule Quality Report Cards (`.quality-report-card`, `.quality-report-section`, `.quality-report-divider-dashed`, `.quality-pill`)**:
  - Unified post-transformation diagnostic card across Point, Text, and Files tabs that eliminates ad-hoc inline `.error-banner` alerts:
    - **Top Section (`.coords-section`)**: Micro-label `Coordinates` (`Coordonate` in RO) displaying input coordinate structure pills: `Total` (for batch/files), `2D`, and `3D`.
    - **Dashed Divider (`.quality-report-divider-dashed`)**: Subtle `1px dashed var(--border-subtle)` line separating input traits from transformation outcomes.
    - **Bottom Section (`.transform-section`)**: Micro-label `Transformation` (`Transformare` in RO) displaying engine geodetic diagnostic outcomes: `Valid` (green), `Outside Grid` (amber), `No Grid Data` (amber), `Invalid Values` (red), and `Invalid Angle` (red).
    - Labels read naturally as `[Section] + [Label]` (e.g., `Coordinates: 3D`, `Transformation: Outside Grid`).
- **`C3.4` Live API Health Dot & Global Offline Alert Banner (`.health-dot.online|.offline`, `.brand-stage-badge`, `#global-offline-banner`)**:
  - Stable `8px × 8px` semantic status indicator and monospace version tag (`v1.0.0-beta`) aligned beneath the brand title. The dot never expands or shifts neighboring header elements.
  - When the backend API is unreachable, a persistent high-visibility `#global-offline-banner` renders at the top of the workbench with a direct issue reporting action (`#global-offline-report-btn`) opening the Feedback dialog.

### Category C4: Input Controls & Output Readout Fields
- **`C4.1` Asymmetrical 3-Column Coordinate Inputs (`.coord-grid-3`, `.coord-input-box`, `.coord-col-label`)**:
  - Monospace numeric text fields with `input.select()` on focus and `1.35fr / 1.35fr / 0.8fr` column proportions.
- **`C4.2` Persistent Read-Only Coordinate Outputs (`.coord-output-box`, `.computed-coords-wrapper`)**:
  - Non-shifting monospace output boxes displaying computed tuples or `—` when idle.
- **`C4.3` Batch Textareas & Streaming File Dropzone (`.batch-textarea`, `.dropzone`, `.file-queue-list`, `.file-queue-item`)**:
  - Multi-line monospace text editors with vertical height resize persistence across sessions via `localStorage['pytransdat_batch_textarea_h']`.
  - Native scrollbar hover cursor detection applying `cursor: default` over the scrollbar gutter and thumb to prevent unwanted text `I-beam` cursor interference.
  - Non-destructive batch error formatting omitting dummy coordinates on unconvertible lines (`PointId, # [reason]`) to prevent downstream GIS failures.
  - File queue and direct client-side downloads without premature completion toasts prior to native OS save dialogs.
- **`C4.4` Native Select Dropdowns & Date Inputs (`renderSelect`, `.form-select`, `.telemetry-date-input`)**:
  - Custom-chevron `<select>` controls and `color-scheme`-aware native `<input type="date">` pickers.
- **`C4.5` Action Buttons (`.btn-primary`, `.btn-secondary`, `.btn-outline`, `.btn-icon`, `renderCopyButton` / `.copy-code-btn`)**:
  - Primary transformation button (`Subtle Surface + Border`), secondary map/export triggers, immediate-action icon buttons (such as the Theme Toggle), and top-right code-block copy buttons.

### Category C5: Floating Windows, Popovers & Web GIS Map Widgets
- **`C5.1` Non-Modal Draggable & Resizable Floating Panels (`panels.css`: `.floating-grid-panel`, `.floating-help-panel`, `.floating-feedback-panel`)**:
  - Desktop-class floating windows for **Grid Metadata**, **PyTransDatRO Help & Benchmark**, and **Feedback** with draggable headers, minimize/close controls, and corner resize handles.
- **`C5.2` Modal Guardrail Dialogs (`.modal-backdrop`, `.modal-card`)**:
  - High-volume point warning prompt (`VolumeWarningModal.ts`) and detailed batch diagnostic table (`ReportModal.ts`).
- **`C5.3` Web GIS Map Docks, Drawers, Popups & Scale Bar (`map.css`)**:
  - Top-Left `.map-controls-dock`, Top-Right `.map-top-right-dock` (Expandable Legend + `.map-custom-scale` GIS graphic scale bar), Bottom-Left `.map-layers-drawer`, Bottom-Center `.map-coords-widget`, and `.geodetic-popup` attribute cards.

### Category C6: Analytics & Telemetry Visualizations
- **`C6.1` KPI Summary Cards (`.telemetry-kpi-grid`, `.telemetry-kpi-card`, `.kpi-value`)**:
  - Aggregate volume cards (`Calls`, `Points`, `Points / Call`, `Active Period`) at the top of the Usage view.
- **`C6.2` Dual-Metric SVG Timeline Chart & Interactive Legend Color Picker (`.timeline-chart-container`, `.chart-controls-deck`, `.chart-legend-bar`, `.interactive-legend-item`, `.legend-color-popover`)**:
  - Interactive multi-source SVG time-series plot with neutral grey axes and high-contrast dual-theme default series colors (`REST API`: `#818cf8` / `#4f46e5`; `Web App`: `#34d399` / `#059669`; `Desktop / CLI`: Coral/Crimson Rose `#fb7185` / `#e11d48`, avoiding low-contrast pale yellow on white backgrounds).
  - Clicking any line item in the legend opens a 10-swatch `.legend-color-popover` (matching the Web Map layer widget UX) that updates the SVG series stroke, tooltip dot, and legend line sample immediately and persists custom colors in `localStorage['pytransdat_timeline_colors']`.
- **`C6.3` 0.1° Spatial Matrix Canvas & Quadrant Cards (`.spatial-matrix-deck`, `.spatial-quadrant-card`, `.spatial-cell-tooltip`)**:
  - Theme-reactive HTML5 `<canvas>` punchcard heatmap of Romania (`0.1° × 0.1°` cells) with NW/NE/SW/SE quadrant cards.

---

## 5. Header Architecture & Global Navigation Invariants (`Header.ts`)

The top navigation bar (`52px` fixed height) serves as the unified application chrome:

1. **Single-Line Compactness (`height: 52px`)**:
   - Preserves a strict `52px` single-line profile across desktop viewports, organized into Left (Brand + Repository Attribution + Primary View Tabs) and Right (Help Menu + Language Selector Menu + Theme Toggle) clusters.

2. **Brand Identity & Metadata Subrow Alignment**:
   - The brand logo (`36px × 36px` inline SVG inside `.brand-logo-wrap`) sits next to a stacked text column where the lower metadata subrow (`v1.0.0-beta` badge on the left, `7.5px` semantic `.health-dot` on the right) spans `width: 100%; justify-content: space-between;` to match the exact width of `PyTransDatRO`.

3. **2-Line Repository Attribution (`powered by` / `PyTransDatRO`)**:
   - Located in the **left header zone** (`.header-repo-link.zone-left-attribution`) immediately to the right of the `PyTransDatRO` brand block, separated by a `1px solid var(--border-subtle)` left divider, pointing to `https://github.com/bogdanmorosanu/PyTransDatRO`:
     - **Line 1 (`.repo-link-prefix`)**: `powered by` (lowercase) in smaller, muted typography (`9.5px`, `color: var(--text-dim)`).
     - **Line 2 (`.repo-link-name`)**: `PyTransDatRO` (preserving `"Py"` in front to match the GitHub repository name) in `12px` semi-bold (`font-weight: 600`, `color: var(--text-muted)`, brightening to `--color-primary-light` on hover).

4. **External Link Indicator Rule (Menu Items Only)**:
   - The oblique external link arrow (`↗` SVG icon, `class="ext-icon"`) applies **strictly to external items inside dropdown menus** (such as the `Help` menu links to the ROMPOS portal, RomGEO desktop downloads, RomGEO 3D Grid, RomGEO Web API, and `/docs` REST API reference).
   - Top-level header brand and attribution links (`powered by PyTransDatRO`) **never** display the oblique `↗` arrow, keeping the primary header chrome clean and uncluttered.

5. **Help Menu Structure & `RomGEO` Capitalization**:
   - The `Help` menu (`#header-help-details`) contains:
     - **Interactive Guides & Diagnostics**: `PyTransDatRO Help & Benchmark`, `Grid Metadata Inspector`, `REST API Docs (/docs)` (`isExternal: true`).
     - **ANCPI / ROMPOS Resources**:
       - `ROMPOS Portal` (`https://rompos.ro`, `isExternal: true`)
       - `RomGEO Desktop (Win)` (`isExternal: true`)
       - `RomGEO Desktop (Mac)` (`isExternal: true`)
       - `RomGEO 3D Grid (.txt)` (`isExternal: true`)
       - `RomGEO Web API` (`https://api.romgeo.ro/api/v1/`, `isExternal: true`) — placed at the bottom of the ANCPI / ROMPOS list.
     - **Support**: `Send Feedback / Issue`.
   - All references to the official ANCPI transformation suite and API use capital **`GEO`** (**`RomGEO`**).

6. **Distinction Between Status-Displaying Menus (`Language`) and Immediate-Action Buttons (`Theme Toggle`)**:
   - **Language Selector Dropdown Menu (`#header-lang-details`)**: Because the language control displays the **current active state** (`RO` or `EN`), it is rendered as a disclosure **Dropdown Menu** (`renderDropdownMenu` with `.nav-menu-trigger`, a `14px` globe SVG icon, the active language code, and a rotating chevron `▾`). Opening the menu displays `Română (RO)` and `English (EN)` (`role="menuitemradio"`, `aria-checked`) with an active checkmark (`✓`, `.menu-check-icon`) next to the currently selected language.
   - **Theme Toggle Button (`#header-theme-btn`)**: Because the theme toggle is an **immediate action** that flips the visual mode according to the sun/moon icon shown, it is rendered as a single-click icon button (`.btn-icon`).

7. **Brand Logo Reset (`resetToFreshState`) vs. Brand Title Navigation**:
   - **Brand Title (`#header-brand-title` / `PyTransDatRO`)**: Standard SPA route activation (`navigateTab('point')` or default view), preserving currently entered coordinates, batch text, files, and layer states.
   - **Brand Logo (`#header-brand-logo` / `.brand-logo-wrap`)**: Explicit First-Load Reset Action. Clicking the cartographic logo resets the application to a pristine state identical to a fresh page load: clears single-point inputs/outputs, batch textareas, file drop queues, map markers, and resets direction to default `ETRS89 ➔ Stereo 70`.

---

## 6. Dual-Theme Brand Logo (`Map Projection Graticule`)

The official application brand logo is the **Map Projection Graticule** (`cartography-2`), defined in [`webapp/src/assets/logos.ts`](src/assets/logos.ts), [`webapp/src/assets/logo.svg`](src/assets/logo.svg), and [`webapp/public/favicon.svg`](public/favicon.svg).
- **Cartographic Visual Concept**: Curved latitude/longitude graticule arcs (`--logo-accent-1` and `--logo-accent-2`) intersecting a rounded stereographic tangent frame (`--logo-grid`) with a high-contrast central origin node (`--logo-node-main`) on a rounded `40×40` tile (`--logo-bg`, `--logo-border`).
- **Automatic Dual-Theme Reactivity**:
  - In the application header, `getLogoInlineSvg(36)` renders an inline `<svg class="brand-logo-svg">` parameterized by CSS variables derived via `color-mix()` in `variables.css`. Switching between Dark and Light themes re-themes the logo instantaneously with zero contrast inversion bugs.
  - Simultaneously, `syncFaviconWithTheme(theme)` updates the browser tab `<link rel="icon">` favicon using `getLogoStandaloneSvg('dark' | 'light')` (`#141416` dark tile vs `#e9eaec` light tile).

---

## 7. Web GIS Map Widgets & Telemetry Visual Synchronization

### 7.1 Map Overlay Widgets (`webapp/src/style/map.css`)
All floating widgets inside the Leaflet map canvas (`.map-controls-dock`, `.map-legend-widget`, `.map-custom-scale`, `.map-layers-drawer`, `.map-coords-widget`, and `.geodetic-popup`) are bound to the active theme tokens:
- **Zero Hardcoded Slate/Navy**: Map widgets use `var(--bg-glass)`, `var(--bg-surface)`, `var(--border-medium)`, and `var(--radius-md)` / `var(--radius-lg)` so their background blackness and corner curvature match the left workbench panel in both Dark and Light modes.
- **Segmented GIS Graphic Scale Bar (`.map-custom-scale`)**: High-contrast alternating monochrome segments (`--scale-seg-dark` / `--scale-seg-light`) positioned beneath the persistent top-right Legend dock.

### 7.2 Telemetry Charts & Spatial Matrix (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`)
1. **Neutral Axes & Labels**: Both X and Y axes, tick marks, and axis title labels use neutral grey typography (`var(--text-muted)` / `var(--border-subtle)`), never saturated series colors.
2. **Visual Symbology Brightness Hierarchy**:
   - **Calls (Invocations)**: Lighter, luminous/pastel tints of the source hue family.
   - **Points (Volume)**: Darker, deeper tones of the same hue family.
   - **Total Curve**: High-visibility dashed curve (`stroke-dasharray="7 5"`, `2.6px`) in neutral theme text (`var(--text-main)`).
3. **2-Line Grouped Symbology Legend**: Line 1 groups active sources (`[Source]: [Calls] Calls, [Points] Points`); Line 2 isolates `Total`.
4. **Theme-Reactive Spatial Matrix Canvas**: `SpatialMatrixChart.ts` observes `data-theme` on `document.documentElement` via `MutationObserver`, switching ghost cells (`rgba(255,255,255,0.04)` in dark mode vs `rgba(0,0,0,0.05)` in light mode), hover outlines, and the 5-tier color ramp automatically.

---

## 8. Checklist for Future UI & Map Changes

Before committing any frontend UI or map changes, verify:
- [ ] Are all colors derived from the 5 foundational tokens (`--text`, `--background`, `--primary`, `--secondary`, `--accent`) and native CSS `color-mix()` variables in `variables.css` (zero dark navy-blue `#090d16` / `#0f172a` / `#1e293b`)?
- [ ] Does `.btn-primary` use `--primary` (`#8c969e` dark / `#636c74` light) with `--background` text for balanced contrast?
- [ ] Is the left workbench panel width maintained at `480px` and paired with the `1.35fr 1.35fr 0.8fr` 3-column coordinate grid?
- [ ] Does the header display the 2-line lowercase `"powered by"` / `"PyTransDatRO"` attribution in the left header zone without an oblique `↗` arrow, while external links inside dropdown menus retain `↗` (`.ext-icon`)?
- [ ] Is the Language selector rendered as a status-displaying Dropdown Menu (`#header-lang-details`) with active checkmarks, distinct from the single-click Theme Toggle icon button (`#header-theme-btn`)?
- [ ] Does clicking the brand logo trigger a clean first-load reset (`resetToFreshState`) while clicking the brand title performs soft SPA tab navigation?
- [ ] Does the Quality Report adhere to the Two-Section layout (`Coordinates` top, dashed divider, `Transformation` bottom) with zero ad-hoc inline `.error-banner` elements?
- [ ] Do map markers preserve their layer color outside the grid with continuous soft pulse animation?
- [ ] Does the Help menu include `RomGEO Web API` (`https://api.romgeo.ro/api/v1/`) at the bottom of the ANCPI / ROMPOS list and use capital `GEO` (`RomGEO`) across all links?
- [ ] Does the `Map Projection Graticule` brand logo use inline CSS-variable SVG rendering (`getLogoInlineSvg`) and `syncFaviconWithTheme` so it adapts cleanly to both Dark and Light themes?
- [ ] Does Quasigeoid Height ($\zeta = h - H$) use the dedicated technical accent metric palette (`--color-metric-text: var(--accent)`) and never status green/red/yellow?
- [ ] Do folder bookmark tabs merge seamlessly with the `.studio-card` container below (zero gap)?
- [ ] Do all 80 automated Vitest unit tests (`npm --prefix webapp test`) and 25 Playwright E2E browser tests (`.\test-e2e.ps1`) pass with 0 errors?
- [ ] Does `npm --prefix webapp run build` compile with 0 TypeScript/bundler errors?
