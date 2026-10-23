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

### 3.3 Angle Formatting Display Rules
- **Unspaced DMS**: `DD°MM'SS.SSSSS"` (e.g., `46°00'00.00000"`), with zero whitespace after `°` or `'`.
- **Degrees**: `DD.DDDDDDDD` (8 decimal places).
- **Radians**: `R.RRRRRRRRRR` (10 decimal places).
- **Pure Numeric Readouts**: Input and output boxes display pure numbers only (never append `m`, `°`, or `rad` inside the text field). Cleared or uncomputed output boxes display an em-dash (`—`) to prevent vertical layout shift.

---

## 4. Header Architecture & Navigation Invariants (`Header.ts`)

The top navigation bar (`52px` fixed height) serves as the unified application chrome:

1. **Single-Line Compactness (`height: 52px`)**:
   - Preserves a strict `52px` single-line profile across desktop viewports, organized into Left (Brand + Repository Attribution + Primary View Tabs) and Right (Help Menu + Language Selector Menu + Theme Toggle) clusters.
2. **Brand Identity & Metadata Subrow Alignment**:
   - The brand logo (`36px × 36px` inline SVG inside `.brand-logo-wrap`) sits next to a stacked text column where the lower metadata subrow (`v1.0.0-beta` badge on the left, `7.5px` semantic `.health-dot` on the right) spans `width: 100%; justify-content: space-between;` to match the exact width of `PyTransDatRO`.
3. **2-Line Repository Attribution (`powered by` / `PyTransDatRO`)**:
   - Located in the **left header zone** immediately to the right of the brand block, separated by a `1px solid var(--border-subtle)` left divider, pointing to GitHub:
     - **Line 1 (`.repo-link-prefix`)**: `powered by` (lowercase) in smaller, muted typography (`9.5px`, `color: var(--text-dim)`).
     - **Line 2 (`.repo-link-name`)**: `PyTransDatRO` in `12px` semi-bold (`font-weight: 600`, `color: var(--text-muted)`).
4. **External Link Indicator Rule (Menu Items Only)**:
   - The oblique external link arrow (`↗` SVG icon) applies **strictly to external items inside dropdown menus** (such as ROMPOS portal, RomGEO desktop, etc.).
   - Top-level header brand and attribution links (`powered by PyTransDatRO`) **never** display the oblique `↗` arrow, keeping the primary header chrome clean and uncluttered.
5. **Help Menu Structure & `RomGEO` Capitalization**:
   - Contains PyTransDatRO guides/diagnostics and official ANCPI/ROMPOS resource links.
   - All references to the official ANCPI transformation suite and API use capital **`GEO`** (**`RomGEO`**).
6. **Distinction Between Status-Displaying Menus (`Language`) and Immediate-Action Buttons (`Theme Toggle`)**:
   - **Language Selector Dropdown Menu (`#header-lang-details`)**: Displays current active state (`RO` or `EN`) via disclosure Dropdown Menu (`renderDropdownMenu`) with globe icon and active checkmarks (`✓`).
   - **Theme Toggle Button (`#header-theme-btn`)**: Immediate single-click icon button (`.btn-icon`) flipping dark/light mode according to the sun/moon icon shown.
7. **Brand Logo Reset vs. Brand Title Navigation**:
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
