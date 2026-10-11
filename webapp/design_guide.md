# PyTransDatRO Web Application & Map — Design System Guide

This document defines the core styling principles, visual hierarchy, color palette, typography, layout architecture, and component specifications for the **PyTransDatRO** web application (`webapp/`).

Any future UI additions, enhancements, or component refactors must strictly adhere to the guidelines and invariants specified here.

---

## 1. Core Design Philosophy: The Trinity & Singularity of Purpose

### 1.1 Singularity of Purpose
PyTransDatRO has a dedicated, singular mission: **accurate, reliable geodetic coordinate transformations** between Romania's national projection (Stereo 70 / Black Sea 1975) and European/global standards (ETRS89 / GRS80).
- Avoid non-essential feature bloat, extraneous data validation, or auxiliary utilities that detract from this core purpose.
- Prioritize professional engineering workflows: speed, precision, uncluttered layout, and effortless copy-paste interoperability with CAD and GIS software.

### 1.2 The Trinity Balance

```
                 [ Simplicity ]
                 /            \
                / PyTransDatRO \
               /  Design Center \
        [ Optimal ] ---------- [ Modern ]
```

1. **Simplicity (Clarity over Decoration)**:
   - **Zero superfluous emojis or decorative icons**: Tab labels, action buttons, and segmented options must never contain emojis (no ⚡, ⬍, ⬄, 🧭, etc.). Clean typography communicates authority.
   - **Subtle micro-headers**: Label sections with quiet, high-contrast micro-labels (`Input coordinates (or click on map)`, `Computed coordinates`, `Download As:`).
   - **Visual quietness**: Secondary configuration toggles remain subtle so the primary CTA button commands clear focus.

2. **Optimal (Task Efficiency & Ergonomics)**:
   - **480px left workbench panel**: Provides ample breathing room for multi-column coordinate tuples without encroaching on the map.
   - **Asymmetrical 3-column coordinate layout**: Northing/Lat (`1.35fr`), Easting/Lon (`1.35fr`), and Altitude/Height (`0.8fr`) prevent long DMS strings (`46°00'00.00000"`) from wrapping or clipping.
   - **Focus auto-select (`input.select()`)**: Focusing any coordinate box automatically selects its text for instant replacement.
   - **Top-right code-block copy buttons**: Icon-based copy buttons with instant tooltip feedback (`✓ Copied`) provide frictionless extraction.
   - **Dynamic coordinate order synchronization**: Switching between `NE(H)` and `EN(H)` dynamically swaps Column 1 and Column 2 across inputs, outputs, and clipboard exports.

3. **Modern (Precision Geospatial Engineering)**:
   - Inspired by modern developer and geospatial tools (Linear, Raycast, Mapbox Studio).
   - **Strict semantic separation of colors**: Physical data measurements (e.g., Quasigeoid Height $\zeta$) must never use status colors (green/red/yellow). Status colors are reserved exclusively for system health and grid coverage validity.
   - **High-contrast anchor button**: Primary CTA uses an engineered theme-native anchor surface rather than generic cobalt blue.

---

## 2. Design Tokens & Color System

All styles derive from CSS custom properties defined in [`webapp/src/style/variables.css`](src/style/variables.css).

### 2.1 Theme Palettes

#### Dark Theme (Default)
- **App Background Canvas**: `#090d16` (Deep Obsidian)
- **Card Surface**: `#0f172a` (Slate Surface)
- **Elevated Surface (Active tabs, badges)**: `#1e293b`
- **Input Field Background**: `#0a0f1d`
- **Hover Background**: `#27354f`
- **Subtle Border**: `#1e293b`
- **Medium Border**: `#334155`
- **Focus Border**: `#38bdf8`
- **Text Main**: `#f8fafc`
- **Text Muted**: `#94a3b8`
- **Text Dim**: `#64748b`

#### Light Theme
- **App Background Canvas**: `#f8fafc` (Crisp Porcelain)
- **Card Surface**: `#ffffff`
- **Elevated Surface (Active tabs, badges)**: `#f1f5f9`
- **Input Field Background**: `#ffffff`
- **Hover Background**: `#e2e8f0`
- **Subtle Border**: `#e2e8f0`
- **Medium Border**: `#cbd5e1`
- **Focus Border**: `#0284c7`
- **Text Main**: `#0f172a`
- **Text Muted**: `#475569`
- **Text Dim**: `#94a3b8`

### 2.2 Primary CTA Anchor Button (`.btn-primary`)
Rather than relying on generic saturated blue, the primary action button uses a commanding, theme-native anchor style:
- **Dark Mode**: Crisp off-white / porcelain surface (`#f8fafc`) with deep charcoal text (`#090d16`), bold weight (`600`), and subtle elevation.
- **Light Mode**: Deep slate / charcoal surface (`#0f172a`) with crisp white text (`#f8fafc`).
- **Label**: Pure text `Transform` / `Transformă` (no emojis, no icons).

### 2.3 Semantic & Metric Colors (Strict Isolation)
- **Primary Technical Metric (`--color-metric-text`)**:
  - Dark Theme: `#38bdf8` (Cyan Blue)
  - Light Theme: `#0284c7`
  - Background Tint: `rgba(56, 189, 248, 0.06)`
  - Border Tint: `rgba(56, 189, 248, 0.25)`
  - **Purpose**: Strictly reserved for geodetic data measurements, notably the **Quasigeoid Height** ($\zeta = h - H$).
- **Success (`#10b981`)**: Strictly for verified system health and valid in-grid coordinates (`✓ In Grid`).
- **Warning (`#f59e0b`)**: Strictly for cautionary notices (e.g. `⚠ Outside Grid`, `⚠ No Grid Data`, or large volume warning modals >2,000 points).
- **Danger / Error (`#ef4444`)**: Strictly for syntax errors (`✕ Syntax Error`), network failures, or unparseable input.

> [!IMPORTANT]
> **Data Values Rule**: Quasigeoid height $\zeta$, coordinate shifts $\Delta E, \Delta N$, and coordinate tuples must **never** be rendered in green, red, or yellow. They are physical scalar measurements and must use either neutral typography or the dedicated technical metric cyan palette.

---

## 3. Typography, Curvature Scale & Formatting

### 3.1 Typography
- **Interface Font**: `'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif`
- **Data & Coordinate Font**: `'JetBrains Mono', 'Fira Code', Menlo, Monaco, Consolas, monospace`
  - Used for all coordinate inputs, outputs, Quasigeoid values, format badges, and code blocks.

### 3.2 Curvature Scale (Rounded Radii)
- `--radius-xs` (`4px`): Micro badges, unit tags, corner accents.
- `--radius-sm` (`6px`): Inner sliding toggle gliders, segmented buttons, tooltips.
- `--radius-md` (`8px`): Form inputs, main buttons, folder tabs, segmented trays.
- `--radius-lg` (`12px`): Studio cards, floating map control containers, modal cards.
- `--radius-xl` (`16px`): Large dialogs and overlay sheets.
- `--radius-full` (`9999px`): Circular status pills and map pin markers.

### 3.3 Angle Formatting Conventions
- **Unspaced DMS**: `DD°MM'SS.SSSSS"` (e.g., `46°00'00.00000"`), without whitespace after `°` or `'`. Matches official Romgeo and pure-Python library conventions, preventing copy-paste splitting across input boxes.
- **Degrees**: `DD.DDDDDDDD` (8 decimal places).
- **Radians**: `R.RRRRRRRRRR` (10 decimal places).

---

## 4. Component Architecture & Specifications

### 4.1 Folder Bookmark Navigation Tabs (`.tabs-bar`, `.tab-btn`)
- **Seamless Physical Continuity**: Folder bookmark tabs must seamlessly connect with the active card container below. Zero 16px background gap; the active tab visually and physically merges with the top border of the card surface below.
- **Prefix Label**: Subtle prefix `Transform:` (`Transformă:`) precedes tab buttons.
- **Concise Labels**: Text-only, title case: `Point`, `Text`, `Files` (RO: `Punct`, `Text`, `Fișiere`).
- **Header Direct Links**: `Grid Info` and `Geodetic Help` are direct text links in the top navbar (`Header.ts`), keeping the transformation tabs strictly focused on input modes.

### 4.2 3-Tier Shared Transformation Configuration Bar
Located at the top of the left workbench panel, unifying configuration across all tabs:
1. **Tier 1: Sliding Direction Toggle**:
   - Pure CSS sliding toggle with 6px rounded-rectangle glider (`var(--radius-sm)`).
   - Default direction: `ETRS89 ➔ Stereo 70` (first button); `Stereo 70 ➔ ETRS89` (second button).
2. **Tier 2: 50/50 Balanced Split**:
   - Left: `Angular unit: [select]` (Native dropdown with `Degrees`, `DMS`, `Radians`).
   - Right: `Coordinate order: [ NE(H) | EN(H) ]` (Segmented toggle with vertical divider).
   - Sentence case: `Angular unit:` and `Coordinate order:` (never uppercase `Angular Unit:`).
3. **Tier 3: Copyable Angle Format Example**:
   - Dynamic badge updating reactively with chosen unit (`45.84123450`, `45°50'28.44420"`, or `0.7997813210`).
   - Clicking the badge highlights and copies strictly the numeric angle value to clipboard (excluding label text) with toast confirmation.

### 4.3 Single Point Transformation (`PointConverter.ts`)
- **Projection Center Defaults**: Initial input coordinates set to Romania's projection center origin:
  - Latitude $\phi = 46.0^\circ\text{ N}$ (`46°00'00.00000"` in DMS)
  - Longitude $\lambda = 25.0^\circ\text{ E}$ (`25°00'00.00000"` in DMS)
  - Ellipsoidal Height $h = 100.0\text{ m}$
- **Asymmetrical 3-Column Inputs & Outputs Grid**:
  - Col 1 & 2: `minmax(0, 1.35fr)` for Northing/Lat and Easting/Lon to fit long DMS strings without clipping.
  - Col 3: `minmax(0, 0.8fr)` for Altitude/Height.
- **Auto 2D/3D Detection**:
  - Empty 3rd coordinate $\to$ 2D mode (Quasigeoid Height displays `—`).
  - Provided 3rd coordinate $\to$ 3D mode (computes $\zeta = h - H$).
- **Output Hygiene**:
  - Pure numeric values only (strip `m`, `°`, `rad` from text boxes).
  - Persistent output boxes display `—` when cleared or invalid, preventing vertical layout shift.
  - Top-right code-block style copy icon button on computed coordinates header with instant `✓ Copied` feedback.
  - Auto-select on focus: `input.select()` across all coordinate inputs.
- **Dedicated Quasigeoid Height Field**:
  - Dedicated field on its own row below coordinates.
  - Styled in technical metric cyan-blue (`--color-metric-text`, `--color-metric-bg`).
  - Standardized label: `Quasigeoid Height (ζ = h - H)`.
- **Paste Ergonomics**:
  - Single-DMS paste guard (`isSingleDmsString`): prevents pasting a DMS string into box 1 from erroneously splitting at spaces into boxes 2 and 3.

### 4.4 Batch Text Transformation (`BatchConverter.ts`)
- **Strict Vertical Layout**:
  - Input textarea $\to$ Primary Transform button $\to$ Output textarea $\to$ Quality Report $\to$ Download actions $\to$ Map actions.
  - No redundant pre-calculation diagnostic cards.
- **Strict Numeric Validation**:
  - `isValidDecimal` regex prevents non-numeric typos (e.g. `5001a50.000`) from being partially evaluated.
- **Strict 1-to-1 Line Correspondence**:
  - Input row $N$ maps directly to output row $N$. Empty lines, comments, and invalid rows maintain their exact position.
- **Human 1-Based Indexing**:
  - Unnamed points are referred to as `Point 1`, `Point 2`, etc.
- **Out-of-Grid Suppression**:
  - Suppress bogus approximate transformed coordinates; emit only `[Point_ID], # Out of grid`.
- **Relocated Export Actions**:
  - Grouped below output textarea under `Download As:` (`CSV`, `GeoJSON`, `DXF`).

### 4.5 Bulk File Processing (`FileUpload.ts`)
- **Interactive Queue Pipeline**:
  - Dropzone $\to$ Queued files list (file name, size, status, individual `✕` remove button) $\to$ `Transform Files` manual trigger $\to$ Cumulative Quality Report $\to$ Download actions (`Download All` / individual links) $\to$ Map actions.
  - Raw 5-line text preview removed.
- **Robust Line Streaming Parser**:
  - Auto-detects Point ID column (3 vs 4 tokens).
  - Parses DMS angles with sexagesimal radian conversion.
  - Decodes byte streams using UTF-8 with automatic CP1252 / Latin-1 fallback for degree symbols.

### 4.6 Web GIS Map Layer Architecture & Decoupled Viewports (`MapView.ts`)
- **3 Independent Leaflet LayerGroups**:
  - `pointLayer`, `textLayer`, `filesLayer`, independently managed and toggleable in Leaflet layer control.
- **Decoupled Viewports (No Auto-Panning)**:
  - Transformations never automatically pan or jump the user's map viewport.
- **Dual Map Navigation Actions**:
  - Placed in the bottom-right action footer of each tab:
    - **Center on Map**: Pans viewport to centroid of valid in-grid coordinates, preserving current zoom.
    - **Zoom to Extents**: Pans and calls `fitBounds` to encompass all points, or zooms to level 14 for a single point.
- **Boundary-Safe Guard**:
  - Strictly guarded within the Romanian grid envelope (`[43.5, 20.2]` to `[48.3, 29.8]` and coarse SPG boundaries from `gridBounds.ts`).
  - Never pan or zoom away from Romania on out-of-grid coordinates.
- **Smart Map Click**:
  - In `ETRS89ToStereo70`: populates clicked lat/lon formatted in active angular units.
  - In `Stereo 70ToETRS89`: reverse-calculates clicked point into Stereo 70 projection meters.

### 4.7 Standardized Quality Report Terminology
Concise, bilingual status badges across all tabs:
- `In Grid` (`În grid`)
- `Outside Grid` (`În afara gridului`)
- `No Grid Data` (`Fără date în grid`)
- `Syntax Error` (`Eroare sintaxă`)

---

## 5. Checklist for Future UI & Map Changes

Before committing any frontend UI or map changes, verify:
- [ ] Are all tabs, buttons, and segmented options free of informal emojis and decorative icons?
- [ ] Is the left workbench panel width maintained at `480px`?
- [ ] Do folder bookmark tabs merge seamlessly with the card container below (zero gap)?
- [ ] Does the primary action button use the high-contrast theme-native anchor style (`.btn-primary`)?
- [ ] Are secondary options (direction, units, order) contained in the 3-tier shared config bar?
- [ ] Are all coordinates and numerical values rendered in monospace typography (`JetBrains Mono`)?
- [ ] Are output text boxes pure numbers (no unit suffixes `m`, `°`, `rad` inside inputs)?
- [ ] Does Quasigeoid Height use the dedicated technical cyan metric palette (`--color-metric-text`) and never green/red?
- [ ] Are map viewports decoupled from coordinate transforms (no auto-panning)?
- [ ] Do all tabs provide separate `Center on Map` and `Zoom to Extents` actions, strictly guarded within Romania's envelope?
- [ ] Does `npm --prefix webapp run build` compile cleanly with 0 TypeScript/bundler errors?
