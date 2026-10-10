# PyTransDatRO Web Application & Map — Design System Guide

This document defines the core styling principles, visual hierarchy, color palette, typography, and component specifications for the **PyTransDatRO** web application (`webapp/`).

Any future additions, UI enhancements, or component refactors must adhere to the design rules specified here.

---

## 1. Core Design Philosophy: The Trinity

The visual identity and user experience of PyTransDatRO sit at the intersection of three non-negotiable principles:

```
          [ Simplicity ]
         /              \
        /   PyTransDatRO \
       /   Design Center  \
[ Optimal ] ------------ [ Modern ]
```

1. **Simplicity (Clarity over Decoration)**:
   - **Zero superfluous visual artifacts**: Avoid decorative emojis or icons next to standard text tabs/menus (`Punct Unic`, `Lot Coordonate`, etc.). Clean typography communicates authority.
   - **Focused visual hierarchy**: Secondary selectors (direction, coordinate order, angle units) must remain visually quiet so the primary call-to-action stands out effortlessly.
   - High readability for numerical coordinate data without visual clutter.

2. **Optimal (Task Efficiency & Screen Ergonomics)**:
   - **Ergonomic coordinate handling**: Inputs and outputs must comfortably accommodate multi-column coordinate strings (Northing, Easting, Elevation, Point ID) without awkward line wrapping or horizontal squishing.
   - **Tactile segmented controls**: Options and toggles must feature subtle vertical separation dividers and crisp elevated active states so users can instantly parse where one option ends and the next begins.
   - Flexible layouts: Provide adaptable views (e.g., Stacked vs. Side-by-side for batch processing) based on user workflow needs.

3. **Modern (Precision Engineering Standard)**:
   - Inspired by modern geospatial and developer tools (e.g., Linear, Raycast, Mapbox Studio).
   - Cohesive curvature system using harmonized rounded corners across all components.
   - **Strict semantic separation of colors**: Data measurements (such as Quasigeoid undulation $\zeta$) must never use status colors (like green or red) that mislead the user regarding correctness or health.

---

## 2. Design Tokens & Color System

All styles derive from CSS custom properties defined in [`webapp/src/style/variables.css`](variables.css).

### 2.1 Theme Palettes

#### Dark Theme (Default)
- **App Background Canvas**: `#090d16` (Deep Obsidian)
- **Card Surface**: `#0f172a` (Slate Surface)
- **Elevated Surface (Active pills, cards)**: `#1e293b`
- **Input Field Background**: `#0a0f1d`
- **Hover Background**: `#27354f`
- **Subtle Border**: `#1e293b`
- **Medium Border**: `#334155`
- **Focus Border**: `#3b82f6`
- **Text Main**: `#f8fafc`
- **Text Muted**: `#94a3b8`
- **Text Dim**: `#64748b`

#### Light Theme
- **App Background Canvas**: `#f8fafc` (Crisp Porcelain)
- **Card Surface**: `#ffffff`
- **Elevated Surface (Active pills, cards)**: `#f1f5f9`
- **Input Field Background**: `#ffffff`
- **Hover Background**: `#e2e8f0`
- **Subtle Border**: `#e2e8f0`
- **Medium Border**: `#cbd5e1`
- **Focus Border**: `#2563eb`
- **Text Main**: `#0f172a`
- **Text Muted**: `#475569`
- **Text Dim**: `#94a3b8`

### 2.2 Semantic & Status Colors (Strict Isolation)
- **Primary Action (CTA)**: `#2563eb` (Hover: `#1d4ed8`, Glow: `rgba(37, 99, 235, 0.22)`). Reserved exclusively for main action buttons (e.g., `⚡ Transformă`).
- **Success (`#10b981`)**: Strictly for verified system health, successful file exports, and valid point counts.
- **Warning (`#f59e0b`)**: Strictly for large batch warning modals (>2,000 points) or near-boundary notifications.
- **Danger / Error (`#ef4444`)**: Strictly for points outside the official TransDatRO grid bounding box, network failures, or invalid syntax.

> [!IMPORTANT]
> **Data Values Rule**: Quasigeoid undulation $\zeta$, coordinate shifts $\Delta E, \Delta N$, and coordinate tuples must **never** be rendered in green, red, or yellow. They are physical scalar measurements and must be presented in neutral monospace typography (`var(--text-main)`).

---

## 3. Typography & Curvature Scale

### 3.1 Typography
- **Interface Font**: `'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif`
- **Data & Numerical Font**: `'JetBrains Mono', 'Fira Code', Menlo, Monaco, Consolas, monospace`
  - Used for all coordinate inputs, outputs, azimuths, tolerances, point IDs, and code blocks.

### 3.2 Curvature Scale (Rounded Radii)
- `--radius-xs` (`4px`): Micro badges, unit tags, corner accents.
- `--radius-sm` (`6px`): Inner segmented buttons, tooltips.
- `--radius-md` (`8px`): Form inputs, main buttons, tab buttons, segmented control trays.
- `--radius-lg` (`12px`): Studio cards, floating map control containers, modal cards.
- `--radius-xl` (`16px`): Large dialogs and overlay sheets.
- `--radius-full` (`9999px`): Circular status pills and map pin markers.

---

## 4. Component Guidelines

### 4.1 Navigation Tabs (`.tabs-bar`, `.tab-btn`)
- Use plain, professional text labels without emoji icons (e.g., `Punct Unic`, `Lot Coordonate`, `Fișier`, `Statistici`).
- Each tab button has a rounded geometry (`border-radius: var(--radius-md)`).
- Active state: uses an elevated surface (`--bg-surface-elevated`) with subtle border and text contrast, rather than flat bottom border lines.

### 4.2 Segmented Option Controls (`.segmented-control`, `.segmented-btn`)
- Used for mutually exclusive settings:
  - Transformation Direction (`Stereo 70 -> ETRS89` / `ETRS89 -> Stereo 70`)
  - Dimension (`2D` / `3D`)
  - Angle Units (`Grade zecimale` / `DMS` / `Radiani`)
  - Coordinate Order (`N, E` / `E, N`)
- **Container**: Recessed track (`var(--bg-surface-input)`), rounded (`var(--radius-md)`).
- **Dividers**: Inactive adjacent buttons must be visually separated by subtle vertical divider lines (`1px solid var(--border-subtle)`).
- **Active State**: Must use the elevated neutral surface (`var(--bg-surface-elevated)`), high-contrast text (`var(--text-main)`), and a subtle border. **Do not use saturated blue for segmented options** to avoid visual competition with the primary action button.

### 4.3 Batch Coordinates Converter (`BatchConverter.ts`)
- **Workbench Width**: The utility panel grid column must maintain a minimum width of `560px` to `640px` to avoid horizontal coordinate truncation.
- **Layout Modes**:
  - **Stacked (`⬍ Suprapus`)**: Recommended default for dense geodetic coordinate sets (Input on top, Output below), giving 100% width to coordinate tuples.
  - **Split (`⬄ Alăturat`)**: Optional side-by-side mode available via the header layout toggle button for wider desktop setups.
- Both textareas must use `font-family: var(--font-mono)` with `font-size: 12px`.

### 4.4 Interactive Leaflet Web GIS Map (`MapView.ts`)
- **Map Controls**:
  - Zoom widgets and layer switcher controls must be styled to match the dark/light card surfaces (`var(--bg-surface)`), with rounded corners (`var(--radius-md)`) and refined borders.
- **Bounding Box Overlay**:
  - Displays the official Romanian TransDatRO grid boundary (`[43.5, 20.2]` to `[48.3, 29.8]`).
  - Subtle styling: dashed border with low opacity fill (10%) to preserve basemap legibility.
- **Geodetic Markers**:
  - Single point: Pulse pin marker in primary accent with white center dot. Out-of-grid points display in `--color-danger`.
  - Batch clusters: Styled with `.marker-cluster` using primary tech blue and monospace count numbers.

---

## 5. Checklist for Future UI Changes

Before committing any new UI component or modifying existing views, verify:
- [ ] Are tabs and buttons free of informal emoji decorations?
- [ ] Do segmented controls have vertical dividers and use neutral elevated active states?
- [ ] Are all coordinates and numerical values rendered in monospace typography?
- [ ] Are green and red reserved strictly for system validation/health and never used for geodetic data measurements?
- [ ] Are border radii consistent with the scale (`8px` for inputs/buttons/tabs, `12px` for cards)?
- [ ] Does the component adapt cleanly to both Dark (`data-theme="dark"`) and Light (`data-theme="light"`) modes?
