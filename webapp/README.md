# PyTransDatRO Web Application (`webapp/`)

Modern Web GIS Single Page Application (SPA) serving as the browser interface for **PyTransDatRO** (Romania's official geodetic coordinate transformation platform).

Built with **Vite**, **TypeScript**, **Vanilla CSS**, and **Leaflet**.

---

## 1. Quick Start

### Prerequisites
- **Node.js**: v18+ (LTS recommended)
- **FastAPI Backend**: Running at `http://127.0.0.1:8000` (see repo root instructions)

### Development Workflow

```bash
# Navigate to the frontend directory
cd webapp

# Install dependencies
npm install

# Start development server (with hot module replacement)
npm run dev

# Run automated utility unit tests (Vitest)
npm test

# Run tests in interactive watch mode
npm run test:watch
```

The application will be accessible at:
👉 **`http://localhost:5173/`**

> **Note on API Proxying**: Vite is configured (`vite.config.ts`) to transparently proxy all `/api`, `/health`, and `/transdatonline` calls to `http://127.0.0.1:8000`. Keep the Python backend running during frontend development.

### Production Build & Integration

```bash
# Type-check TypeScript and build production bundle
npm run build
```

This compiles optimized, minified assets into `webapp/dist/`.

When the FastAPI server starts (`python -m uvicorn api.main:app`), it automatically detects `webapp/dist/` and mounts it on the root URL (`http://127.0.0.1:8000/`) using FastAPI's `StaticFiles(html=True)`.

---

## 2. Directory Structure

```
webapp/
├── design_guide.md         # Single source of truth design system & styling manual
├── index.html              # HTML entrypoint with preloaded fonts and layout shells
├── package.json            # Node dependencies and build scripts
├── tsconfig.json           # Strict TypeScript configuration
├── vite.config.ts          # Vite build options and dev server API proxy
└── src/
    ├── main.ts             # Application orchestrator & component lifecycle manager
    ├── assets/             # Vector graphic assets and brand icons
    │   └── logo.svg        # Official Geodetic Shift Lattice application logo
    ├── components/         # Interactive UI workbench components
    │   ├── Header.ts       # 3-Zone Header (Brand, Navigation, Native Details Help Dropdown, Health)
    │   ├── GridMetadataPanel.ts # Non-modal, draggable & resizable SPG grid metadata tree panel
    │   ├── GeodeticHelpPanel.ts # Non-modal, draggable & resizable geodetic reference & benchmark panel
    │   ├── PointConverter.ts  # Single-point transformation workbench
    │   ├── BatchConverter.ts  # Multi-line text conversion workbench
    │   ├── FileUpload.ts      # Drag-and-drop file processor with streaming progress
    │   ├── MapView.ts         # Leaflet map, 3 independent layer groups, dual map navigation
    │   ├── StatsView.ts       # Active grid metadata and operational status modal
    │   ├── TelemetryDashboard.ts # Full-page usage telemetry dashboard (Computation & Spatial views)
    │   ├── TimelineChart.ts   # Interactive SVG dual-metric timeline chart with dynamic source curves
    │   ├── SpatialMatrixChart.ts # 0.1° punchcard canvas matrix with quadrant lines and ghost cells
    │   ├── VolumeWarningModal.ts # Guardrail prompt for transforms > 2,000 points
    │   ├── ReportModal.ts     # Diagnostic calculation details modal
    │   └── FeedbackPanel.ts   # Non-modal draggable, resizable feedback & discrepancy panel
    ├── i18n/               # Internationalization system
    │   ├── index.ts        # Reactive translation store & persistent state
    │   ├── ro.ts           # Romanian translation dictionary (default)
    │   └── en.ts           # English translation dictionary
    ├── style/              # Design system & styles (Pure Vanilla CSS)
    │   ├── variables.css   # Color palette tokens, dark/light themes, typography
    │   ├── base.css        # CSS reset, 480px workbench panel layout
    │   ├── components.css  # Buttons, pills, form inputs, folder tabs, modals, toasts
    │   └── map.css         # Leaflet container, custom popups & cluster styling
    └── utils/              # Client-side geodetic & formatting helpers
        ├── apiClient.ts    # Type-safe REST client (passes X-Client-Source: web_app)
        ├── dmsParser.ts    # Unspaced DMS <-> Decimal Degrees converter
        ├── pointParser.ts  # Auto-delimiter detector & strict coordinate validator
        ├── exporters.ts    # Client-side CSV, GeoJSON, and AutoCAD DXF generator
        ├── gridBounds.ts   # Coarse SPG grid envelope union validation
        └── toast.ts        # Non-intrusive floating toast notifications
```

---

## 3. Key Architecture & Design Invariants

1. **Decoupled REST API Integration**:
   - The frontend communicates with the transformation engine exclusively over HTTP via `src/utils/apiClient.ts`.
   - All requests send the HTTP header `X-Client-Source: web_app`, allowing the backend `SqliteUsageLogger` to attribute telemetry appropriately.

2. **Vanilla CSS Design System & Visual Ergonomics**:
   - Built with pure CSS custom properties (`src/style/variables.css`).
   - Detailed specifications, visual tokens, and component guidelines are documented in [design_guide.md](design_guide.md).
   - Supports **Dark Mode**, **Light Mode**, and **System Preference** with no heavy UI framework dependencies.
   - Clean, modern geodetic aesthetic: Deep obsidian/slate neutrals, high-contrast theme-native anchor CTA buttons, monospace coordinate typography (`JetBrains Mono`).
   - Fixed `480px` left workbench panel width with seamless folder bookmark tabs physically merging with the card container below (zero gap).
   - Dedicated Quasigeoid Height box ($\zeta = h - H$) styled in technical metric cyan (`--color-metric-text: #38bdf8;`), strictly isolated from status colors.

3. **3-Tier Shared Transformation Configuration**:
   - Elevated above tabs in `main.ts` to globally manage `Point`, `Text`, and `Files` tabs:
     - **Tier 1**: Direction sliding toggle (`ETRS89 ➔ Stereo 70` default / `Stereo 70 ➔ ETRS89`).
     - **Tier 2**: 50/50 split between `Angular unit:` (Degrees, DMS, Radians) and `Coordinate order:` (`NE(H)` vs `EN(H)`).
     - **Tier 3**: Dynamic, selectable angle format example badge with click-to-copy functionality.

4. **Bilingual Internationalization (`i18n`)**:
   - Romanian (`ro`) is the primary default language; English (`en`) is fully supported.
   - User language preference is stored in browser `localStorage`.
   - UI re-renders reactively upon language change without full page reloads.

5. **Map Performance, Controls & Geodetic Ergonomics**:
   - Built with **Leaflet** and **Leaflet.markercluster**.
   - **3 Independent LayerGroups**: `pointLayer`, `textLayer`, and `filesLayer` with live point counts and 6-swatch customizable color pickers (Cyan, Cobalt, Emerald, Amber, Violet, Rose) with `localStorage` persistence.
   - **Full Marker & Cluster Color Synchronization**: Dynamic cluster icons (`custom-cluster-circle`) inherit the active swatch color, synchronizing individual dots, zoom-out cluster bubbles, Layers drawer swatches, and Legend dots.
   - **Ergonomic Dock Layout**:
     - **Top-Left Controls Dock**: Zoom In/Out, **Zoom to Romania (Home)** extents button, **Single-Distance Measurement Ruler** (`Distance: [value] km`), and **Fullscreen Toggle**.
     - **Top-Right Expandable Legend & GIS Graphic Scale**:
       - **Legend Widget**: Positioned top-right with an icon-only 32x32px toggle button and collapsible translucent panel displaying active layers with their synchronized color dots and live point counts. Operates as a persistent state dock that stays open during map navigation and automatically pushes the scale bar down without covering it.
       - **Custom Segmented GIS Graphic Scale Bar**: Positioned immediately beneath the Legend, featuring an alternating high-contrast monochrome segmented bar (pure white/obsidian in dark mode, solid black/white in light mode) with top division ticks (0, mid, max) and uncrowded start/end labels (`0` and `[max] m/km`).
     - **Bottom-Left Layers Drawer**: Dedicated 32x32px icon-only drawer button opening a structured panel organized as:
       1. **Data Layers**: Layer visibility checkboxes, point counts, and 6-color swatch pickers.
       2. **Shift Grids (Coming Soon)**: Preview cards for upcoming $\Delta E$, $\Delta N$, Total Vector, and $\zeta$ raster services.
       3. **Base Maps**: Radio list items matching data layer styling (OpenStreetMap default, CartoDB Dark Matter, CartoDB Positron, ESRI World Imagery) with user selection saved to `localStorage['pytransdat_basemap']`.
     - **Bottom-Center Coordinates Widget**: Unobtrusive live readout displaying compact abbreviations (`Lat: [val]  Long: [val]`) reactively formatted per the active `Angular Unit` (Degrees, DMS, Radians) without unit clutter or attribution overlap.
   - **Uniform Themed Geodetic Popups**: Custom Dark Obsidian and Light Porcelain cards presenting consistent attribute tables across Point, Text, and Files layers: Layer Name header, `Point Name/Number`, source `File` (for bulk files), unified coordinates distinguishing Stereo 70 Normal Altitude ($H$, Black Sea 1975) and ETRS89 Ellipsoidal Height ($h$, GRS80), Quasigeoid height in cyan, bottom `Grid Coverage` status, and robust 1-click "Copy Attributes" (`Copiază atributele`) clipboard action formatting key-value pairs with in-place `✓ Copied` feedback.
   - **CAD-Style Point Labels at High Zoom**: Subtle point ID labels rendered next to markers when zoomed in to level $\ge 15$.
   - **Decoupled Non-Disruptive Viewports**: Transformations never auto-pan; map clicks only populate coordinates when the Single Point tab is active.
   - **Large volume guardrail**: Prompts via `VolumeWarningModal` when transforms exceed 2,000 points.
   - Base maps: OpenStreetMap (default), CartoDB Dark Matter, CartoDB Positron, and ESRI World Imagery.
   - Min zoom reduced to 3 for full regional/global zoom range.

6. **Telemetry & Usage Analytics Architecture (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`)**:
   - Full-page analytics view accessible via the primary header's **Usage** (`Date Utilizare`) tab, featuring GitHub-style underline sub-navigation:
     - **Computation View (`TimelineChart.ts`)**:
       - Dynamic discovery of logging client sources (`rest_api`, `web_app`, `desktop_cli`) from `logs/usage.db` via `/api/v1/telemetry/timeline` with exact fallback for unmapped sources and exclusion of `unknown`.
       - Native themed `<input type="date">` calendar pickers with strict boundary clamping (`2023-11-01` to current date) and a range limits indicator pill.
       - Daily, Weekly, and Monthly temporal aggregation with preset intervals (`All Time`, `12M`, `YTD`, `30D`, `7D`).
       - Strictly additive Total curve ($\text{Total} \ge \text{each source}$) with high-visibility dashed styling.
       - Independent multi-metric filtering (`Points` and `Calls` toggleable individually or together) with neutral-toned Y-axes and a 2-line grouped symbology legend.
       - Luminous/pastel tints for Calls and darker tones for Points across client color palettes.
     - **Spatial View (`SpatialMatrixChart.ts`)**:
       - Macroscopic $0.1^\circ$ (~10 km) punchcard canvas matrix representing geographic distribution across Romania with subtle ghost cells for unvisited coordinates.
       - Mathematical quadrant summary cards (**NW**, **NE**, **SW**, **SE**) calculated around the national center ($45.9^\circ\text{ N}, 24.9^\circ\text{ E}$) with volume and percentage share.
       - Deduplicated full-width overview banner displaying Total volume and Active Cells count.
       - Cartographic coordinate notation (`°N` and `°E`) with the Latitude axis positioned on the left margin.
       - Discrete 5-step color ramp and Logarithmic / Linear scaling toggle in the canvas footer bar.
       - Synchronized time horizon presets, date clamping, and standalone source toggle buttons matching the Computation view.
     - **Theme-Native Visual Controls**:
       - Filter toggles use standalone buttons with neutral borders (`--border-medium`), elevated backgrounds (`--bg-surface-elevated`), 8px spacing, and high-contrast theme-native active states (`--btn-primary-bg` / `--btn-primary-text`), removing saturated electric blue accents.

7. **Client-Side Exporting**:
   - Results can be downloaded directly in the browser as:
     - **CSV** (Comma, Semicolon, or Tab delimited)
     - **GeoJSON** (`Point` features with coordinate and height properties)
     - **AutoCAD DXF** (R12 ASCII format with `POINT` entities for immediate CAD import)

---

## 4. Verification & Testing

To verify code quality and build integrity:

```bash
# Type check and verify production bundle compiles cleanly
npm run build
```
