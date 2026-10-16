# PyTransDatRO Web Application Architecture & Specification

## 1. Overview & Architecture

The `webapp/` package provides a modern, high-performance Web GIS Single Page Application (SPA) designed to bring Romania's official geodetic coordinate transformation platform (**TransDatRO**) into modern web browsers.

It replaces the legacy Java/GWT TransDatOnline web interface with a responsive, zero-framework, client-side application built with **Vite**, **TypeScript**, **Vanilla CSS**, and **Leaflet**.

```mermaid
flowchart TD
    subgraph Browser["Web Browser Client (webapp/)"]
        SharedConfig["Shared Config Bar<br/>(Direction, Angular Unit, Coordinate Order, Format Badge)"]
        UI["Workbench Interface<br/>(Point / Text / Files)"]
        HeaderNav["Header Navigation<br/>(Grid Info, Geodetic Help, Feedback, Theme/Lang)"]
        State["i18n & Theme State<br/>(RO / EN, Dark / Light / System)"]
        Parsers["Client Parsers & Exporters<br/>(Unspaced DMS, SPG Grid Bounds, DXF / CSV / GeoJSON)"]
        Map["Leaflet GIS Map<br/>(3 Layers: Point/Text/Files, Romanian Bounds Guard, Dual Nav)"]
        Client["Type-Safe API Client<br/>(X-Client-Source: web_app)"]
        
        HeaderNav --> State
        SharedConfig --> UI
        UI --> State
        UI --> Parsers
        UI --> Map
        UI --> Client
    end

    subgraph Backend["FastAPI Web Server (api/)"]
        Static["StaticFiles Mount<br/>GET / -> webapp/dist/"]
        APIRoutes["Modern REST API<br/>/api/v1/transform/*<br/>/api/v1/grid/info<br/>/api/v1/telemetry/stats"]
    end

    subgraph CoreEngine["Pure Python Core (pytransdatro/)"]
        TransRO["TransRO RAM Singleton<br/>(Bicubic Spline + Helmert2D)"]
        Logger["SqliteUsageLogger<br/>(logs/usage.db)"]
    end

    Client -->|HTTP / JSON / Chunked Stream| APIRoutes
    Browser -.->|Initial Load| Static
    APIRoutes --> TransRO
    TransRO -.-> Logger
```

---

## 2. Technology Stack & Key Invariants

1. **Vanilla CSS Design System & Modular Stylesheets**:
   - Zero bulky CSS frameworks (no Bootstrap, no Tailwind runtime overhead).
   - Driven entirely by custom CSS design tokens defined in [`webapp/src/style/variables.css`](../webapp/src/style/variables.css) and documented in [`webapp/design_guide.md`](../webapp/design_guide.md).
   - **Modular Stylesheet Architecture (`webapp/src/style/`)**:
     - `layout.css`: 480px workbench panel, top navigation bar, folder bookmark tabs, studio cards, modals, and toasts.
     - `forms.css`: Text inputs, selects, sliding direction toggle, segmented switches, and dropzone.
     - `buttons.css`: High-contrast theme-native CTA (`.btn-primary`), secondary, outline, copy buttons, and action buttons.
     - `badges.css`: Status pills, Quasigeoid height metric field, angle format badges, and quality report cards.
     - `menus.css`: Top navigation menu item triggers (`.nav-menu-trigger`), dropdown panels, menu items, and dividers.
     - `telemetry.css`: KPI cards, SVG timeline chart styling, and spatial matrix canvas deck.
   - **The Trinity Balance**: Simplicity (zero decorative emojis), Optimal ergonomics (480px workbench panel, asymmetrical 3-column coordinate layout), and Modern engineering precision.
   - **Theme-Native High-Contrast CTA**: Porcelain `#f8fafc` on dark obsidian; slate `#0f172a` on light porcelain.
   - **Metric Color Isolation**: Physical scalar measurements (e.g., Quasigeoid Height $\zeta = h - H$) use dedicated technical cyan (`--color-metric-text: #38bdf8;`), strictly isolated from status colors (never green or red).
   - **Seamless Folder Bookmark Tabs**: Physical continuity between the `.tabs-bar` and the active card container below with zero background gap.

2. **Reusable TypeScript UI Factory System (`webapp/src/ui/`)**:
   - `renderDropdownMenu(props)`: Generates accessible dropdown menus with `.nav-menu-trigger` menu-item aesthetics, rotating chevron animations, ARIA semantics, and keyboard navigation (`Escape`, `ArrowDown`, `ArrowUp`).
   - `renderSegmentedControl(props)`: Accessible radio-switch controls with keyboard arrow navigation.
   - `renderSelect(props)`: Type-safe dropdown selects with focus ring management.
   - `renderCopyBadge(props)`: Interactive copyable text badges with clipboard integration.
   - `renderCopyButton(props)`: Standard code-block copy buttons with animated checkmark feedback.

3. **Strict TypeScript & Vite**:
   - Strict TypeScript compilation (`tsconfig.json`) ensuring robust domain typing for coordinates, operations, units, and API payloads.
   - Bundled with Vite for instant Hot Module Replacement (HMR) during development and optimized tree-shaken production bundles (~105 KB gzipped).

4. **Leaflet Web GIS Engine & 3-Layer Architecture**:
   - High-performance, lightweight mapping powered by `leaflet` and `leaflet.markercluster`.
   - **3 Independent LayerGroups**: `pointLayer`, `textLayer`, and `filesLayer`, independently managed with live point counters and 6-color customizable swatch pickers (Cyan, Cobalt, Emerald, Amber, Violet, Rose) with `localStorage` persistence.
   - **Full Marker & Cluster Color Synchronization**: Dynamic cluster icons (`custom-cluster-circle`) inherit the exact active swatch color, synchronizing markers, zoom-out cluster bubbles, Layers drawer swatches, and Legend dots.
   - **Decoupled Viewports & Dual Navigation**: Coordinate transformations never auto-pan or jump the map; dedicated `Center on Map` and `Zoom to Extents` actions are provided across Point, Text, and Files tabs.
   - **Boundary-Safe Guard**: Map centering and bounds fitting are strictly guarded within Romania's geographic grid envelope (`[43.5, 20.2]` to `[48.3, 29.8]`), preventing map jumps on invalid coordinates.
   - **Custom Ergonomic Docks**:
     - *Top-Left Controls Dock*: Grouped buttons for Zoom In/Out, **Zoom to Romania (Home)**, **Metric Distance Measurement Ruler**, and **Fullscreen Toggle** (with `fullscreenchange` synchronization).
     - *Top-Right Expandable Legend & Scale Bar*: Stateful persistent Legend widget (persisted in `localStorage['pytransdat_legend_open']`) that automatically reflows above a high-contrast segmented **GIS Graphic Scale Bar** (with start `0` and end `[max] m/km` labels).
     - *Bottom-Left Layers Drawer*: 32x32px icon-only button opening a structured drawer organized into Data Layers, Shift Grids (Coming Soon), and Base Maps (OpenStreetMap default, with user choice persisted in `localStorage['pytransdat_basemap']`).
     - *Bottom-Center Coordinates Badge*: Abbreviated `Lat: [val]  Long: [val]` responsive to active angular units (Degrees, DMS, Radians), elevated on mobile screens.
   - **Uniform Geodetic Popups**: Unified attribute tables across all layers distinguishing Stereo 70 Normal Altitude ($H$, Black Sea 1975) from ETRS89 Ellipsoidal Height ($h$, GRS80), Quasigeoid height $\zeta$ in cyan, and 1-click **Copy Attributes** (`Copiază atributele`) key-value clipboard action with in-place feedback.

4. **Decoupled API Contract**:
   - The frontend communicates with the backend exclusively via standard REST endpoints (`/api/v1/...`).
   - Every outbound request injects the HTTP header `X-Client-Source: web_app`, allowing backend telemetry to accurately differentiate browser traffic from automated scripts or desktop CAD plugins.

---

## 3. Core Functional Modules

### A. 3-Tier Shared Transformation Configuration Bar (`main.ts`)
Elevated above all tabs, governing `Point`, `Text`, and `Files` globally:
- **Tier 1 (Direction Toggle)**: Sliding directional toggle with a sleek 6px rounded-rectangle glider (`var(--radius-sm)`). Default: `ETRS89 ➔ Stereo 70` (first button); `Stereo 70 ➔ ETRS89` (second button).
- **Tier 2 (50/50 Balanced Split)**: `Angular unit: [select]` (Degrees, DMS, Radians) | `Coordinate order: [ NE(H) | EN(H) ]`. Sentence-case labeling.
- **Tier 3 (Dynamic Format Example)**: Selectable badge displaying copy-ready formatted coordinates (`45.84123450`, `45°50'28.44420"`, `0.7997813210`). Clicking the badge copies strictly the numeric angle value without label text.

### B. Single Point Transformation (`PointConverter.ts`)
- **Projection Center Defaults**: Initial input coordinates set to Romania's projection center origin:
  - Latitude $\phi = 46.0^\circ\text{ N}$ (`46°00'00.00000"` in DMS)
  - Longitude $\lambda = 25.0^\circ\text{ E}$ (`25°00'00.00000"` in DMS)
  - Ellipsoidal Height $h = 100.0\text{ m}$
- **Asymmetrical 3-Column Layout**: Northing/Lat (`1.35fr`), Easting/Lon (`1.35fr`), and Altitude/Height (`0.8fr`) so long DMS strings never wrap or clip.
- **Auto 2D/3D Dimensionality Detection**: Empty 3rd coordinate $\to$ 2D mode; provided 3rd coordinate $\to$ 3D mode (calculates $\zeta = h - H$).
- **Output Hygiene**: Pure numeric values only (strips `m`, `°`, `rad` from text boxes). Displays persistent `—` when cleared or invalid, preventing vertical layout shift.
- **Top-Right Code-Block Copy Action**: Header copy icon with instant `✓ Copied` feedback.
- **Focus Auto-Select**: All coordinate inputs execute `input.select()` on focus.
- **Dedicated Quasigeoid Height Field**: Rendered in its own dedicated box on its own row using technical metric cyan-blue (`--color-metric-text`).
- **Unspaced DMS Formatting & Paste Guard**: `DD°MM'SS.SSSSS"` formatting with single-DMS paste splitting protection.
- **Client-Side SPG Pre-Flight Check**: Instant validation against Romania's grid envelope ([`gridBounds.ts`](../webapp/src/utils/gridBounds.ts)) before dispatching network calls.
- **Unified Quality Report**: Positioned directly below outputs, displaying `In Grid`, `Outside Grid`, and calculation mode.

### C. Batch Text Transformation (`BatchConverter.ts`)
- **Strict Vertical Layout**: Input textarea $\to$ Primary Transform button $\to$ Output textarea $\to$ Quality Report $\to$ Download actions $\to$ Map actions.
- **Strict Numeric Validation**: Regex check (`isValidDecimal`) prevents non-numeric typos (e.g. `5001a50.000`) from being partially evaluated into false coordinates.
- **Strict 1-to-1 Line Correspondence**: Input row $N$ maps directly to output row $N$. Empty lines, comments, and invalid rows maintain their exact sequence index.
- **Human 1-Based Indexing**: Unnamed coordinates are referenced as `Point 1`, `Point 2`, etc.
- **Out-of-Grid Suppression**: Suppresses bogus approximate transformed coordinates, emitting only `[Point_ID], # Out of grid`.
- **Top-Right Copy Action**: Output header includes a dedicated `Copy Results` button.
- **Relocated Download Group**: Export actions (`CSV`, `GeoJSON`, `DXF`) grouped below output under `Download As:`.

### D. Bulk File Streaming Pipeline (`FileUpload.ts`)
- **Interactive Queue**: Multi-file dropzone with file list, file size indicators, individual `✕` remove buttons, and manual `Transform Files` trigger.
- **Raw Preview Elimination**: Replaced 5-line text preview with clean queue metadata and post-transformation Quality Report cards.
- **Line Streaming Robustness**: Auto-detects Point ID column (3 vs 4 tokens), parses DMS angles with sexagesimal radian conversion, outputs unspaced DMS, and handles Windows-1252 / Latin-1 degree symbol fallback decoding.
- **Downloads**: Individual transformed file downloads plus batch `Download All (ZIP)`.

### E. Interactive Web GIS Map (`MapView.ts`)
- **Desktop Split-Screen**: Left-hand 480px workbench panel paired with a full-height interactive Leaflet map.
- **3-Layer Independent Architecture**: `pointLayer`, `textLayer`, and `filesLayer` allow users to toggle visibility independently, customize colors via 6-swatch popovers, and inspect live point counters.
- **Full Cluster Synchronization**: Cluster grouping circles (`custom-cluster-circle`) dynamically inherit their layer's active swatch color, keeping individual markers, zoom clusters, drawer swatches, and legend dots in perfect visual synchronization.
- **Decoupled Viewports & Dual Navigation**:
  - `Center on Map`: Pans to centroid of valid in-grid coordinates while preserving user zoom level.
  - `Zoom to Extents`: Pans and fits bounds to encompass all points (or zooms to level 14 for a single point).
- **Smart Map Click**:
  - In `ETRS89ToStereo70`: Populates clicked lat/lon formatted in active angular units.
  - In `Stereo 70ToETRS89`: Calls backend API to reverse-calculate clicked point into Stereo 70 projection meters. Map clicks only capture coordinates when the Single Point tab is active.
- **Ergonomic Map Controls**:
  - *Top-Left Dock*: Zoom In/Out, **Zoom to Romania (Home)** extents button, **Single-Distance Measurement Tool** (with clear action), and **Fullscreen Toggle** (synchronized via `fullscreenchange`).
  - *Top-Right Dock*: Stateful persistent Legend widget (persisted in `localStorage['pytransdat_legend_open']`) positioned above a high-contrast segmented **GIS Graphic Scale Bar** with start `0` and end `[max] m/km` labels.
  - *Bottom-Left Drawer*: 32x32px icon-only button opening a structured drawer containing Data Layers, Shift Grids (Coming Soon), and Base Maps (persisting user choice to `localStorage['pytransdat_basemap']`).
  - *Bottom-Center Badge*: Abbreviated coordinates readout (`Lat: [val]  Long: [val]`) responsive to active angular units, elevated on mobile viewports.
- **Uniform Geodetic Popups**:
  - Popups across Point, Text, and Files layers share a structured attribute table displaying: Layer Name, `Point Name/Number`, source `File` (for bulk files), Northing, Easting, Normal Altitude ($H$, Black Sea 1975), Latitude, Longitude, Ellipsoidal Height ($h$, GRS80), Quasigeoid height $\zeta$ in cyan, and bottom `Grid Coverage` status.
  - **Copy Attributes Action**: 1-click clipboard button copying all point attributes in structured key-value pairs (`Key: Value`) with temporary `✓ Copied` feedback.
- **High-Volume Guardrail Modal (`VolumeWarningModal.ts`)**:
  - Uncapped file ingestion without row truncation. Datasets exceeding **2,000 points** prompt the user before plotting:
    1. *Render All Points* (with Marker Clustering).
    2. *Sample Preview* (render first 1,000 points only).
    3. *Skip Map* (display table/export results only to preserve responsiveness).
  - Clearing or removing files in the upload queue automatically invokes `mapView.clearFilePoints()` to keep map layers synchronized.

### F. Client-Side Multi-Format Exporters (`exporters.ts`)
- **CSV**: Configurable delimiter (Comma, Semicolon, Tab) with clean header metadata, respecting active coordinate order (`NE` vs `EN`) and suppressing coordinates for out-of-grid points.
- **GeoJSON**: Standard WGS84/ETRS89 `Point` FeatureCollections with coordinate properties.
- **AutoCAD DXF**: R12 ASCII DXF generating native `POINT` entities for immediate drag-and-drop import into AutoCAD and Civil 3D.

### G. Header Navigation & Documentation Modals (`Header.ts`)
- **Grid Info (`StatsView.ts`)**: Operational telemetry, daily point volume charts, and active `.spg` grid metadata.
- **Geodetic Help (`GeodeticHelpModal.ts`)**: Interactive reference guide documenting Romanian datums (Stereo 70, ETRS89, MN75, GRS80), Quasigeoid undulation ($\zeta = h - H$), and syntax conventions.
- **Feedback Dialog (`FeedbackModal.ts`)**: Discrepancy reporting against desktop TransDat or national GNSS stations.

---

## 4. Internationalization (`i18n`) & Aesthetics

The application is fully bilingual:
- **Romanian (`ro`)**: Default system language, matching official Romanian geodetic terminology (e.g. *Cota elipsoidală*, *Stereografic 1970*, *Înălțimea cvasigeoidului*).
- **English (`en`)**: Complete localization for international users and cross-border projects.

### Reactive State Flow:
```
User clicks RO/EN Switcher
   │
   ▼
i18nStore.setLanguage('en')
   │
   ├── Persists to localStorage['pytransdat_lang']
   └── Triggers all registered subscriber callbacks (re-renders labels in-place)
```

---

## 5. Deployment & Production Serving

In production, the compiled frontend bundle is hosted directly by the FastAPI web service:

```
FastAPI Server (api/main.py)
   ├── /api/v1/*                --> JSON & Streaming Endpoints
   ├── /transdatonline/*        --> Legacy Compatibility Endpoints
   ├── /health, /docs, /redoc   --> Observability & OpenAPI
   └── /                        --> StaticFiles(directory="webapp/dist", html=True)
```

### Build Instructions:

```bash
# 1. Compile frontend assets
cd webapp
npm install
npm run build

# 2. Start the unified web service
cd ..
D:\anaconda3\envs\pytransdat\python.exe -m uvicorn api.main:app --host 0.0.0.0 --port 8000
```

The web app is then immediately served at `http://127.0.0.1:8000/`.

---

## 6. Verification & Quality Assurance

- **Frontend Build Verification**:
  ```bash
  npm --prefix webapp run build
  ```
  Ensures 0 TypeScript type errors and successful Vite bundling into `webapp/dist/`.
- **Backend Test Suite**:
  ```bash
  D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/
  ```
  Ensures all 76 core, API, streaming, and telemetry tests pass without regressions.
