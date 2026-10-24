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
├── index.html              # HTML entrypoint with preloaded fonts and layout shells
├── package.json            # Node dependencies and build scripts
├── playwright.config.ts    # Chromium Playwright E2E test runner configuration & webServer orchestration
├── tsconfig.json           # Strict TypeScript configuration
├── vite.config.ts          # Vite build options and dev server API proxy
├── e2e/                    # Playwright end-to-end browser test suite (70 tests across 12 specs)
│   ├── fixtures.ts         # Offline tile & font stubs and page helpers
│   ├── helpers/            # Shared locator and assertion helpers
│   ├── smoke.spec.ts       # Critical DOM shell and initial layout smoke tests
│   ├── point.spec.ts       # Single-point transformation, direction toggle, angle units, map popups
│   ├── batch.spec.ts       # Batch text parsing, comments, out-of-grid suppression, 1-to-1 lines
│   ├── file.spec.ts        # File drag-and-drop, queue, streaming progress, export download
│   ├── map.spec.ts         # Map zoom constraints, layer toggling, centering vs fit-bounds
│   ├── i18n.spec.ts        # Dynamic RO/EN language switching, UI translation reactivity
│   ├── errors.spec.ts      # Syntax errors, out-of-grid badges, API rate limits
│   ├── computation.spec.ts # Computation timeline dashboard, KPI cards, Active Period deck, UK date format
│   ├── spatial.spec.ts     # Spatial density heatmap matrix, quadrant distribution, Active Period deck, UK date format
│   ├── lifecycle.spec.ts   # Client SPA lifecycle, modal management, and deep links
│   ├── legacy.spec.ts      # Legacy TransdatOnline historical telemetry subpage
│   └── referencePoints.ts  # Reference coordinate fixtures across Romania
└── src/
    ├── main.ts             # Application orchestrator & component lifecycle manager
    ├── assets/             # Vector graphic assets and dual-theme brand icon
    │   ├── logo.svg        # Official Map Projection Graticule standalone logo
    │   └── logos.ts        # Dual-Theme Map Projection Graticule SVG generator & favicon synchronizer
    ├── components/         # Interactive UI workbench components
    │   ├── Header.ts       # Header (Brand, 2-Line Attribution, Help & Language Dropdown Menus, Theme Toggle)
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
    └── en.ts               # English translation dictionary
    ├── style/              # Design system & styles (Pure Vanilla Modular CSS)
    │   ├── variables.css   # GitHub Primer Dark/Light palette & Subtle Surface + Border (#8) primary button tokens
    │   ├── base.css        # CSS reset, body defaults, custom scrollbars
    │   ├── layout.css      # Workbench layout (480px), navbar, tabs bar, studio cards, modals, toasts
    │   ├── forms.css       # Inputs, selects, sliding toggles, segmented buttons, dropzone
    │   ├── buttons.css     # Primary CTA buttons, icons, copy buttons, map actions
    │   ├── badges.css      # Status pills, Quasigeoid height metric field, quality report cards
    │   ├── menus.css       # Top-navigation menu triggers, dropdown menus, active checkmarks, dividers
    │   ├── panels.css      # Draggable & resizable floating windows (Grid Metadata, Help, Feedback)
    │   ├── telemetry.css   # KPI cards, timeline charts, interactive legend color popovers, spatial matrix deck
    │   └── map.css         # Theme-synchronized Leaflet widgets, popups & cluster styling
    ├── ui/                 # Reusable TypeScript UI factory controls
    │   ├── menu.ts         # renderDropdownMenu with nav-menu-trigger, active checkmarks, and full ARIA
    │   ├── segmentedControl.ts # renderSegmentedControl with arrow-key keyboard navigation
    │   ├── selectControl.ts # renderSelect for uniform dropdown lists
    │   ├── copyBadge.ts    # renderCopyBadge for interactive copyable chips
    │   ├── copyButton.ts   # renderCopyButton with animated checkmark feedback
    │   ├── index.ts        # Unified re-export entrypoint
    │   └── ui.test.ts      # Unit tests for UI factory helpers & dual-theme brand logo
    └── utils/              # Client-side geodetic & formatting helpers
        ├── apiClient.ts    # Type-safe REST client (passes X-Client-Source: web_app)
        ├── dmsParser.ts    # Unspaced DMS <-> Decimal Degrees converter
        ├── pointParser.ts  # Auto-delimiter detector & strict coordinate validator
        ├── exporters.ts    # Client-side CSV, GeoJSON, and AutoCAD DXF generator
        ├── gridBounds.ts   # Coarse SPG grid envelope union validation
        ├── numberFormat.ts # Locale-aware number & compact formatter for telemetry (RO comma, EN dot)
        └── toast.ts        # Non-intrusive floating toast notifications
```

---

## 3. Key Architecture & Design Invariants

1. **Decoupled REST API Integration**:
   - The frontend communicates with the transformation engine exclusively over HTTP via `src/utils/apiClient.ts`.
   - All requests send the HTTP header `X-Client-Source: web_app`, allowing the backend `SqliteUsageLogger` to attribute telemetry appropriately.

2. **GitHub Primer Dual-Theme Palette & Subtle Surface + Border (`#8`) Primary Button**:
   - Built with pure CSS custom properties (`src/style/variables.css`) parameterized by `data-theme="dark|light"`.
   - Detailed design specifications and the 27 Component Archetypes across 9 functional categories are documented in [`docs/webapp/design_system.md`](../docs/webapp/design_system.md) and [`docs/webapp/archetypes.md`](../docs/webapp/archetypes.md).
   - **GitHub Primer Neutral Palette**:
     - **Dark Theme (`data-theme="dark"`)**: `--bg-app: #0d1117`, `--bg-surface: #161b22`, `--bg-surface-elevated: #21262d`, `--border-subtle: #21262d`, `--text-main: #e6edf3`, `--btn-primary-bg: #1f2630`, `--btn-primary-border: #48586a`.
     - **Light Theme (`data-theme="light"`)**: `--bg-app: #f6f8fa`, `--bg-surface: #ffffff`, `--bg-surface-elevated: #f3f4f6`, `--border-subtle: #d0d7de`, `--text-main: #1f2328`, `--btn-primary-bg: #eaeef2`, `--btn-primary-border: #8c959f`.
   - **Official Dual-Theme Brand Logo (`Map Projection Graticule`, `src/assets/logos.ts`)**: Reactive CSS-variable inline SVG (`getLogoInlineSvg(36)`) and automatic favicon theme synchronization (`syncFaviconWithTheme`).
   - Fixed `480px` left workbench panel width with seamless folder bookmark tabs physically merging with the card container below (zero gap).
   - Dedicated Quasigeoid Height box ($\zeta = h - H$) styled in the technical metric palette (`--color-metric-text`), strictly isolated from status colors.

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

To verify code quality, test suite passes, and build integrity:

```bash
# Run automated utility and UI unit tests (Vitest)
npm test

# Run unit tests in interactive watch mode during development
npm run test:watch

# Run Playwright End-to-End (E2E) browser tests in headless mode
npm run test:e2e
# Alternatively from repo root via PowerShell runner:
.\test-e2e.ps1

# Run Playwright E2E tests in visual interactive UI runner
npm run test:e2e:ui
# Alternatively from repo root via PowerShell runner:
.\test-e2e.ps1 -Ui

# Type check TypeScript and verify production bundle compiles cleanly
npm run build
```

### Unit Tests (Vitest)
The unit test suite executes 80 tests in $< 400\text{ ms}$ covering:
- **`dmsParser.test.ts`**: Unspaced DMS parsing (`DD°MM'SS.SSSSS"`), whitespace handling, degree symbol variations, sexagesimal conversion.
- **`pointParser.test.ts`**: Delimiter detection (tabs, spaces, commas, semicolons), point IDs, coordinate ordering (`NE(H)` vs `EN(H)`).
- **`numberFormat.test.ts`**: Romanian comma `,` decimal formatting across telemetry and UI display surfaces while strictly preserving dot-decimal notation for CAD/GIS coordinate outputs.
- **`gridBounds.test.ts`**: Romania envelope pre-flight bounding box validation.
- **`exporters.test.ts`**: Pure client-side CSV, GeoJSON, and AutoCAD DXF R12 generator output fidelity.
- **`ui.test.ts`**: Accessible UI factory helper functions (`renderDropdownMenu` with active checkmarks, `renderSegmentedControl`, etc.) and the dual-theme `Map Projection Graticule` brand logo (`getLogoInlineSvg`, `getLogoStandaloneSvg`).

### End-to-End Browser Tests (Playwright)
The Playwright suite executes 70 browser specs across 12 test suites against an isolated FastAPI backend (`http://127.0.0.1:8100`) and compiled bundle (`dist-e2e/`), with offline network stubs for map tiles and fonts:
- **`smoke.spec.ts`**: DOM shell, layout containers, and initial component mounting.
- **`point.spec.ts`**: Single-point coordinate conversions, direction toggling, angular units, copy buttons, map pin popups, untransformable map click coordinate erasing, and layer color preservation with pulse animation.
- **`batch.spec.ts`**: Delimiter detection, coordinate swapping, comment/blank line preservation, error row formatting (`PointId, # [reason]`), textarea vertical height persistence, scrollbar hover cursor, and export triggers.
- **`file.spec.ts`**: File dropzone, queue management, streaming transformation progress, and batch downloads without premature toasts.
- **`map.spec.ts`**: Zoom constraints, Home button, layer visibility toggles, measurement ruler, and center vs fit-bounds actions.
- **`i18n.spec.ts`**: Reactive Romanian (`RO`) / English (`EN`) dropdown switching, attribute updates, and localStorage persistence.
- **`errors.spec.ts`**: Pre-flight validation, Two-Section Quality Report card layout (`Coordinates` top, dashed divider, `Transformation` bottom), impossible angle detection, and rate-limit error displays.
- **`legacy.spec.ts`**: Legacy TransdatOnline historical telemetry subpage, SVG curve rendering, timeline series toggles, and date presets.
- **`lifecycle.spec.ts`**: Brand Logo fresh state reset, soft title navigation, offline alert banner & health dot reactivity, floating panels (Grid Metadata search filter, Geodetic Help live benchmark card), and bot challenge verification.
- **`computation.spec.ts`**: Usage Computation subpage, centered Archetype 8.1 KPI cards, calendar date range clamping, and segmented source curves.
- **`spatial.spec.ts`**: Usage Spatial subpage, quadrant breakdown cards, and 0.1° punchcard matrix.



---

## 9. Single-Source-of-Truth Versioning & Build Injection

- The frontend package version is declared in `webapp/package.json` and synchronized across repository metadata via `tools/bump_version.py`.
- At compile time, Vite injects `__APP_VERSION__` via `define` in `vite.config.ts`.
- Exposed centrally in `src/config.ts` (`APP_CONFIG.appVersion`) and rendered reactively across the brand header badge (`v${APP_CONFIG.appVersion}`) and telemetry/feedback reporting payloads.
