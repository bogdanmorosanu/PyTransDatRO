# Changelog - PyTransDatRO Web Application (`webapp/`)

All notable changes to the PyTransDatRO web application and its interactive Web GIS interface will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
- **Legacy TransdatOnline Historical Statistics Subpage (TelemetryDashboard.ts, LegacyTimelineChart.ts)**:
  - Added new "Legacy" subpage under the "Usage" page presenting 15 years of TransdatOnline historical transformation volumes (2011 to 2026).
  - Built custom SVG LegacyTimelineChart supporting points-only visualization, smooth cubic bezier curve interpolation, gridlines, axis text, interactive crosshairs, and hover tooltips.
  - Interactive legend swatches with color picker popovers allowing live customization of line colors, persisted in localStorage.
  - Date presets (7d, 30d, ytd, 12m, ll) calculating dynamically relative to the archive boundary (2026-10-14).
  - Independent faceted toggles for client channels (**Web App**, **REST API**) and operations (**Stereo 70 ↔ ETRS89**, **Stereo 30 ↔ ETRS89**), plus an additive Total curve.
  - End-to-end Playwright browser test coverage in webapp/e2e/legacy.spec.ts.
- **FastAPI Endpoints for Legacy Telemetry (pi/routes/info.py, pi/schemas/info.py)**:
  - GET /api/v1/telemetry/legacy/stats: Overall summary of points transformed and timeline limits.
  - GET /api/v1/telemetry/legacy/timeline: Time-series data supporting daily, weekly, and monthly bucket expressions with date bounds slicing and faceted filtering.
- **Data Integrity & Deduplication Tooling (	ools/)**:
  - 	ools/audit_legacy_stats_integrity.py: Audits timeline continuity across 5,758 calendar days and outputs structured reports into eports/legacy_stats_integrity_report.txt.
  - 	ools/extract_legacy_stats.py: Resolves MySQL concurrency race conditions by aggregating records per calendar day and enforcing unique constraints in SQLite.

### Changed
- Refined archive timeline pill labels in English (Archive Earliest) and Romanian (Cea Mai Veche Dată) in webapp/src/i18n/.
- Updated Header.ts to dispatch a global 	hemechange event for instant chart theme responsiveness.

---

## [1.1.0] - 2026-10-18

### Added
- **Playwright End-to-End (E2E) Browser Test Suite (`webapp/e2e/`, `playwright.config.ts`, `test-e2e.ps1`)**:
  - Comprehensive Chromium E2E suite covering single-point coordinate transformations, batch text parsing, file streaming uploads/downloads, map centering and zoom constraints, bi-directional i18n (`EN` / `RO`) persistence, and API error/rate-limit display.
  - Implemented offline request stubbing in `e2e/fixtures.ts` for external map tiles and Google Fonts.
  - Added read-only `window.__e2e.mapView()` diagnostic snapshot in `MapView.ts` (active only in `VITE_E2E=true` builds).
- **GitHub Primer Neutral Dual-Theme Palette & Subtle Surface + Border (`#8`) Primary Button (`webapp/src/style/variables.css`)**:
  - Adopted GitHub Primer neutral Dark (`#0d1117` app canvas, `#161b22` card surface, `#21262d` elevated controls) and Light (`#f6f8fa` app canvas, `#ffffff` card surface, `#f3f4f6` elevated controls) tokens.
  - Established Option `#8` (**Subtle Surface + Border**) as the primary action button style (`#1f2630` / `#48586a` in Dark mode; `#eaeef2` / `#8c959f` in Light mode).
- **Official Dual-Theme Brand Logo (`Map Projection Graticule`, `webapp/src/assets/logos.ts`, `logo.svg`, `public/favicon.svg`)**:
  - Adopted **C2: Map Projection Graticule** as the official dual-theme vector logo, featuring curved stereographic graticule arcs intersecting a projection frame with a central origin node, rendered inline via CSS variables (`getLogoInlineSvg(36)`) and synchronized dynamically with the browser tab favicon (`syncFaviconWithTheme`).
- **Interactive Computation Chart Legend Color Picker & Dual-Theme Line Palette (`TimelineChart.ts`, `TelemetryDashboard.ts`, `telemetry.css`)**:
  - Replaced low-contrast tones with high-contrast Coral Rose / Crimson Rose for Desktop/CLI client sources.
  - Added click-to-change 10-swatch color popovers on every line item in the Computation chart legend, persisting custom line colors in `localStorage['pytransdat_timeline_colors']`.
- **Status-Displaying Language Selector Dropdown Menu (`Header.ts`, `ui/menu.ts`, `menus.css`)**:
  - Converted the header Language control into a disclosure Dropdown Menu (`#header-lang-details`, `renderDropdownMenu`) with a globe SVG icon, active language code (`RO` / `EN`), rotating chevron (`▾`), and `Română (RO)` / `English (EN)` menu items (`role="menuitemradio"`, `aria-checked`) with an active checkmark (`✓`, `.menu-check-icon`).
- **Dedicated Floating Windows Stylesheet (`webapp/src/style/panels.css`)**:
  - Restored and unified full CSS rules for `.floating-grid-panel`, `.floating-help-panel`, and `.floating-feedback-panel` (draggable headers, JSON tree inspector, benchmark comparison table, and feedback templates).
- **Modular Stylesheet Architecture (`webapp/src/style/`)**:
  - Deconstructed monolithic styling into category-focused stylesheets: `layout.css`, `forms.css`, `buttons.css`, `badges.css`, `menus.css`, `panels.css`, `telemetry.css`, and `map.css`.
- **Reusable UI Factory System (`webapp/src/ui/`)**:
  - Implemented typed UI helpers for standard controls: `menu.ts` (`renderDropdownMenu`), `segmentedControl.ts` (`renderSegmentedControl`), `selectControl.ts` (`renderSelect`), `copyBadge.ts` (`renderCopyBadge`), and `copyButton.ts` (`renderCopyButton`), with unit tests in `ui.test.ts`.
- **Frontend Utility Unit Test Infrastructure (`webapp/src/utils/*.test.ts`)**:
  - Integrated `vitest` v5.0.5 for fast, zero-config TypeScript unit testing.
  - Added unit test suites for `dmsParser.test.ts`, `pointParser.test.ts`, `gridBounds.test.ts`, and `exporters.test.ts`.
- **Dynamic Geodetic Help & Benchmark Floating Panel (`GeodeticHelpPanel.ts`)**:
  - Non-modal, draggable floating window presenting geodetic references, official engine replication benchmarks from `/api/v1/grid/comparison`, datum definitions (Stereo 70, ETRS89), and coordinate syntax conventions.
- **Refined Floating Feedback System & REST API Health Synchronization (`FeedbackPanel.ts`, `api/routes/feedback.py`)**:
  - Non-modal draggable feedback panel with 5 report categories (coordinate discrepancy, grid update, bug, feature request, API downtime), client-side bot checks, and optional SMTP email dispatch with local logging fallback.
- **Usage Telemetry Section Harmonization (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`)**:
  - Strict database boundary clamping (`2023-11-01` to current date), neutral Y-axes, 2-line grouped symbology legend, source fallback transparency, and additive Total curve logic.
  - Macroscopic $0.1^\circ$ punchcard canvas matrix with 4 mathematical quadrant summary cards (NW, NE, SW, SE) and interactive hover tooltips.
- **Web Map Ergonomics & Multi-Layer GIS Architecture (`MapView.ts`, `FileUpload.ts`, `map.css`)**:
  - Decoupled map markers into 3 independent layers (`pointLayer`, `textLayer`, `filesLayer`) with live point counts and customizable color swatches.
  - Persistent state-preserving legend dock, segmented GIS graphic scale bar, single-distance measurement ruler, and unified dual-datum coordinate popups with "Copy Attributes" action.

### Changed
- **2-Line Lowercase Header Attribution & External Link Rule (`Header.ts`, `layout.css`, `design_guide.md`)**:
  - Replaced `"by PyTransDatRO"` with stacked lowercase `"powered by"` / `"PyTransDatRO"` link without an oblique arrow (`↗`), reserving `↗` strictly for external links inside dropdown menus.
- **Web GIS Map Widget & Spatial Matrix Theme Synchronization (`map.css`, `SpatialMatrixChart.ts`)**:
  - Synchronized all Leaflet map overlay docks, drawers, scale bars, and geodetic popups with active theme tokens.
  - Added `MutationObserver` theme reactivity to `SpatialMatrixChart.ts` so ghost cells, hover borders, and the 5-tier color ramp adapt cleanly to both Dark and Light themes.
- **Design System Manual Overhaul (`webapp/design_guide.md`, `docs/web_app.md`)**:
  - Rewrote `webapp/design_guide.md` with the GitHub Primer dual-theme palette, Subtle Surface + Border (`#8`) button, 24-item Component Catalogue (`C1`–`C6`), enhanced ON/OFF toggle states, and C2 Graticule logo.
  - Relocated pure algorithmic and parser specifications to `docs/web_app.md`.
- **Decoupled Data Exporters (`webapp/src/utils/exporters.ts`)**:
  - Extracted pure generator functions (`generateCSVContent`, `generateGeoJSONObject`, `generateDXFContent`) from DOM-dependent download triggers, enabling pure headless testing.
- **Coordinate Order & Angular Unit Synchronization (`main.ts`, `PointConverter.ts`, `BatchConverter.ts`)**:
  - Promoted coordinate order (`NE(H)` vs `EN(H)`) to the shared configuration bar, dynamically swapping columns across inputs, outputs, and exports.
  - Decoupled coordinate plotting from auto-panning; explicit Romania-bounded centering (`Center on Map` and `Zoom to Extents`) prevents disruptive viewport jumps.

### Fixed
- **Untransformable Point Handling for "No data on grid" Zones (`PointConverter.ts`, `BatchConverter.ts`, `exporters.ts`)**:
  - Resolved bug where points falling inside the rectangular bounding box but in a no-data zone (`warning: "No data on grid"`) were treated as valid coordinates.
  - Updated `PointConverter.ts` to blank computed outputs (`—`), suppress map plotting, and display `⚠ No Grid Data`.
  - Updated `BatchConverter.ts` to append `# No data on grid` without emitting raw untransformed numbers, counting points as outside grid.
  - Updated `exporters.ts` to preserve `# No data on grid` rows in exported CSVs.
- **Unspaced DMS Parsing Seconds vs South Collision (`webapp/src/utils/dmsParser.ts`)**:
  - Resolved bug where lowercase/uppercase `s` designating seconds in letter DMS notation (e.g., `45d 59m 58.99s`) was misinterpreted as South, incorrectly inverting the coordinate sign.
- **Point Parser Label Classification Heuristic (`webapp/src/utils/pointParser.ts`)**:
  - Fixed regex heuristic where point labels containing `d`, `m`, or `s` (such as `GPS_BASE`, `STATION`, or `Base`) were erroneously treated as coordinate angle symbols.
- **Strict Decimal Validation & 1-to-1 Line Correspondence (`pointParser.ts`, `BatchConverter.ts`)**:
  - Replaced lenient `parseFloat` with strict numeric validation to catch malformed tokens as syntax errors.
  - Guaranteed exact 1-to-1 line correspondence between input and output, preserving empty lines and comments.
- **Single Point Map Render Infinite Recursion Fix (`MapView.ts`)**:
  - Resolved recursive call stack overflow in single-point marker color updates using CSS custom properties `--point-marker-color` and `--point-pulse-color`.
- **File Upload Queue Reset Map Synchronization & 2,000-Point Truncation Fix (`FileUpload.ts`, `main.ts`)**:
  - Synchronized file queue clearing with map layers and removed premature 2,000-point parsing truncation, correctly triggering the `VolumeWarningModal` guardrail.
- **Popup Copy Coordinate HTML Character Escaping (`MapView.ts`)**:
  - Fixed DMS quote character (`'`, `"`) corruption in popup HTML onclick handlers using URI-safe data attributes.

---

## [1.0.0] - 2026-10-10

### Added
- **Interactive Leaflet Web GIS (`MapView.ts`)**:
  - Basemap selector supporting OpenStreetMap, CartoDB Positron, CartoDB Dark Matter, and ESRI World Imagery.
  - Automatic boundary overlay displaying the official Romanian TransDatRO grid bounding box (`[43.5, 20.2]` to `[48.3, 29.8]`).
  - Marker clustering via `leaflet.markercluster` for performant point rendering.
  - Automatic map panning and zooming to bounds when points are transformed.
  - Reverse geodetic click inspection: clicking on the map displays geographic (WGS84/ETRS89) and projected (Stereo 70) coordinates.
- **Single Point Converter (`PointConverter.ts`)**:
  - Live transformation between Stereo 70 and ETRS89 (both directions).
  - Flexible coordinate input format support (Decimal Degrees, Radians, DMS string parsing).
  - High-precision results display with quasigeoid undulation $\zeta$ and grid distortion shifts $\Delta E, \Delta N$.
  - "Show on Map" button and direct marker placement.
- **Batch Converter & Table Preview (`BatchConverter.ts`)**:
  - Multi-point text input with auto-delimiter detection (space, tab, comma, semicolon).
  - Interactive results table with instant map synchronization.
  - Client-side export to CSV, GeoJSON, and AutoCAD DXF.
- **File Upload (`FileUpload.ts`)**:
  - Drag-and-drop file processor with streaming progress bar.
- **Guardrails & Modals**:
  - `VolumeWarningModal.ts`: Performance warning prompt for transformations with > 2,000 points (full vs 1,000-point sample vs skip map).
  - `ReportModal.ts`: Diagnostic calculation report generator.
  - `FeedbackModal.ts`: Discrepancy reporting.
  - `StatsView.ts`: Telemetry visualization and grid model metadata viewer.
- **Design System & Internationalization**:
  - Vanilla CSS design system with Dark / Light / Auto mode support (`variables.css`, `base.css`, `components.css`, `map.css`).
  - Bilingual i18n support for Romanian (`ro`, default) and English (`en`).
- **REST Client & Proxy (`apiClient.ts` / `vite.config.ts`)**:
  - Type-safe REST client communicating with FastAPI endpoints under `/api/v1/`.
  - Vite reverse proxy forwarding `/api`, `/health`, and `/docs` to `http://127.0.0.1:8000`.
