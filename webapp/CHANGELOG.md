# Changelog - PyTransDatRO Web Application (`webapp/`)

All notable changes to the PyTransDatRO web application and its interactive Web GIS interface will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
- **Playwright End-to-End (E2E) Browser Test Suite (`webapp/e2e/`, `playwright.config.ts`)**:
  - Added a comprehensive Chromium Playwright E2E suite covering single-point coordinate transformations, batch text parsing (1-to-1 blank/comment line preservation and `Point N` numbering), file streaming upload/download, map centering and zoom constraints, bi-directional i18n (`EN` / `RO`) persistence, and API error/rate-limit display.
  - Implemented offline request stubbing in `e2e/fixtures.ts` for external map tiles and fonts.
  - Added read-only `window.__e2e.mapView()` diagnostic snapshot in `MapView.ts` (active only in `VITE_E2E=true` builds).

### Fixed
- **Untransformable Point Handling for "No data on grid" Zones (`PointConverter.ts`, `BatchConverter.ts`, `exporters.ts`)**:
  - Resolved an application bug where points falling inside the rectangular grid bounding box but in a no-data zone (`warning: "No data on grid"`) were mistakenly treated as valid transformations by only checking `warning === 'Out of grid'`.
  - Updated `PointConverter.ts` to blank computed output fields (`—`), omit `quasigeoidZeta`, suppress plotting on the map, and display `⚠ No Grid Data` without the false green `✓ Valid coordinates` badge.
  - Updated `BatchConverter.ts` to append `# No data on grid` (or `pointId, # No data on grid`) without echoing raw untransformed input numbers, and count it as outside grid rather than valid.
  - Updated `exporters.ts` to preserve `# No data on grid` rows in exported CSVs.

### Added
- **GitHub Primer Neutral Dual-Theme Palette & Subtle Surface + Border (`#8`) Primary Button (`webapp/src/style/variables.css`)**:
  - Adopted the GitHub Primer neutral Dark (`#0d1117` app canvas, `#161b22` card surface, `#21262d` elevated controls) and Light (`#f6f8fa` app canvas, `#ffffff` card surface, `#f3f4f6` elevated controls) theme tokens, eliminating dark navy-blue backgrounds (`#090d16` / `#0f172a` / `#1e293b`).
  - Set Option `#8` (**Subtle Surface + Border**) as the permanent primary action button style (`--btn-primary-bg: #1f2630`, `--btn-primary-hover: #28313e`, `--btn-primary-text: #e6edf3`, `--btn-primary-border: #48586a` in Dark mode; `--btn-primary-bg: #eaeef2`, `--btn-primary-hover: #dde3ea`, `--btn-primary-text: #1f2328`, `--btn-primary-border: #8c959f` in Light mode).
- **Interactive Computation Chart Legend Color Picker & Dual-Theme Line Palette (`TimelineChart.ts`, `TelemetryDashboard.ts`, `telemetry.css`)**:
  - Replaced the low-contrast pale yellow (`#fde68a` / `#f59e0b`) on `Desktop / CLI` (`desktop_cli`) with high-contrast Coral Rose (`#fb7185`) / Crimson Rose (`#e11d48`), legible on both Dark and Light backgrounds.
  - Added click-to-change 10-swatch color popovers (`.legend-color-popover`, matching the Web Map layer widget UX) on every line item in the Computation chart legend, persisting custom line colors in `localStorage['pytransdat_timeline_colors']`.
- **Enhanced ON/OFF Toggle State Visibility (`forms.css`, `telemetry.css`)**:
  - Added a `6px × 6px` circular ON/OFF status dot (`::before`), elevated `--btn-primary-bg` active surface, and crisp `--btn-primary-border` outline on standalone Usage filter toggles (`REST API`, `Web App`, `Desktop / CLI`, `Points`, `Calls`), while recessed/dimming (`opacity: 0.62`) inactive buttons.
  - Increased ON/OFF contrast on the Transform direction `.sliding-toggle` glider and `.segmented-control` coordinate order pills.
- **Official Dual-Theme Brand Logo (`Map Projection Graticule`, `webapp/src/assets/logos.ts`, `logo.svg`, `public/favicon.svg`)**:
  - Adopted **C2: Map Projection Graticule** as the official dual-theme vector logo, featuring curved stereographic graticule arcs intersecting a projection frame with a central origin node, rendered inline via CSS variables (`getLogoInlineSvg(36)`) and synchronized dynamically with the browser tab favicon (`syncFaviconWithTheme`).
- **Status-Displaying Language Selector Dropdown Menu (`Header.ts`, `ui/menu.ts`, `menus.css`)**:
  - Converted the header Language control into a disclosure Dropdown Menu (`#header-lang-details`, `renderDropdownMenu`) with a globe SVG icon, active language code (`RO` / `EN`), rotating chevron (`▾`), and `Română (RO)` / `English (EN)` menu items (`role="menuitemradio"`, `aria-checked`) with an active checkmark (`✓`, `.menu-check-icon`), visually distinguishing status-displaying menus from the immediate-action Theme Toggle button.
- **`RomGEO Web API` Link & `RomGEO` Capitalization (`Header.ts`, `i18n/en.ts`, `i18n/ro.ts`)**:
  - Added `RomGEO Web API` (`https://api.romgeo.ro/api/v1/`, `isExternal: true`) at the bottom of the `ANCPI / ROMPOS Resources` section in the Help menu and standardized all `RomGeo` references to `RomGEO`.
- **Dedicated Floating Windows Stylesheet (`webapp/src/style/panels.css`)**:
  - Restored and unified full CSS rules for `.floating-grid-panel`, `.floating-help-panel`, and `.floating-feedback-panel` (draggable headers, JSON tree inspector, benchmark comparison table, and feedback templates).
- **Modular Stylesheet Architecture (`webapp/src/style/`)**:
  - Deconstructed the monolithic 3,561-line `components.css` into maintainable, category-focused stylesheets:
    - `layout.css`: Header, workbench layout, tabs bar, studio cards, modals, and toast notifications.
    - `forms.css`: Form groups, labels, inputs, selects, sliding direction toggle, segmented switches, format popovers, and file dropzone.
    - `buttons.css`: Theme-native primary CTA, secondary, outline, copy buttons, action buttons, and format pills.
    - `badges.css`: Status pills, Quasigeoid height metric field, angle format badges, file status pills, health dots, and quality report cards.
    - `menus.css`: Top navigation menu item triggers (`.nav-menu-trigger`), dropdown menu panels, active checkmarks, menu items, category headings, and dividers.
    - `telemetry.css`: KPI cards, controls decks, interactive SVG timeline charts, legend color popovers, and spatial matrix canvas deck.
- **Reusable UI Factory System (`webapp/src/ui/`)**:
  - `menu.ts`: Generic `renderDropdownMenu` helper supporting sleek top navigation menu item triggers (`.nav-menu-trigger`), rotating chevron animations, optional `active` radio checkmarks (`role="menuitemradio"`), custom `menuWidth`, ARIA semantics (`role="menu"`, `role="menuitem"`), keyboard navigation (`Escape`, `ArrowDown`, `ArrowUp`), and automated dismissal.
  - `segmentedControl.ts`: Accessible `renderSegmentedControl` helper for radio-switch toggles with ARIA attributes and arrow-key cycling.
  - `selectControl.ts`: `renderSelect` helper ensuring uniform styling and focus rings.
  - `copyBadge.ts`: `renderCopyBadge` helper for interactive copyable chips with keyboard activation.
  - `copyButton.ts`: `renderCopyButton` helper providing standard code-block copy actions with animated checkmark feedback.
  - `ui.test.ts`: Added unit tests covering all UI factories and the dual-theme `Map Projection Graticule` brand logo (67/67 tests passing overall).

### Changed
- **2-Line Lowercase Header Attribution & External Link Rule (`Header.ts`, `layout.css`, `design_guide.md`)**:
  - Replaced `"by PyTransDatRO"` with a 2-line stacked lowercase `"powered by"` / `"PyTransDatRO"` attribution link without an oblique arrow (`↗`), positioned in the left header zone right after the brand block.
  - Updated the external link indicator rule so `↗` (`.ext-icon`) applies strictly to external links inside dropdown menus (e.g., Help menu links to ROMPOS, RomGEO, and `/docs`).
- **Web GIS Map Widget & Spatial Matrix Theme Synchronization (`map.css`, `SpatialMatrixChart.ts`)**:
  - Synchronized all Leaflet map overlay docks, drawers, scale bars, and geodetic popups with the active theme tokens (`--bg-glass`, `--bg-surface`, `--border-medium`, `--radius-*`).
  - Added `MutationObserver` theme reactivity to `SpatialMatrixChart.ts` so ghost cells, hover borders, and the 5-tier color ramp adapt cleanly to both Dark and Light themes.
- **Design System Manual Overhaul (`webapp/design_guide.md`, `docs/web_app.md`)**:
  - Rewrote `webapp/design_guide.md` with the GitHub Primer dual-theme palette, the `Subtle Surface + Border` (`#8`) primary button, the 24-item Component Catalogue (`C1`–`C6`), enhanced ON/OFF toggle states, interactive Computation chart legend color pickers, header navigation rules, and the `Map Projection Graticule` dual-theme logo while preserving the `480px` workbench width and `1.35fr 1.35fr 0.8fr` coordinate input geometry, and relocated pure algorithmic/parser specs to `docs/web_app.md`.
- **Header Help Menu Redesign (`Header.ts`, `menus.css`)**:
  - Redesigned the Help dropdown trigger from an isolated `.btn.btn-outline` button to a seamless `.nav-menu-trigger` menu item, giving it clean top-navigation aesthetics while maintaining full accessibility and dropdown capability.
- **Component Unification (`main.ts`, `PointConverter.ts`, `BatchConverter.ts`)**:
  - Refactored `main.ts` shared config bar to use `renderSelect`, `renderSegmentedControl`, and `renderCopyBadge`.
  - Refactored `PointConverter.ts` and `BatchConverter.ts` to replace repetitive manual copy button DOM manipulation with `renderCopyButton`.
- **Frontend Utility Unit Test Infrastructure (`webapp/package.json`, `webapp/src/utils/*.test.ts`)**:
  - Integrated `vitest` v5.0.5 for ultra-fast, zero-config TypeScript unit testing directly using the existing Vite setup.
  - Added `dmsParser.test.ts` (13 tests) validating angular conversions, symbol variations (`°`, `'`, `"`), letter notations (`d/m/s` and `D/M/S`), whitespace separation, cardinal directions, and seconds rollover formatting.
  - Added `pointParser.test.ts` (28 tests) validating decimal verification, automatic delimiter detection (comma, semicolon, tab, space), coordinate order normalization (`NE` vs `EN`), 2D vs 3D detection, point ID preservation, and blank/comment line tracking.
  - Added `gridBounds.test.ts` (10 tests) validating the official Stereo 70 and ETRS89 Romania bounding envelopes against national city coordinates, edge thresholds, and exterior points.
  - Added `exporters.test.ts` (8 tests) validating structured exports for CSV (delimiters, coordinate precision), GeoJSON (`FeatureCollection`, `Point` geometries, properties), and AutoCAD ASCII DXF (`POINT` entities on `POINTS` layer, `TEXT` labels on `LABELS` layer).

### Changed
- **Decoupled Data Exporters (`webapp/src/utils/exporters.ts`)**:
  - Extracted pure generator functions (`generateCSVContent`, `generateGeoJSONObject`, `generateDXFContent`) from DOM-dependent download triggers, enabling pure headless testing while retaining backward-compatible wrappers (`exportToCSV`, `exportToGeoJSON`, `exportToDXF`).

### Fixed
- **DMS Parser Seconds vs South Indicator Collision (`webapp/src/utils/dmsParser.ts`)**:
  - Resolved bug where lowercase/uppercase `s` designating seconds in letter DMS notation (e.g. `45d 59m 58.99s`) was misinterpreted as South, incorrectly inverting the coordinate sign.
  - Added case-insensitive splitting for `D`, `M`, `S` components.
- **Point Parser Label Classification (`webapp/src/utils/pointParser.ts`)**:
  - Fixed regex heuristic where point labels containing the characters `d`, `m`, or `s` (such as `GPS_BASE`, `STATION`, or `Base`) were erroneously treated as coordinate angle symbols.
- **Refined Floating Feedback System & REST API Health Synchronization (`FeedbackPanel.ts`, `components.css`, `api/routes/feedback.py`, `api/services/mailer.py`, `en.ts`, `ro.ts`)**:
  - Replaced the static `FeedbackModal.ts` with a non-modal, draggable, resizable floating panel matching `GridMetadataPanel.ts` and `GeodeticHelpPanel.ts`.
  - Dropped Name and Contact Email fields for complete GDPR exemption and user privacy.
  - Removed the standalone "Tested Coordinates" box to streamline form ergonomics.
  - Added 5 specialized report categories:
    - **Coordinate Discrepancy** (`discrepancy`): Pre-fills input, PyTransDatRO result, and expected benchmark coordinates.
    - **New Official Grid Release** (`grid_update`): Auto-detects and injects the currently loaded SPG grid metadata, prompting the user for newer official ANCPI/RomGeo versions.
    - **Bug Report** (`bug`): Pre-fills issue summary, reproduction steps, expected/observed behavior, and browser environment.
    - **Feature / Format Suggestion** (`feature`): Solicits use cases and format specifications.
    - **REST API Offline / Down** (`api_down`): Directly synchronized with the header's API health indicator dot and offline alert banner for 1-click incident reporting.
  - Added lightweight bot protection: off-screen honeypot field, in-browser dynamic math verification challenge, and fast-submit prevention.
  - Added backend `POST /api/v1/feedback` endpoint with standard library SMTP email dispatch and local `logs/feedback.jsonl` logging fallback.
  - Added 1-click "Copy Report" and "Open GitHub Issue" actions for zero-friction user sharing.
- **Dynamic Geodetic Help & Benchmark Floating Panel (`GeodeticHelpPanel.ts`, `components.css`, `apiClient.ts`, `en.ts`, `ro.ts`)**:
  - Replaced the static, modal-based help dialog with a non-modal, draggable floating window matching the UX, dragging mechanics, boundary clamping, and resizing capabilities of `GridMetadataPanel.ts`.
  - Organized geodetic documentation into 4 clean, scannable cards adhering strictly to [`design_guide.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/design_guide.md):
    - **Card 1: About & Disclaimer**: States that PyTransDatRO is an independent open-source MIT project reproducing ANCPI TransDatRO & RomGeo algorithms, accompanied by an explicit non-official disclaimer.
    - **Card 2: Official Engine Replication Benchmark**: Dynamically pulls test metrics across 9,046 test points covering Romania from `/api/v1/grid/comparison`, displaying sub-millimeter residuals ($\le 0.34\text{ mm}$ planimetric, $0.5\text{ mm}$ vertical), dynamic PASS status, and explicit tolerance criteria ($\le 0.003\text{ m}$ / $0.00003''$).
    - **Card 3: Reference Systems & Datums**: Crisp definitions of Stereo 70 (`EPSG:31700`) with normal heights $H$ explicitly in meters, ETRS89 (`EPSG:4258`), and the Quasigeoid $\zeta = h - H$ model.
    - **Card 4: Formats & Conventions**: Quick syntax guide for unspaced DMS (`46°00'00.00000"`), decimal degrees, radians with note that angle format governs input parsing, coordinate display, and exported batch files; `NE(H)` vs `EN(H)` coordinate order; and 2D/3D full transformation auto-detection.
- **ANCPI Comparison Report Pipeline & Backend Integration (`tools/generate_ancpi_comparison_reports.py`, `api/routes/info.py`, `api/schemas/info.py`)**:
  - Renamed `tools/generate_reports.py` to `tools/generate_ancpi_comparison_reports.py`.
  - Preserved human-readable `.txt` summaries and interactive `.html` maps in `reports/`.
  - Added export of machine-readable benchmark JSON written strictly once to `pytransdatro/grids/ancpi_comparison_summary.json`.
  - Exposed `GET /api/v1/grid/comparison` REST endpoint for dynamic frontend consumption.

### Changed
- In the Header and Navigation chrome:
  - Consolidated actions: Moved standalone **REST API Reference** (`/docs` with `↗`) and **Feedback** into the **Help** dropdown under the `PyTransDatRO` group.
  - Positioned **Feedback** above **REST API Reference**, separated by a clean dashed line (`.dropdown-divider-dashed`) to visually distinguish internal modal/panel actions from external documentation.
  - Formatted **REST API Reference** with left-aligned text, no internal icon, and a right-aligned external link arrow (`↗`), maintaining 100% visual consistency with the other external links in the dropdown.
  - Refined the brand subrow health indicator: Guaranteed a bright, circular green dot (`7.5px`) with glow when online, turning red with an inline alert message `Report REST API failure` when an issue is detected, linking directly to the feedback panel with `api_down` pre-selected.
  - Aligned Feedback form with [`design_guide.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/design_guide.md): uppercase micro-headers (`REPORT TYPE`, `DETAILED DESCRIPTION`), dedicated template note beneath the title, and placed "Reset Template", "Copy Report", and "Open GitHub Issue" on their own dedicated row aligned to the right.
  - Renamed group title from **HELP** / **Ajutor** to **PyTransDatRO**.
  - Renamed link text from **Help** / **Ajutor** to **PyTransDatRO Help** / **Ajutor PyTransDatRO**.
  - Replaced the circular icon with a crisp standalone SVG question mark (`?`) help icon without an enclosing circle for maximum legibility.
- Removed badge tag from Help floating window header for a cleaner aesthetic.

### Removed
- Removed legacy `FeedbackModal.ts` in favor of floating, draggable `FeedbackPanel.ts`.
- Removed legacy `GeodeticHelpModal.ts` in favor of `GeodeticHelpPanel.ts`.
- Removed legacy `tools/generate_reports.py` in favor of `tools/generate_ancpi_comparison_reports.py`.

### Added
- **Usage Section Iteration 2 Harmonization & Ergonomics (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`, `components.css`, `en.ts`, `ro.ts`)**:
  - **Source Fallback Transparency**: Unmapped database client sources display their exact database string (`sourceName`) directly in toggles and legends without transformation.
  - **Visual Symbology Brightness Hierarchy**: Standardized brightness across all chart lines and legend symbols: **Calls** use lighter tints (pastel/luminous), while **Points** use darker tones of the same hue family.
  - **Unified Deck Ergonomics & Filter Layout**:
    - Aligned time presets to the left and granularity toggles (`Daily`, `Weekly`, `Monthly`) to the right (`justify-content: space-between`).
    - Merged **SOURCES** and **METRIC** into a single cohesive row across both views (sources on left, metric on right; concise labels `SOURCES` and `METRIC`).
    - Harmonized filter controls: SOURCES (both Computation and Spatial) and METRIC (Computation) are styled as individual standalone toggle buttons with theme-native neutral borders and backgrounds (`--bg-surface-elevated`), generous 8px spacing, and zero continuous frame enclosure, completely removing blue borders and checkboxes. Active states use high-contrast theme-native `--btn-primary-bg`/`--btn-primary-text`. METRIC on Spatial retains the mutually exclusive segmented pill toggle.
    - Updated active button styling to theme-native high contrast (`--btn-primary-bg` / `--btn-primary-text`), removing electric blue tints.
  - **Strict Database Boundary Clamping**: Clamped calendar date pickers (`startVal < minHistoryDate` and `endVal > maxHistoryDate`) directly in input handlers and DOM attributes, preventing selection of dates prior to earliest database record (`2023-11-01`).
  - **Neutral Computation Chart Axes**: Neutralized left and right Y-axes tick marks and title labels to match the neutral grey of the X-axis (`var(--text-muted)`).
  - **2-Line Grouped Symbology Legend**: Redesigned the Computation graph legend into a clean 2-line layout: Line 1 groups active sources by client name (`Web App: Calls, Points`), Line 2 displays `Total: Calls, Points` on its own line labeled simply "Total".
  - **Deduplicated Spatial Top Banner**: Kept solely the clean numerical Active Cells metric in the top overview banner, eliminating duplicate cell area text.
  - **Stacked Spatial Heatmap Footnotes**: Split and stacked the bottom notes into two distinct, clear items: plain grid cell dimensions (`0.1° latitude by 0.1° longitude (~10 km, ~86 km² avg. cell area at 46°N)` without the word "Macroscopic") and quadrant divider explanation.
  - **Heatmap Latitude Axis Alignment**: Relocated the Latitude axis to the left side with right-aligned coordinate text, eliminating the clipped letter "N" at the canvas border.
- **Usage Section Visual & Architectural Refinement (`TelemetryDashboard.ts`, `TimelineChart.ts`, `en.ts`, `ro.ts`, `components.css`, `variables.css`)**:
  - **Header Structure & Navigation**:
    - Removed the redundant "Back to Geodetic Studio" button; changed page title to "Usage Data" (`Date de Utilizare`).
    - Styled top view switcher tabs with GitHub-style underline tabs ("Computation" and "Spatial") adjacent to the title with a subtle separator, matching the primary `Transform` and `Usage` header navigation. Completely removed the obsolete "Grid Specs" tab.
  - **Computation View**:
    - Replaced the `noUiSlider` timeline scrubber with native themed `<input type="date">` calendar controls matching the Spatial section, including range limits pill (`Earliest Available` • `Current Date`).
    - Applied `color-scheme: dark | light` ensuring native browser calendar pickers and form controls render natively in dark or light theme.
    - Dynamically discovered sources from the database mapped via canonical lookup (`REST API`, `Web App`, `Desktop / CLI`), cleanly falling back to capitalized database names for future sources, and filtering out `unknown`.
    - Implemented strictly additive Total curve logic across active selected sources ($\text{Total} \ge \text{each active source}$).
    - Enhanced Total curve symbology with wider dash spacing (`stroke-dasharray="7 5"`), thicker stroke (`2.6px`), and visible 24px repeating dashed sample in the legend.
    - Added 3-way Metric selector: `Points`, `Calls`, or `Both`, dynamically adapting curves, left/right Y-axes with explicit title labels (`Calls` on left, `Points` on right), and interactive tooltips.
    - Standardized terminology to "Points" and "Calls". Converted top KPI cards to neutral technical styling, removing saturated green accents.
    - Separated source toggle checkboxes (clean minimalist checkboxes without attached line symbols) from the dedicated graph symbology legend.
    - Dynamically generated graph symbology legend synchronized with the visible curves and active metric filter.
  - **Spatial View**:
    - Harmonized active period controls into a unified compact preset bar matching the Computation view.
    - Redesigned top overview banner into a full-width summary card spanning from Total Points/Calls (left) to Active Cells (right) with the average cell size indicator.
    - Optimized quadrant summary cards with vertical stacking (quadrant name above, volume below) and compact adjacent percentage (`1,250 (34.2%)` without redundant "of total" text).
    - Verified axis coordinate labels use `°N` and `°E` cartographic notation.
    - Relocated scaling mode selector (`Logarithmic` | `Linear`) to the bottom-right of the heatmap canvas viewport alongside the 5-step discrete color ramp.
    - Combined grid resolution, average cell area, and quadrant divider notes into a single cohesive footnote below the heatmap.
- **Refactored 3-Zone Application Header (`Header.ts`, `GridInfoModal.ts`, `logo.svg`, `components.css`, `ro.ts`, `en.ts`)**:
  - **Zone 1: Brand Cluster Hierarchy & Favicon**:
    - Enlarged vector logo to `36px × 36px` to harmonize with a two-line text cluster.
    - Title `PyTransDatRO` stacked above a balanced metadata subrow: `v1.0.0-beta` left-aligned and the glowing API health indicator dot (`7.5px`) right-aligned to match the title boundary, omitting the redundant word "Online".
    - Updated `webapp/public/favicon.svg` to match the geodetic shift lattice logo in the browser tab.
  - **Zone 2: GitHub-Style Underline App Navigation**:
    - Relocated `Transform` and `Usage` directly adjacent to the brand cluster, delimited by a subtle vertical separator.
    - Converted buttons into clean text tabs featuring an active underline indicator highlight and hover transitions.
  - **Zone 3: Consolidated Help Dropdown, Repository Attribution & Uniform Buttons**:
    - **Native `<details><summary>` Disclosure Menu**: Replaced fragile JavaScript event listener class toggling with the browser's native HTML5 `<details><summary>` disclosure pattern, guaranteeing 100% reliable opening/expansion on click across all devices and input methods with seamless click-outside, item click, and Escape-key closure.
    - Moved `by PyTransDatRO` external repository link (`https://github.com/bogdanmorosanu/PyTransDatRO`) into the right action cluster.
    - Standardized outline button appearance across `Help`, `API`, and `Feedback`, removing the text underline glitch on the `API` link.
  - **Zone 4: Zero-Emoji System Utilities**:
    - Replaced emoji theme toggle with clean SVG geometric icons (sun / moon) adhering to `design_guide.md`.
    - Added clean `RO` / `EN` segmented language toggle with full reactive i18n support.
- **Spatial Overview Matrix Visual (`SpatialMatrixChart.ts`, `TelemetryDashboard.ts`, `apiClient.ts`, `components.css`)**:
  - **Macroscopic $0.1^\circ$ Punchcard Canvas**: Renders a fixed, non-zoomed macroscopic grid of Romania on HTML5 canvas (~10 km bins), featuring faint ghost cells for unvisited coordinates and a high-contrast 5-tier discrete color ramp.
  - **Mathematical Quadrant Summary Cards**: Four dynamic summary cards (**NW**, **NE**, **SW**, **SE**) calculated around the national center ($45.9^\circ\text{ N}, 24.9^\circ\text{ E}$), displaying absolute volumes and relative percentage shares without administrative boundary bias.
  - **Subtle Quadrant Divider Grid Lines**: Rendered subtle dashed crosshair lines along Latitude $45.9^\circ\text{ N}$ and Longitude $24.9^\circ\text{ E}$ across the heatmap canvas, directly visualizing the 4 quadrants.
  - **Physical Cell Area Indicator**: Integrated an informational badge calculating average cell coverage at projection center ($46^\circ\text{ N}, 25^\circ\text{ E}$): $\approx 86\text{ km}^2$ ($11.1 \times 7.7\text{ km}$).
  - **Vanilla Graphical Date Window Scrubber**: Added dependency-free date range controls with synchronized start and end dates bound between earliest database history (`2023-11-01`) and the current date.
  - **Design System Alignment & Strict English i18n**: Eliminated all decorative emojis from tabs and action buttons adhering to `design_guide.md`, and completed 100% pure English localization across all labels, quadrant shares, notes, and discrete legend steps.
  - **Metric Selection & Scale Mode Toggles**: Live toggling between **Coordinates Transformed (Points)** and **Service Invocations (Calls)**, as well as **Logarithmic (default)** vs. **Linear** scaling with an interactive color ramp legend.
  - **Independent Temporal & Source Filtering**: Dedicated time horizon presets (`All Time`, `12M`, `YTD`, `30D`, `7D`) and dynamic source filter chips decoupled from the timeline visual.
  - **High-Contrast Cell Inspection Tooltip**: Highlighting hovered cells with exact coordinates, bounds, quadrant, and volume breakdown.
- **Interactive Usage Data & Telemetry Dashboard (`TelemetryDashboard.ts`, `TimelineChart.ts`, `noUiSlider`)**:
  - **Draggable Time Window Slider (`noUiSlider`)**: Replaced date pickers with a dual-handle range slider with `behaviour: 'drag'`. Users can slide a custom time window across the entire timeline bounded between the earliest log timestamp and the current date, emulating GIS-style timeline scrubbers.
  - **Dynamic Source Discovery**: Sources are discovered dynamically from the SQLite `sources` table via `/api/v1/telemetry/timeline`. Any new logging client source added to the database automatically appears as an interactive toggle checkbox without frontend code changes.
  - **Dynamic Total & KPI Recalculation**: The Total line on the chart and the top KPI summary cards (Total Invocations, Total Coordinates, Batch Intensity) dynamically recalculate based strictly on the currently active/checked client sources.
  - **Weekly Granularity (`Săptămânal`)**: Added a weekly aggregation button between Daily and Monthly on both backend (SQLite Monday week-start grouping) and frontend.
  - **Accurate Period Count**: Switching between Daily, Weekly, and Monthly correctly recomputes the exact count (`X zile`, `X săptămâni`, `X luni`) corresponding to the active date interval.
  - **Visual Distinction for Metrics**: Clean dual-metric visualization separating Invocations Count (solid primary line) and Coordinates Processed (lighter secondary line) alongside the dashed dynamic Total curve.
- **Persistent Basemap Selection (`MapView.ts`)**: User-selected basemap is now saved to `localStorage['pytransdat_basemap']` and restored across page reloads (defaulting to OpenStreetMap on initial launch).
- **Fullscreen Synchronization (`MapView.ts`)**: Added a global `fullscreenchange` listener on `document` to synchronize the controls dock button state and trigger Leaflet `invalidateSize()` whenever exiting or entering fullscreen mode (including via Escape key).
- **Z-Index Layer Token (`variables.css`, `map.css`)**: Introduced `--z-map-controls: 500;` to ensure floating map docks and badges are consistently layered below application modals and dropdowns.
- **Custom Segmented GIS Graphic Scale Bar (`MapView.ts`, `map.css`)**: Replaced Leaflet single-line scale with a GIS-standard alternating segmented graphic bar. Switched to a high-contrast monochrome cartographic scheme (crisp white & deep obsidian in dark mode, solid black & clean white in light mode), and simplified distance labels to show only start (`0`) and end (`[max] m/km`) values to avoid text crowding while preserving the midpoint division tick.
- **Persistent State-Preserving Legend & Scale Bar Reflow (`MapView.ts`, `map.css`)**: Converted the Legend widget to a persistent stateful dock: expanding it keeps it open indefinitely while navigating, zooming, or transforming coordinates without auto-closing on map clicks. Embedded the legend panel in normal flex flow above the scale bar widget so expanding the legend naturally pushes the scale bar down rather than covering it. Added an in-panel close button, active state button styling, and persisted the user's open/collapsed preference in `localStorage['pytransdat_legend_open']`.
- **Layers Drawer Structure & Basemap Overhaul (`MapView.ts`)**: Reordered drawer groups into: 1. Data Layers, 2. Shift Grids (Coming Soon), 3. Base Maps. Redesigned Base Maps list using standard radio items (`.map-layer-item`) identical in styling to data layers. Set OpenStreetMap as the default active basemap. Updated the Layers drawer toggle to a clean 32x32px icon-only button matching the Legend toggle.
- **Compact Geographic Coordinates Display (`MapView.ts`)**: Abbreviated geographic coordinate labels to `Lat: [val]  Long: [val]` in the bottom-center badge for a tighter, unobtrusive footprint.
- **Uniform Multi-Layer Popups (`MapView.ts`, `FileUpload.ts`)**: Standardized popup tables across Point, Text, and Files layers: Point Name/Number, (Source File for Files layer), Northing, Easting, Altitude, Latitude, Longitude, Height, Quasigeoid Height, and Grid Coverage. Paired input and transformed output coordinates in `FileUpload.ts` (`parsePairedPoints`) to populate both Stereo 70 and ETRS89 values on every file point. Ensured geographic coordinates strictly honor the active `Angular Unit` (Degrees, DMS, Radians) at the time the transformation was executed across all layers.

### Changed
- **Geodetic Altitude ($H$) vs Ellipsoidal Height ($h$) Separation (`MapView.ts`, `PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)**: Popups and data pipelines now distinguish Normal Altitude ($H$, Black Sea 1975 datum) for Stereo 70 coordinates and Ellipsoidal Height ($h$, GRS80 ellipsoid) for ETRS89 coordinates ($h = H + \zeta$). Values are rendered independently and no longer duplicate the same number under both labels.
- **Standardized "Copy Attributes" Action (`MapView.ts`, `ro.ts`, `en.ts`)**: Renamed the popup clipboard action from "Copy Coordinates" to "Copy Attributes" (`Copiază atributele`), formatting copied strings as consistent key-value pairs (`Key: Value`) including `Point Name/Number` for all layers.
- **Leaflet Attribution Styling (`map.css`)**: Styled `.leaflet-control-attribution` with a dark frosted glass appearance matching the application's Obsidian theme, eliminating the glaring white background box.
- **Mobile Responsive Coordinates Badge (`map.css`)**: Elevated the bottom-center coordinates badge on mobile screens ($\le 768$px) to prevent collision with the bottom-left Layers button and map attribution.

### Fixed
- **Workbench Panel Fixed Width Stabilization (`base.css`)**: Imposed strict `480px` dimensions (`width: 480px; min-width: 480px; max-width: 480px; box-sizing: border-box;`) on `.utility-panel` and constrained the split workbench grid column (`grid-template-columns: 480px minmax(0, 1fr)`). This eliminates column stretching and prevents control overlap or UI shifting when switching between the Point, Text, and Files tabs.
- **Layers Drawer Event Listener Leak (`MapView.ts`)**: Moved the global document click listener for `#map-layers-panel` and color popovers into one-time initialization (`initCustomDocks()`), preventing listener accumulation on layer renders.
- **File Upload Queue Reset Map Synchronization (`FileUpload.ts`, `main.ts`)**: Clearing the file queue or removing individual queued files now invokes `onClearFiles` / `mapView.clearFilePoints()`, keeping the map layers and cluster groups synchronized with the file list.
- **File Upload 2,000-Point Truncation Fix (`FileUpload.ts`)**: Removed a hardcoded `Math.min(..., 2000)` limit in `parsePairedPoints` that prematurely truncated uploaded coordinate sets before map plotting. All points are now fully parsed and passed to `MapView.plotFilePoints`, properly triggering the `VolumeWarningModal` guardrail when datasets exceed 2,000 points and enabling full rendering of large datasets (e.g. 9,000+ points).
- **Popup Copy Coordinate HTML Attribute Glitch (`MapView.ts`)**: Resolved an issue where DMS coordinates containing `'` and `"` characters corrupted inline `onclick` string literals and spilled raw JS text into the popup. Implemented safe URI encoding (`data-copy-text`) coupled with a Leaflet `popupopen` event listener and temporary `✓ Copied` feedback.
- **Single Point Map Render Infinite Recursion Fix (`MapView.ts`)**: Resolved an infinite recursive loop (`Maximum call stack size exceeded`) triggered when clicking "Transform" on the Point tab, where `updateSinglePoint` invoked `applyLayerColors('point')` which in turn re-invoked `updateSinglePoint`. CSS custom properties `--point-marker-color` and `--point-pulse-color` now dynamically restyle the single-point pin without re-triggering marker updates.

### Changed
- **Web Map Standardization & Refinements - Iteration 2 (`MapView.ts`, `FileUpload.ts`, `map.css`, `main.ts`)**:
  - **Separate Persistent Legend & Collapsible Layers Widgets**: Decoupled the legend from the layers panel; the Legend widget was separated into its own dock showing active layers and live color dots, while the Layers widget provides a dedicated drawer for toggling visibility, selecting basemaps, and picking colors.
  - **Full Marker & MarkerCluster Color Synchronization**: Dynamic cluster icons (`custom-cluster-circle`) now inherit the exact active swatch color for that layer, fully synchronizing individual markers, cluster grouping bubbles, Layers drawer swatches, and Legend dots.
  - **Shift Grids Section (Coming Soon)**: Restored the 4 planned shift grid layers ($\Delta E$, $\Delta N$, Total Vector, $\zeta$) in the Layers drawer with coming soon badges.
  - **Restructured Themed Popups**: Popups now feature the Layer Name as the prominent header, `Point Name/Number` as an attribute, source `File` name for uploaded files, unified coordinate names without unit suffixes (`Northing`, `Easting`, `Altitude`, `Latitude`, `Longitude`, `Height`, `Quasigeoid Height`), and bottom-placed `Grid Coverage` status.
  - **Smooth Inline Copy Feedback**: Replaced disruptive `alert()` dialogs in popups with an in-place label transition to `✓ Copied` for 1.5 seconds.
  - **Bottom-Center Reactive Coordinates Widget**: Repositioned the coordinates badge to bottom-center (avoiding attribution overlap) and wired it to reactively format coordinates per the active `Angular Unit` (Degrees, DMS, Radians) without unit clutter or click actions.
  - **Streamlined Single-Distance Measurement Tool**: Simplified the ruler tool to measure one clean distance at a time (`Distance: [value] km`), auto-resetting on subsequent clicks without cumulative total confusion.

### Removed
- **TransDatRO Official Grid Boundary Layer**: Removed the approximated grid boundary polygon overlay and reporting from the web map, eliminating misleading representations.
- **Web Map Ergonomics & Controls Overhaul (`MapView.ts`, `map.css`, `main.ts`)**:
  - **Themed Geodetic Popups**: Fully custom Dark Obsidian / Light Porcelain styled popups with structured dual-datum coordinate tables (Stereo 70 Northing/Easting and ETRS89 Latitude/Longitude), orthometric/ellipsoidal height, and technical cyan Quasigeoid height $\zeta$, plus 1-click clipboard copy. Single-point marker also gains a rich interactive popup.
  - **Top-Left Custom Controls Dock**: Grouped floating dock housing Zoom In/Out, **Zoom to Romania (Home)**, **Metric Distance Measurement Ruler**, and **Fullscreen Toggle**.
  - **Distance Measurement Ruler**: Native Leaflet geodesic measuring tool calculating baseline distance in meters/kilometers with dynamic bearing and a dedicated clear action.
  - **Bottom-Left Grouped Layers & Legend Dock**: Relocated layer panel grouping Basemaps (Dark Matter, Positron, OSM, Satellite), Data Layers with live point count badges, Grid Boundary overlay, and a **6-swatch Color Picker** per data layer (Cyan, Cobalt, Emerald, Amber, Violet, Rose) with `localStorage` persistence.
  - **CartoDB Positron Basemap**: Clean, minimalist light-mode basemap preventing background visual noise.
  - **Bottom-Right Status & Metric Scale Dock**: Relocated coordinates badge with 1-click clipboard copy, aligned alongside an official Leaflet **Metric Scale Bar**.
  - **CAD-Style Point Labels at High Zoom**: Subtle point ID labels displayed at street/field zoom levels ($\ge 15$) for instantaneous point identification without opening popups.
  - **Global Zoom Range**: Reduced `minZoom` to 3, allowing full regional/global zoom-out with 1-click return to Romania via the Home button.
- **Non-Disruptive Map Click Interaction (`main.ts`)**:
  - Guarded map click coordinate capture strictly to when the **Point** tab is active; clicking the map while reviewing Batch Text or Files tabs no longer causes disruptive tab switching.
- **Mobile Viewport Invalidate Fix (`main.ts`, `MapView.ts`)**:
  - Added `mapView.invalidateSize()` invocation upon toggling between workbench drawer and map panel on mobile screens, eliminating grey or blank unrendered tile glitches.
- **Separate Center on Map & Zoom to Extents Actions Across All Tabs (`MapView.ts`, `PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)**:
  - Implemented dual map navigation controls on Point, Text, and Files tabs: `Center on Map` pans to the center of valid in-grid coordinates preserving the user's current zoom level; `Zoom to Extents` pans and fits bounds to encompass all points (or zooms to level 14 for single point).
- **SPG Envelope Pre-Flight Validation (`gridBounds.ts`)**:
  - Added coarse pre-flight envelope boundaries derived from the union of 2D planimetric and 1D vertical SPG grids (`[213600, 785640]` m Northing, `[101260, 890790]` m Easting; `[43.30°, 48.58°]` Lat, `[19.70°, 30.30°]` Lon) preventing false calls or incorrect panning.
- **Backend File Streaming Robustness (`api/routes/file.py`, `api/schemas/point.py`)**:
  - Added `DMS = "dms"` to backend `AngleUnitEnum`.
  - Implemented automatic Point ID column detection (handles 3 or 4 columns with optional ID), DMS coordinate parsing with sexagesimal radian conversion, unspaced DMS formatting, and Windows-1252 (ANSI) degree symbol fallback decoding.
- **Automated 200-Point Test Datasets & Verification (`tests/data/`, `tests/test_api_modern.py`)**:
  - Created 200-point test files and added regression tests `test_modern_file_streaming_200_points_st70_to_etrs89` and `test_modern_file_streaming_200_points_etrs89_to_st70_dms`.
- **3 Distinct Leaflet Map Layers & Decoupled Centering (`MapView.ts`)**:
  - Split map markers into 3 independent Leaflet layers (`pointLayer`, `textLayer`, `filesLayer`), exposed in the Leaflet layer control for granular visibility toggling.
  - Decoupled coordinate plotting from auto-panning/zooming; explicit Romania-bounded centering methods (`centerOnPoint`, `zoomToTextExtents`, `zoomToFilesExtents`) are triggered only via `Center on Map` buttons.
- **Theme-Native High-Contrast Primary Buttons (`variables.css`, `components.css`)**:
  - Implemented high-contrast porcelain/charcoal (dark mode) and slate/white (light mode) styling for `.btn-primary` actions.
- **Dedicated Quasigeoid Height Metric Box & Persistent Output Layout (`PointConverter.ts`, `components.css`)**:
  - Promoted Quasigeoid Height ($\zeta = h - H$) to a dedicated metric field row (`.metric-field`) with dedicated metric token styling.
  - Standardized pure numeric output fields without trailing unit suffixes and ensured persistent `—` layout before initial calculation.
  - Default initial coordinates set to Romania's projection center origin ($\phi = 46.0^\circ\text{ N}, \lambda = 25.0^\circ\text{ E}, h = 100.0\text{ m}$) with default direction `ETRS89ToStereo70`.
  - Smart map click captures coordinates directly for ETRS89 or reverse-calculates Stereo 70 meters.
- **Batch Text Copy Results Header Button (`BatchConverter.ts`)**:
  - Added dedicated `Copy Results` icon button on the output card micro-header with instant clipboard feedback.
- **Bulk File Queue Pipeline & Unified Quality Reports (`FileUpload.ts`)**:
  - Streamlined workflow (Dropzone $\to$ Queue $\to$ Transform $\to$ Cumulative Quality Report $\to$ Downloads $\to$ Bottom-Right Center on Map) removing legacy raw line preview.
- **Global Coordinate Order Promotion (`main.ts`, `components.css`)**:
  - Promoted `Coordinate order:` (`NE(H)` vs `EN(H)`) to the top Shared Configuration Bar, globally managing Point, Text, and Files tabs and synchronizing coordinate order across all components.
- **Code-Block Copy Action & Focus Auto-Select (`PointConverter.ts`, `components.css`)**:
  - Added top-right code-block style copy icon button above computed coordinates with tooltip and `✓ Copied` visual feedback.
  - Added auto-selection of full text on focus across all coordinate inputs.
- **Boundary-Safe Map Centering Across All Tabs (`MapView.ts`, `PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)**:
  - Guarded map centering and bounds fitting strictly within Romania's bounding box (`[43.5, 20.2]` to `[48.3, 29.8]`), preventing panning to oceans on out-of-grid coordinates.
  - Added dedicated `Center on Map` button in the diagnostic/report sections of Point, Text, and Files tabs.
- **Strict Numeric Validation & 1-to-1 Line Correspondence (`pointParser.ts`, `dmsParser.ts`, `BatchConverter.ts`)**:
  - Replaced lenient `parseFloat` parsing with strict numeric checks (`Number(token)` and regex `/^[+-]?(?:\d+\.?\d*|\.\d+)$/`) to immediately catch typos such as `5001a50.000` as syntax errors instead of partial coordinate evaluations.
  - Maintained strict 1-to-1 row correspondence between input lines and output lines in identical sequence.
  - Introduced human 1-based indexing (`Point 1`, `Point 2`, ...) for map pins and popups when input rows omit point IDs.
- **Bulk Files Queue Item Removal (`FileUpload.ts`)**:
  - Added individual `✕` remove button per queued file item in the drag-and-drop file queue.
- **Asymmetrical 3-Column Coordinate Layout (`components.css`, `PointConverter.ts`)**:
  - Allocated asymmetrical column ratios (`1.35fr 1.35fr 0.8fr`) for Northing/Lat, Easting/Lon, and Altitude/Height, preventing text clipping for long DMS angle strings (`46° 15' 27.88517"`).
- **Tab Prefix Label (`components.css`, `main.ts`, `en.ts`, `ro.ts`)**:
  - Added subtle prefix `Transform:` before folder bookmark navigation tabs (`Point`, `Text`, `Files`).
- **Selectable Format Example Badge (`components.css`, `main.ts`)**:
  - Live example badge displaying copy-ready formatted coordinates (e.g., `Example: 45° 50' 28.44420"`) with selectable text and copyable degree symbol.
- **Pre-Flight Input Validation (`PointConverter.ts`)**:
  - Instant validation when converting from Stereo 70 to ETRS89, checking for invalid angle characters (`°`, `'`, `"`) or letters and displaying inline feedback without making extraneous backend requests or triggering false "Out of grid" warnings.
- **Inline Diagnostic Reports (`PointConverter.ts`, `BatchConverter.ts`)**:
  - Integrated diagnostic cards reporting grid coverage status (`In grid coverage` / `Outside official grid envelope`), Quasigeoid undulation $\zeta$, detected delimiter/columns, and syntax error lines, eliminating disruptive modal popups.
- **Shared Transformation Configuration Bar (`main.ts`, `components.css`)**:
  - Animated pure-CSS sliding directional toggle with sleek 6px rounded-rectangle glider (`Stereo 70 ➔ ETRS89` $\leftrightarrow$ `ETRS89 ➔ Stereo 70`).
  - Native select dropdown for angular units with dynamic target/source labeling.
- **Geodetic Help Reference Modal (`GeodeticHelpModal.ts`)**:
  - Dedicated reference documentation for Romanian datums (Stereo 70, ETRS89, MN75, GRS80), Quasigeoid undulation ($\zeta = h - H$), and syntax conventions.
- **Header Direct Text Navigation (`Header.ts`)**:
  - Added text-only buttons for `Grid Info` and `Geodetic Help` directly into the navbar, removing `Stats & Grid` from input tabs.
- **Multi-File Queueing & Batch Streaming (`FileUpload.ts`)**:
  - Multi-file dropzone with manual `Transform Files` trigger, queue status tracking, per-file download links, and `Download All` batch export.
- **Auto 2D/3D Dimensionality Detection (`PointConverter.ts`)**:
  - Automatically runs 2D or 3D transformations based on presence of the vertical coordinate (Altitude/Height), eliminating the manual 2D/3D selector.

### Changed
- **Unspaced DMS Angle Formatting & Paste Handling (`dmsParser.ts`, `PointConverter.ts`)**:
  - Removed spaces after degree (`°`) and minute (`'`) characters (e.g. `46°00'00.00000"`), matching official Romgeo and PyTransDat conventions.
  - Added `isSingleDmsString` guard preventing single DMS coordinates from splitting across multiple input boxes on paste.
- **Shared Config Dynamic Angle Format Example Reactivity (`main.ts`)**:
  - Ensured dynamic reactivity of angle format example badge across Degrees (`45.84123450`), DMS (`45°50'28.44420"`), and Radians (`0.7997813210`).
- **Standardized Quality Report Warnings (`en.ts`, `ro.ts`, all tabs)**:
  - Standardized concise warning labels: `In Grid`, `Outside Grid`, `No Grid Data`, and `Syntax Error` (RO: `În grid`, `În afara gridului`, `Fără date în grid`, `Eroare sintaxă`).
- **Removed Redundant Input Diagnostic in Batch Text (`BatchConverter.ts`)**:
  - Removed `#batch-input-diagnostic` container, streamlining the vertical layout to input textarea $\to$ Transform button $\to$ results $\to$ Quality Report.
- **Dynamic Coordinate Order Swapping & Symmetry (`PointConverter.ts`)**:
  - Column 1 and Column 2 dynamically swap between Northing/Latitude and Easting/Longitude across inputs, outputs, and clipboard copy when toggling `NE` vs `EN`.
- **Geodetic Term Standardization & Localization (`en.ts`, `ro.ts`)**:
  - Standardized elevation difference metric to `Quasigeoid Height (ζ = h - H)` (RO: `Înălțimea cvasigeoidului (ζ = h - H)`).
  - Adopted sentence case `Angular unit:` and micro-headers `Input coordinates` / `Computed coordinates`.
- **Batch Out-of-Grid Coordinate Suppression (`BatchConverter.ts`)**:
  - When a point is outside the official TransDatRO grid envelope, approximate transformed coordinates are suppressed from the output text (emitting only `[Point_ID], # Out of grid`), preventing accidental import of unverified coordinates.
- **Batch Text Export Buttons Relocation (`BatchConverter.ts`)**:
  - Relocated `CSV`, `GeoJSON`, and `DXF` export actions below the output textarea grouped under label `Download As:`.
- **Symmetrical Read-Only Output Fields (`PointConverter.ts`)**:
  - Replaced plain text output row with symmetrical read-only form inputs matching input styling, column widths, and Title Case labels.

### Fixed
- **Radians Parsing & Roundtrip Conversion (`PointConverter.ts`, `BatchConverter.ts`, `pointParser.ts`)**:
  - Fixed exponential coordinate truncation bug when toggling angular units to/from Radians by correctly identifying previous radian state and converting via `radToDeg`.
  - Added native radian coordinate recognition to `PointConverter.calculate()` and `pointParser.parseInputLines()`, preventing false "Out of grid" failures and enabling `Center on Map` in radians mode.
- **Direction Switching Format Inversion (`PointConverter.ts`, `BatchConverter.ts`)**:
  - Replaced rigid string prefix checks (`startsWith('46.')`) with smart inversion using previous valid output coordinates or fallback to projection center origin, fixing corruption when converting non-center coordinates (e.g. Bucharest).
  - Synchronized default sample coordinates in Batch Text when switching direction.
- **Out-of-Grid Quasigeoid Height & Output Guard (`PointConverter.ts`)**:
  - Prevented false `+0.000` Quasigeoid Height calculation and masked coordinate outputs with persistent `—` when coordinates fall outside grid coverage.
- **Bulk File Error Row Detection False Positives (`FileUpload.ts`)**:
  - Fixed false-positive valid row counts by properly recognizing server `# Wrong dimension`, `# Invalid coordinate data`, and `# Limit exceeded` line errors.
- **Smart Map Click Out-of-Grid Guard & Center on Map Feedback (`main.ts`)**:
  - Guarded reverse Stereo 70 calculation against out-of-grid clicks on map edges, and provided explanatory toast warnings when centering on completely out-of-grid point sets.
- **CSV Export Out-of-Grid Bug (`exporters.ts`)**:
  - Ensured `exportToCSV` suppresses approximate coordinates for points flagged with `# Out of grid`, matching textarea output safety invariants and preventing unverified coordinates from being imported into CAD/GIS.
- **Diagnostic Cards & Delimiters Localization (`BatchConverter.ts`, `FileUpload.ts`, `en.ts`, `ro.ts`)**:
  - Replaced hardcoded English template literals with fully localized keys for detected points, columns, delimiters, in-grid/out-of-grid counts, syntax error summaries, and toast feedback.
- **Bidirectional Pre-Flight Validation (`PointConverter.ts`)**:
  - Enforced latitude `[-90°, 90°]` and longitude `[-180°, 180°]` checks on inverse `ETRS89ToStereo70` mode with safe exception handling, and prevented non-numeric 3D altitude inputs from silently evaluating as zero.
- **Dropzone Drag & Drop Flicker (`FileUpload.ts`, `components.css`)**:
  - Replaced inline style mutations and DOM rebuilding on dragleave with lightweight `.drag-over` CSS class toggling.

### Changed
- **Folder Bookmark Tab Styling (`components.css`, `main.ts`)**:
  - Replaced floating pill tabs with classic connected folder bookmark tabs seamlessly merged into the active card.
- **Workbench Panel Geometry (`base.css`)**:
  - Slimmed utility panel width from `640px` to `480px`, expanding the interactive Web GIS map viewport.
- **Symmetrical 3-Column Coordinate Layout (`PointConverter.ts`, `components.css`)**:
  - Symmetrical 3-column inputs and outputs for direct before-and-after visual alignment.
- **Standardized Geodetic Terminology (`en.ts`, `ro.ts`)**:
  - Standardized terms to Northing, Easting, Altitude / Latitude, Longitude, Height; removed verbose unit and datum suffixes.
  - Standardized primary action verb to `Transform` across all views with theme-aligned solid styling and no emojis.

### Changed
- **Navigation Tabs Modernization (`main.ts`, `components.css`)**:
  - Removed decorative icons/emojis from main tabs (`Punct Unic`, `Lot Coordonate`, `Fișier`, `Statistici`).
  - Switched to rounded pill/dock segmented design with elevated active states.
- **Segmented Option Controls (`components.css`)**:
  - Added subtle vertical divider bars between adjacent unselected options for clearer text parsing.
  - Switched active option pills from saturated blue to elevated neutral surfaces (`--bg-surface-elevated`) to restore clean visual hierarchy.
- **Geodetic Metrics Neutral Styling (`PointConverter.ts`)**:
  - Changed Quasigeoid undulation $\zeta$ value styling from green to neutral monospace typography (`var(--text-main)`), preventing confusion with success/validation states.
- **Workbench Panel Geometry (`base.css`)**:
  - Expanded utility panel width range to `560px`–`640px` to accommodate long coordinate tuples without wrapping.
- **Map Controls & Marker Clusters (`map.css`)**:
  - Styled Leaflet zoom controls with rounded corners and card theme backgrounds.
  - Updated cluster styling to match precision technical blue tokens.

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
