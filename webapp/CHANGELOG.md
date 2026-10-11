# Changelog - PyTransDatRO Web Application (`webapp/`)

All notable changes to the PyTransDatRO web application and its interactive Web GIS interface will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
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
