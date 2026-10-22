# PyTransDatRO Web Application UI Actionable Elements Inventory & Living Test Specification

This document provides a comprehensive, structured inventory of all actionable UI controls, reactive surfaces, floating panels, and event mechanics across the **PyTransDatRO** web application. It acts as the single source of truth bridging Product Requirements, Functional Interaction Specifications, and Automated E2E/Integration Testing (Playwright / Vitest).

> **Scope Note**: This document focuses strictly on the **functional mechanics and behaviors** of actionable elements. Visual styling, aesthetic tokens, and design decorations are documented separately in the webapp design guide.

---

## 1. Global Taxonomy & Interaction Classification

Every actionable element is categorized according to functional pattern, DOM selector, interaction type, scope, and test mapping:

| Interaction Pattern | Code / Control Type | Functional Mechanics |
| :--- | :--- | :--- |
| **Brand Action / Link** | Anchor `<a>` / Button `<button>` | Navigates to external link or resets view to Home while preserving theme/session state. |
| **Status Indicator / Alert** | Badge / Dot / Toast | Displays live status (e.g. API health, grid warning). May be clickable for diagnostics/feedback. |
| **Segmented Switch** | `renderSegmentedControl` | Mutually exclusive sliding/pill tabs (ARIA `tablist` / `tab`). |
| **Dropdown Menu / Select** | `renderDropdownMenu` / `renderSelect` | Expanding menu (`<details>` or `<select>`) with options, checkmarks, and keyboard navigation. |
| **Numeric / Text Input** | `<input type="text">` / `<textarea>` | Keyboard input, auto-selection on focus (`input.select()`), copy/paste parsing. |
| **Action CTA Button** | `.btn-primary` / `.btn-secondary` | Triggers asynchronous computation, file processing, or modal dialog. |
| **Map Widget / Control Dock** | Leaflet Custom Control (`L.Control`) | Interacts with the Web GIS canvas (zoom, pan, measure, layer drawer, scale bar). |
| **Floating Window / Modal** | Movable / Resizable Surface | Draggable header, corner resizer, minimize/restore, and close mechanics. |
| **Aggregator KPI Card** | Telemetry Summary Card | Reactive display surface aggregating totals, averages, and rates over active filters. |

---

## 2. Section I: Top Navigation Bar & Global State (`Header.ts`, `main.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Brand Logo** | Brand Action | `#header-brand-logo-wrap` | Click | App / Navigation | **Fresh First-Load Reset Action**: Clicking the logo icon completely resets all application runtime state: clears Point inputs and outputs, clears Batch text area and outputs, clears Uploaded Files queue, purges all map layers/pins, cancels active distance measurements, resets coordinate direction/order to default (`ETRS89 ➔ Stereo 70`, `NE(H)`), closes all open floating panels (Grid info, Help guide, Feedback), resets map view to Romania overview bounds, and returns to the `Transform` Studio tab in `Point` mode. Preserves user preferences (`data-theme` and active language `ro`/`en`). | Add E2E assertion verifying full state purge on logo click. |
| **App Title (`pytrans.ro`)** | Brand Action | `.header-brand-link .brand-title-text` | Click | App / Navigation | **Soft Home Navigation**: Returns view to application root `/` (Studio/Transform tab) via client SPA navigation without wiping uncommitted coordinate inputs or active working session data. | Covered in smoke specs. |
| **Version Badge** | Dynamic Badge | `.brand-stage-badge` | Read / Inspect | Header | Displays version string (e.g. `1.0.0-beta`). Populated dynamically from release metadata. | Validate dynamic version injection in build. |
| **API Health Indicator Dot** | Status Dot / Non-Disruptive | `#header-health-dot-btn` | Hover / Inspect | Header / Diagnostics | Compact 8px status dot (`.online` / `.offline`, green/red) situated stably next to the version badge. Green when REST API `/health` returns 200 OK (polled periodically every 30s); red when offline. Pure status indicator: does not expand header width, shift brand text, or displace navigation elements. Shows informative tooltip. | E2E error suite tests offline mock state. |
| **Global Offline Alert Banner** | System Alert Banner / Modal Trigger | `#global-offline-banner`, `#global-offline-report-btn` | Click | System Alert / Feedback | Prominent, non-disruptive notification banner docked directly beneath the header (or above the active viewport) when the REST API is offline. Displays explicit warning text that the REST API is unavailable and coordinate transformations cannot be computed, accompanied by an explicit actionable prompt: *“Click here to report this issue ↗”*. Clicking anywhere on the banner or action link launches the Feedback modal with `api_down` preselected. Auto-dismisses when the backend recovers. | Add E2E spec verifying offline banner visibility and report modal launch. |
| **Repository Attribution** | External Link | `.header-repo-link` | Click | External | Opens GitHub repository in a new tab (`target="_blank"`, `rel="noopener noreferrer"`). Label: `powered by PyTransDatRO`. | Test link target and rel attributes. |
| **Primary Tabs (`Transform` / `Usage`)** | Tablist | `#nav-btn-studio`, `#nav-btn-telemetry` | Click | Global App View | Switches between Geodetic Transformation Studio (`#workbench-layout` + `#map`) and Usage Telemetry Dashboard (`.telemetry-page-view`). Active tab receives `.active` class. Automatically triggers map size invalidation on return. | Covered in smoke/telemetry specs. |
| **Help Menu Trigger** | Dropdown Menu | `#header-help-details summary` | Click / Keyboard | Modal Triggers | Expands dropdown with two grouped sections: **PyTransDatRO** (Help Guide, Grid Info, Discrepancy Feedback) and **ANCPI / ROMPOS** external resource links. Selecting items launches their respective floating window (`GridMetadataPanel`, `GeodeticHelpPanel`) or feedback modal (`FeedbackPanel`, prefilled with current coordinates). Closes on click-outside or `Escape`. | Covered in UI unit & E2E tests. |
| **Language Selector** | Dropdown Menu / Radio Item | `#header-lang-details`, `[data-lang="ro"]`, `[data-lang="en"]` | Select option | Global (All Pages) | Toggles active locale (`ro` $\leftrightarrow$ `en`). Updates UI labels and re-renders active tabs. **Number formatting rule**: Usage/Telemetry dashboard numbers, cards, and labels follow the active locale (Romanian uses comma `,` decimal separator; English uses period `.`). In contrast, **all coordinate transformation inputs, outputs, angle format badges, and exported files (CSV/GeoJSON/DXF) strictly remain in standard numeric format with dot `.` decimal separator** to guarantee compatibility with CAD, GIS, and scientific software. | Covered in `i18n.spec.ts`. Check locale formatting in Usage statistics. |
| **Theme Toggle** | Action Button | `#header-theme-btn` | Click | Global UI Theme | Toggles `data-theme="dark"` $\leftrightarrow$ `data-theme="light"`. Swaps active icon (shows Sun when dark, Moon when light). Automatically switches Carto Dark $\leftrightarrow$ Carto Positron basemaps if active. Synchronizes browser favicon SVG. | Covered in `smoke.spec.ts`. |
| **Mobile View Switcher** | Responsive Button | `#header-mobile-view-btn` | Click | Viewport Layout | On mobile viewports (<768px), toggles visibility between the Workbench Panel and the full-screen GIS Map canvas. | Covered in smoke responsive tests. |
| **URL Query Parameters** | Deep-Link Hydration | `window.location.search` | Boot / Navigation | App Initialization | Reads `?op=`, `?c1=`, `?c2=`, `?c3=`, `?unit=`, and `?order=` on application load to pre-populate transformation state and calculate initial output. | Add dedicated URL deep-linking test spec. |
| **Global Escape Shortcut** | Keyboard Handler | `document (keydown 'Escape')` | Key press | Global View | Closes open header dropdowns (`details.header-dropdown-wrap`) and cancels active map distance measurement mode. | Add keyboard test assertion. |

---

## 3. Section II: Left Workbench — Transformation Controls (`PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)

### 3.1 Shared Transformation Configuration Bar (Top of Workbench)
*Globally governs Point, Text, and Files tabs.*

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Transformation Direction** | Segmented Switch | `#op-direction-toggle` / `[data-op]` | Click / Arrow keys | Global Transform Setting | Toggles between `ETRS89 ➔ Stereo 70` and `Stereo 70 ➔ ETRS89`. Updates input/output labels, placeholder values, and validation boundaries. **Deferred processing rule**: Switching direction applies automatically to the **Point tab** only. For **Text** and **Files** tabs, it does **not** automatically recompute existing input data; the user must hit `Transform` to re-process with the new direction. | Covered in `point.spec.ts`, `batch.spec.ts`. |
| **Angular Unit Selector** | Dropdown Select | `#angular-unit-select` | Select option | Global Coordinate Format | Options: `Degrees (DD.dddddddd)`, `DMS (DD°MM'SS.SSSSS")`, `Radians (rad)`. Formats ETRS89 angles in inputs, outputs, angle format badge, and map coordinate cursor badge. **Deferred processing rule**: Applies immediately to the **Point tab**; for **Text** and **Files** tabs, the format is applied to output data only after the user hits `Transform`. | Covered in `point.spec.ts`. |
| **Coordinate Order Toggle** | Segmented Switch | `#coord-order-toggle` / `[data-order]` | Click | Global Coordinate Order | Options: `NE(H)` (Northing/Lat first) vs `EN(H)` (Easting/Lon first). Swaps input and output column positions in Point view and determines column parsing order in Batch Text & Files. **Deferred processing rule**: Changes column layout immediately in Point tab, but in Text and Files tabs applies to output data only after hitting `Transform`. | Covered in `point.spec.ts`, `batch.spec.ts`. |
| **Angle Format Example Badge** | Clickable Badge | `#angle-format-example-badge` | Click | Clipboard / Selection | Displays an example coordinate formatted in the active angular unit (e.g., `45°50'28.44420"`). **Click action**: Selects all text and copies strictly the numeric angle value to clipboard for quick paste into inputs. Shows temporary checkmark feedback. | Covered in `ui.test.ts`. Verify E2E clipboard action. |

---

### 3.2 Navigation Mode Tabs
| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Transform Mode Tabs** | Tablist | `.tabs-bar .tab-btn` (`#tab-point`, `#tab-text`, `#tab-files`) | Click | Workbench Surface | Switches active mode panel between `Point`, `Text`, and `Files`. Prefix label: `Transform:`. | Covered in `smoke.spec.ts`. |

---

### 3.3 Single Point Mode (`PointConverter.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Input Coordinate 1** | Text Input | `#pt-input-1` | Keyboard / Paste | Point State | Northing (Stereo 70) or Latitude (ETRS89). Auto-selects text on focus (`input.select()`). Auto-sanitizes input (strips quotes, trailing `m/metri`, compass letters). **Smart Delimiter Split**: Pasting a full 2D or 3D coordinate line (separated by spaces, tabs, commas, or semicolons) automatically splits across Input 1, Input 2, and Input 3. **Smart paste convenience rule**: Pasting a full coordinate line automatically triggers calculation for immediate feedback. Manual keyboard entry continues to require pressing `Enter` or clicking `Transform`. Pasting a single DMS string does not split. Coordinates do not show units of measure, except DMS which explicitly shows sexagesimal symbols (`°`, `'`, `"`). Default origin initialized if empty. **No Inline Error Banners**: All inline validation error banners below coordinate boxes are completely eliminated. Error and quality feedback is communicated exclusively through the structured Quality Reporting Card. | Covered in `point.spec.ts`. |
| **Input Coordinate 2** | Text Input | `#pt-input-2` | Keyboard / Paste | Point State | Easting (Stereo 70) or Longitude (ETRS89). Auto-selects text on focus. Auto-sanitizes input. Units of measure are not shown, except for DMS symbols (`°`, `'`, `"`). Protected by identical error banner elimination rule. | Covered in `point.spec.ts`. |
| **Input Altitude/Height (3rd Coord)** | Text Input | `#pt-input-3` | Keyboard / Paste | Point State | Normal Altitude $H$ or Ellipsoidal Height $h$. Optional: if empty, 2D transformation; if provided, 3D transformation with Quasigeoid height. Auto-selects on focus. Protected by identical error banner elimination rule. | Covered in `point.spec.ts`. |
| **Input Enter Key Listener** | Keyboard Shortcut | `#pt-input-1`, `#pt-input-2`, `#pt-input-3` | Enter key press | Point Execution | Pressing `Enter` in any of the coordinate text inputs executes the transformation immediately without needing to click the button. | Add E2E assertion for Enter key. |
| **Transform CTA Button** | Action Button (`.btn-primary`) | `#pt-calculate-btn` | Click / Enter | Point Output | Executes backend transformation via `/api/v1/transform/point`. Populates output fields and emits map pin. Shows calculating spinner/label during transit. | Covered in `point.spec.ts`. |
| **Output Coordinate Fields** | Read-Only Fields | `.output-grid-3 input` | Select / Copy | Point Output | Displays computed coordinates. **Selectable but not editable**: Users can highlight and copy text. Displays persistent em-dash `—` when cleared or untransformable. No units of measure shown, except DMS symbols (`°`, `'`, `"`). | Covered in `point.spec.ts`. |
| **Quasigeoid Height ($\zeta$)** | Scalar Metric Box | `.metric-field` | Select / Copy | Point Output | Read-only, selectable field showing $\zeta = h - H$ for 3D points. Stripped of unit suffixes. Prefixed with `+` if positive. Persistent em-dash `—` when 2D or invalid. | Covered in `point.spec.ts`. |
| **Reporting & Quality Flags** | Visual Labels / Badges | `.quality-report-badges` | Read-Only (Non-selectable) | Quality Feedback | Visual artifact labels (not selectable) structured in a **Two-Section Card** separated by a discrete **dashed divider line** with subtle sub-titles: **Top Subsection ("Coordinates" / "Coordonate")**: displays dimensional tags (`2D`, `3D`), which read contextually as `[Section] + [Label]` = *Coordinates 2D*, *Coordinates 3D*. **Bottom Subsection ("Transformation" / "Transformare")**: displays computation quality tags (`Valid`, `Invalid Input`, `Invalid Angle`, `Outside Grid`, `No Grid Data`), which read contextually as *Transformation Valid*, *Transformation Outside Grid*, etc. **The Golden Rule**: The UI strictly exposes only classifications understood by the engine; no arbitrary frontend error codes or interval bounds. | Covered in `point.spec.ts`. |
| **Copy Coordinates Button** | Copy Action Button | `#pt-copy-code-btn` | Click | Clipboard | Copies computed output coordinates to clipboard formatted as **comma-separated values** (`val1, val2, [val3]`). Shows feedback checkmark animation and toast notification. Disabled when out-of-grid or uncalculated. | Covered in `point.spec.ts`. |
| **Center on Map** | Action Button | `#pt-center-map-btn` | Click | GIS Map Viewport | Pans map smoothly to the centroid of the point without changing current zoom level. Guarded by Romanian grid bounds. Disabled when out of grid. | Covered in `map.spec.ts`. |
| **Zoom to Extents** | Action Button | `#pt-zoom-extents-btn` | Click | GIS Map Viewport | Fits map to point bounds with a maximum single-point zoom (level 14) so geographic context is maintained. Guarded by Romanian bounds. Disabled when out of grid. | Covered in `map.spec.ts`. |

---

### 3.4 Batch Text Mode (`BatchConverter.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Batch Input Text Area** | Textarea Input | `#batch-input-area` | Paste / Edit | Batch State | Accepts freeform multi-line coordinate text. **Ergonomics**: Vertically resizable to view large point lists, with vertical scrollbar styled with `cursor: default` to guarantee standard arrow hover cursor. **Height Retention Rule**: Clicking `Transform` strictly preserves the user-adjusted vertical height and dimensions of both input and output textareas without collapsing back to initial default sizes. **Mixed Input Parsing**: Must support mixed 2D and 3D points in the same batch, and points with or without point names (3 columns are parsed as 3D coordinates unless the 1st column contains alphabetic/text tokens, in which case it is identified as a Point Name). **Delimiter Detection**: Auto-detects line delimiters (comma, semicolon, tab, space). **Comment & Blank Preservation**: Lines beginning with `#` or blank lines maintain their exact sequence index in the output. Line $N$ in corresponds to Line $N$ out. | Covered in `batch.spec.ts`. |
| **Batch Transform CTA** | Action Button (`.btn-primary`) | `#batch-transform-btn` | Click | Batch Processing | **Input Locking**: Locks/disables both the input textarea (`disabled readonly`) and the transform button during active computation to enforce input/output synchronization. Coordinate mathematics executes swiftly without pre-transformation modal friction; if point volume plotted to Leaflet exceeds the threshold (`batchWarningThreshold` in `webapp/src/config.ts`), `MapView` triggers the Volume Warning dialog before canvas rendering. | Covered in `batch.spec.ts`. |
| **Computed Results Text Area** | Read-Only Textarea | `#batch-output-area` | Select / Copy | Batch Output | Non-editable, selectable output textarea showing transformed coordinates, preserved comments, or standardized quality annotations. **Clean Output Rule for Invalid Lines**: When input rows have syntax errors (`Invalid input`) or physically impossible angles (`Invalid angle`), invalid coordinate numbers are strictly omitted from the output. The output row contains ONLY the Point Name (if provided) followed by the standardized quality tag (e.g. `104, # Invalid input` or `101, # Invalid angle`). If no point name was given, only the quality tag is output (e.g. `# Invalid input`). Symmetrical like-for-like row mapping guarantees line $N$ in corresponds to line $N$ out. | Covered in `batch.spec.ts`. |
| **Copy Results Button** | Copy Action Button | `#batch-copy-results-btn` | Click | Clipboard | Copies entire computed output text to clipboard. Shows feedback tooltip and toast. Disabled when empty. | Covered in `batch.spec.ts`. |
| **Quality Report Summary** | Visual Labels / Badges | `.quality-report-card`, `.quality-report-badges` | Read-Only (Non-selectable) | Quality Feedback | Visual artifact labels (not selectable) arranged in a **Two-Section Card** separated by a discrete **dashed divider line** with subtle sub-titles: **Top Subsection ("Coordinates" / "Coordonate")**: displays volume and breakdown counts (`Total Points: [N]`, `2D: [N]`, `3D: [N]`), reading contextually as *Coordinates 2D (N)*, *Coordinates 3D (N)*. **Bottom Subsection ("Transformation" / "Transformare")**: displays diagnostic counts for calculation results (`Valid: [N]`, `Invalid Input: [N]`, `Invalid Angle: [N]`, `Outside Grid: [N]`, `No Grid Data: [N]`), reading contextually as *Transformation Valid (N)*, *Transformation Outside Grid (N)*, etc. | Covered in `batch.spec.ts`. |
| **Download CSV** | Action Button | `#batch-export-csv` | Click | File Download | Exports exact line-by-line replica of computed results as `.csv` file (including comment and quality check rows). **No Premature Toast**: Clicking export initiates the browser's native Save As file prompt without firing a premature "File Downloaded" notification toast. | Covered in `batch.spec.ts`. |
| **Download GeoJSON** | Action Button | `#batch-export-geojson` | Click | File Download | Exports valid in-grid points as GeoJSON FeatureCollection. **Clean filtering**: Filters out all invalid outputs (out of grid, no data on grid, syntax errors, invalid angles) to avoid invalid geometries in GIS software. **No Premature Toast**: Directly triggers browser Save As dialog without premature notification. | Covered in `batch.spec.ts`. |
| **Download DXF** | Action Button | `#batch-export-dxf` | Click | File Download | Generates standard ASCII AutoCAD DXF file containing POINT entities for valid coordinates. **Clean filtering**: Filters out all invalid outputs (out of grid, no data on grid, syntax errors, invalid angles) to prevent CAD entity corruption. **No Premature Toast**: Directly triggers browser Save As dialog without premature notification. | Covered in `batch.spec.ts`. |
| **Center / Zoom Extents** | Action Buttons | `#batch-center-map-btn`, `#batch-zoom-extents-btn` | Click | GIS Map Viewport | Pans or fits bounds (max zoom 14) to the batch text point cluster. Disabled when point set is empty. | Covered in `map.spec.ts`. |

---

### 3.5 File Upload Mode (`FileUpload.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **File Drag & Drop Zone** | Dropzone / File Picker | `#file-dropzone`, `#file-input-hidden` | Drag & Drop / Click | File Queue | Accepts `.txt`, `.csv`, `.xyz`, `.dat` files up to 20MB per file. Displays drag-over visual indicator. **Staged Workflow**: Supports staging multiple files. Allows adding files, transforming them, changing global transformation settings (e.g. coordinate order or direction), staging additional files, and transforming the untransformed files using the updated app settings. | Covered in `file.spec.ts`. |
| **Transform Files CTA** | Action Button (`.btn-primary`) | `#file-start-btn` | Click | File Queue Processing | Processes untransformed files in the queue sequentially using current active settings. Triggers volume warning modal if cumulative points exceed 1,000 pts. Displays streaming progress. | Covered in `file.spec.ts`. |
| **Map Visualization of Files** | GIS Rendering | `filesLayer` | Auto / Viewport | GIS Map Canvas | Map displays markers and clusters cumulatively from **all processed files** in the queue. | Covered in `file.spec.ts`. |
| **Volume Warning Dialog** | Modal Alert Dialog | `.modal-backdrop`, `#vw-sample-btn`, `#vw-all-btn`, `#vw-skip-btn` | Confirm Choice | Map Rendering Guard | Triggered when point count exceeds configured threshold (`batchWarningThreshold`, default 2,000). User chooses: `Render Sample` (sampled at `batchSampleSize`, default 1,000), `Render All Points`, or `Skip Map Display` to maintain browser performance. Close button defaults to skip. Configured centrally in `webapp/src/config.ts`. | Covered in `file.spec.ts`. |
| **Individual File Download** | File Action Button | `.file-download-btn[data-id]` | Click | File Download | Downloads individual transformed file. Directly invokes browser save dialog without premature notification toasts. | Covered in `file.spec.ts`. |
| **Individual File Remove** | Icon Action Button | `.file-queue-remove[data-id]` | Click | File Queue | Removes individual file item from queue and removes its associated points from the map layer. | Covered in `file.spec.ts`. |
| **Download All Button** | Action Button | `#file-download-all-btn` | Click | Bulk Download | Sequentially triggers downloads for all completed files in the queue (staggered by 300ms to avoid browser download blocking). Visible when more than 1 file is completed. Directly invokes browser save dialogs without premature notification toasts. | Add test for bulk download action. |
| **Clear Queue Button** | Action Button | `#file-clear-queue-btn` | Click | File Queue | Clears all uploaded files from queue and completely purges the files map layer. | Covered in `file.spec.ts`. |
| **File Queue Quality Summary** | Visual Labels / Badges | `.file-queue-stats`, `.quality-report-card` | Read-Only | Queue Diagnostics | Displays aggregate quality report for processed files using the **Two-Section Card** layout: Top: Coordinates (Total Points, 2D, 3D); Bottom: Transformation (Valid, Invalid Input, Invalid Angle, Outside Grid, No Grid Data) separated by dashed line. | Add file quality report verification test. |

---

## 4. Section III: Web GIS Map Canvas (`MapView.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Map Click (Smart Coordinate Picker)** | Canvas Interaction | `#map` | Click on canvas | Point Mode Inputs | In `ETRS89ToStereo70`: Populates clicked lat/lon in active angular units into Point inputs. In `Stereo70ToETRS89`: Reverse calculates Stereo 70 Northing/Easting from map click. **Automatic trigger**: Map click is the only action that automatically computes output coordinates and updates map pin without requiring the Transform button. Only active when Point tab is selected. **Untransformable Point Click Erase Rule**: When clicking in an untransformable location (outside grid or no grid data), **both input and output coordinates are immediately erased/cleared out** (`""` in inputs and `—` in outputs), and the quality report displays `Outside Grid` or `No Grid Data`. This guarantees the user is never left with misleading coordinates from previous or partial calculations. Applies symmetrically to both `ETRS89 ➔ Stereo 70` and `Stereo 70 ➔ ETRS89` directions. | Covered in `map.spec.ts`. |
| **Map Marker Symbology** | Visual Marker | `.custom-point-marker`, `.custom-cluster-circle` | Visual / Click popup | Map Canvas | Point markers strictly retain their assigned layer color from the Layers drawer (**never turn red, flash, or change color on selection or warning/error**). Errors and warnings are reported strictly via the workbench reporting area (`.quality-report-card`). The single point marker features a continuous soft pulsing ripple animation (`markerPulse`) around it at all times to distinguish it from batch points. Clusters group points and inherit the layer color. Center on Map and Zoom to Extents buttons allow pinpointing location. | Covered in `map.spec.ts`. |
| **CAD Point Label Dynamic Visibility** | Canvas Rendering | `.zoom-labels-visible` | Zoom Map | Marker Labels | Point ID text labels on map markers automatically hide below zoom level 15 to prevent clutter, and display at zoom $\ge 15$. | Add zoom threshold visual test. |
| **Cluster Spiderfy** | Cluster Behavior | `.custom-marker-cluster` | Click cluster at max zoom | Map Markers | Marker clusters spiderfy (fan out in a circular spiral) when clicked at maximum zoom levels so identical/overlapping points can be individually inspected. | Covered in `map.spec.ts`. |
| **Map Point Attribute Popup** | Info Popup / Clipboard Action | `.geodetic-popup`, `.geodetic-popup-copy-btn` | Click marker on map | Popover Readout | Displays comprehensive key-value attribute table showing: Point ID, coordinates in both ETRS89 (formatted in angular units active at computation time) and Stereo 70 (meters), Altitude/Height ($H$ and $h$), Quasigeoid height $\zeta$ (in cyan), and quality status (e.g. Valid in grid, Out of grid). Contains a **Copy Attributes** button (`Copiază atributele`) to copy the key-value table to clipboard with in-place feedback. | Covered in `map.spec.ts`. |
| **Zoom In / Zoom Out** | Map Dock Control | `#map-zoom-in-btn`, `#map-zoom-out-btn` | Click | Map Zoom Level | Adjusts map zoom by $\pm 1$ step. | Leaflet custom control dock. |
| **Home (Zoom to Romania)** | Map Dock Control | `#map-home-btn` | Click | Map Viewport | Centers map on Romania bounds (`[43.5, 20.2]` to `[48.3, 29.8]`) at an ergonomic overview zoom (level 7). | Covered in `map.spec.ts`. |
| **Measure Tool (Single Distance)** | Map Dock Toggle | `#map-measure-btn` | Click to activate $\to$ Click Pt 1 $\to$ Pt 2 | Measurement Layer | Enables measurement cursor: clicks 2 points on map, draws polyline, displays geodesic distance in meters/km with tooltip. Measures strictly **one single distance at a time**; starting a third click resets and clears previous measurement. Pressing `Escape` cancels and clears measurement. | Covered in `map.spec.ts`. |
| **Fullscreen Toggle** | Map Dock Toggle | `#map-fullscreen-btn` | Click | Browser Screen | Enters/exits browser fullscreen mode; synchronizes icon state. Automatically invalidates map dimensions. | Covered in `map.spec.ts`. |
| **Legend Widget** | Expandable Map Widget | `.map-legend-dock`, `#map-legend-header` | Click header | Legend Visibility | Shows active transformation layers: Single Point, Batch Text, Uploaded Files. **Color synchronization**: Legend swatches strictly honor the colors selected in the Layers drawer. Can be minimized/expanded (persisted in `localStorage['pytransdat_legend_open']`). Reflows dynamically above graphic scale bar. | Covered in `map.spec.ts`. |
| **Graphic Scale Bar** | GIS Graphic Scale | `.map-scale-dock` | Read-Only | Map Canvas | High-contrast segmented graphic scale bar showing distance in meters or kilometers responsive to zoom. | Covered in `map.spec.ts`. |
| **Layers Drawer Button** | Drawer Toggle Button | `#map-layers-toggle` | Click | Layers Menu | Icon button in bottom-left opening drawer for Base Maps, Shift Grids (Coming Soon), and Data Layers. Closes on click-outside. | Covered in `map.spec.ts`. |
| **Layer Visibility Checkboxes** | Checkboxes | `#layer-chk-point`, `#layer-chk-text`, `#layer-chk-files` | Check / Uncheck | Layer Visibility | Independently toggles visibility of Point, Text, and File layers on the map. Synchronizes layer count badge and Legend. | Covered in `map.spec.ts`. |
| **Base Map Radio Selector** | Radio Group | `input[name="basemap-choice"]` | Click option | Leaflet TileLayer | Mutually exclusive basemap selection: OpenStreetMap (default), Carto Positron, Carto Dark Matter, Satellite. Persisted in `localStorage['pytransdat_basemap']`. Coupled to theme for Carto options. | Covered in `map.spec.ts`. |
| **Data Layer Color Swatch Pickers** | Color Swatch Popover | `.layer-color-trigger`, `.color-swatch-option` | Click | Layer Symbology | **Distinguishable defaults**: Default layer colors are easily distinguishable and avoid red by default (red reserved for selectable custom choice). Changing swatch updates markers, clusters, and legend swatches immediately. Persisted in `localStorage['pytransdat_map_colors']`. | Covered in `map.spec.ts`. |
| **Coordinate Readout Badge** | Status Readout | `.map-coords-badge` | Mousemove on map | Bottom Coordinate Bar | Shows cursor coordinates dynamically formatted in active angular units (Degrees, DMS, Radians) and coordinate order (`Lat/Long` in `NE` mode vs `Long/Lat` in `EN` mode). | Covered in `map.spec.ts`. |

---

## 5. Section IV: Floating Panels & Diagnostic Windows

### 5.1 Grid Metadata Panel (`GridMetadataPanel.ts`)
*Accessed via Header Help Menu $\to$ `Grid Info` (`#menu-grid-info-btn`).*

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Draggable Header Handle** | Drag Handle | `#fgp-drag-handle` | Pointer drag | Panel Position | Pointer-captured dragging allowing moving the panel across the viewport with boundary clamping. | Covered in panel tests. |
| **Window Minimize / Restore** | Window Control Button | `#fgp-minimize-btn` | Click | Window Geometry | Toggles panel between full content view and collapsed title bar (`—` $\leftrightarrow$ `□`). | Covered. |
| **Window Close Button** | Window Control Button | `#fgp-close-btn` | Click | Panel Lifecycle | Closes and unmounts the floating grid metadata window. | Covered. |
| **Corner Resize Handle** | Window Resizer | `#fgp-resize-handle` | Pointer drag | Window Dimensions | Pointer-captured resizing from bottom-right corner with minimum (360x280px) and viewport maximum clamping. | Covered. |
| **Metadata Search Filter** | Search Input | `#fgp-search`, `#fgp-clear-search` | Input / Clear | Metadata Tree | Real-time search input filtering JSON metadata keys and values; clear button `✕` resets search query. | Add search filter assertion. |
| **Expand / Collapse All** | Action Buttons | `#fgp-expand-all-btn`, `#fgp-collapse-all-btn` | Click | Tree State | Recursively expands or collapses all branches in the hierarchical grid metadata tree. | Covered in unit tests. |
| **Copy Raw JSON** | Action Button | `#fgp-copy-json-btn` | Click | Clipboard | Copies complete raw grid metadata JSON structure to clipboard with toast notification. | Add test for JSON clipboard copy. |
| **Tree Branch Expander** | Tree Node Expander | `.fgp-expander` | Click | Branch Node | Toggles individual tree branch open/closed (`▸` $\leftrightarrow$ `▾`). | Covered in unit tests. |

---

### 5.2 Geodetic Help Panel (`GeodeticHelpPanel.ts`)
*Accessed via Header Help Menu $\to$ `Help Guide` (`#menu-help-app-btn`).*

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Draggable / Resizable Window** | Window Mechanics | `#fhp-drag-handle`, `#fhp-minimize-btn`, `#fhp-close-btn` | Drag / Minimize / Close | Panel Geometry | Supports movable header, minimize to title bar, and close actions. | Covered. |
| **Dynamic Benchmark Comparison Card** | Dynamic Surface | `#fhp-benchmark-slot` | Async Load / View | Geodetic Verification | Fetches and renders live mathematical comparison against official ANCPI TransDatRO vectors from `/api/v1/grid/comparison`. | Add test verifying live comparison display. |
| **External ANCPI / ROMPOS Links** | External Hyperlinks | `.fhp-body a` | Click | External Browser Tab | Opens official ANCPI geodetic documentation and ROMPOS services (`target="_blank"`). | Covered. |

---

### 5.3 Feedback & Discrepancy Reporting Modal (`FeedbackPanel.ts`)
*Accessed via Header Help Menu $\to$ `Report Issue` (`#menu-feedback-btn`) or clicking Offline API Health Dot.*

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Report Type Selector** | Dropdown Select | `#fb-type-select` | Select option | Feedback Template | Options: `discrepancy`, `grid_update`, `bug`, `feature`, `api_down`. Dynamically loads contextual pre-filled issue template into the message textarea. | Covered in feedback tests. |
| **Message Textarea** | Textarea Input | `#fb-message-input` | Edit / View | Feedback Payload | Accepts detailed issue description. Pre-populates with client URL, timestamp, active grid file, and environment details. | Covered. |
| **Reset Template Button** | Action Button | `#fb-reset-template-btn` | Click | Template State | Reverts textarea content back to the default pre-filled template for the currently selected report type. | Covered. |
| **Copy Report Button** | Action Button | `#fb-copy-btn` | Click | Clipboard | Copies formatted report text to clipboard with toast notification. | Covered. |
| **Open as GitHub Issue** | Action Button | `#fb-github-btn` | Click | External | Builds pre-filled GitHub issue URL with title and body, opening in new tab (`target="_blank"`). | Covered. |
| **Human Math Verification** | Bot Prevention Input | `#fb-math-input` | Numeric input | Submission Guard | Displays lightweight human challenge (`X + Y = ?`). Required before submitting. Prevents headless automated spam. | Covered. |
| **Submit Feedback Button** | Action Button (`.btn-primary`) | `#fb-submit-btn` | Click | API Submission | Submits report payload to `/api/v1/feedback`. Validates math answer, minimum elapsed time (>1.2s), and honeypot field. Displays progress and success toast. | Covered in feedback tests. |
| **Cancel / Close Buttons** | Window Actions | `#ffp-close-btn`, `#fb-cancel-btn` | Click | Panel Lifecycle | Closes and unmounts the feedback modal. | Covered. |

---

## 6. Section V: Usage & Telemetry Dashboard (`TelemetryDashboard.ts`)

*Accessed via `Usage` header tab. Contains 3 sub-views: Computation, Spatial, Legacy.*

### 6.1 Sub-Navigation Tabs
| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Usage Sub-Tabs** | Tablist Buttons | `#tab-timeline-btn`, `#tab-spatial-btn`, `#tab-legacy-btn` | Click | Telemetry Dashboard View | Switches active sub-page between Computation Timeline, Spatial Grid Matrix, and Legacy Historical Data. Active tab receives `.active` and underline indicator. | Covered in `legacy.spec.ts`. |

---

### 6.2 Sub-View 1: Computation (`TimelineChart.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Aggregator Summary Cards** | KPI Summary Cards | `.telemetry-kpi-banner .telemetry-kpi-card` | Read-Only / Filter Responsive | Computation View | Displays aggregated totals (Total Points, Total Calls, Batch Intensity ratio). Dynamically recalculated from active selected sources. Formats numbers according to active locale (comma `,` in Romanian; period `.` in English). | Covered in `legacy.spec.ts`. |
| **Time Window Presets** | Segmented Control | `[data-preset]` | Click | Time Filtering | Options: `All Time`, `12 Months`, `Year-to-Date (YTD)`, `30 Days`, `7 Days`. Selecting a preset automatically updates the date range input pickers to reflect the start and end dates of that window. Filters data reactively. | Covered in `legacy.spec.ts`. |
| **Granularity Selector** | Segmented Control | `[data-granularity]` | Click | Chart Grouping | Options: `Daily`, `Weekly`, `Monthly`. Groups time-series histogram bars accordingly. | Covered in `legacy.spec.ts`. |
| **Date Range Inputs** | Date Pickers | `#comp-date-start`, `#comp-date-end` | Date change | Custom Date Filter | Start date defaults to earliest DB telemetry date (`2023-11-01`); End date defaults to current date. Strict clamping prevents picking dates outside DB boundaries. Syncs when presets are clicked. | Covered in `legacy.spec.ts`. |
| **Source Filter Toggles** | Toggle Buttons | `[data-source-name]` | Click (Multi-select) | Source Filtering | Independent multi-toggle buttons with status dot indicators for sources: `REST API`, `Web App`, `Desktop / CLI`. Unknown sources filtered out. | Covered in `legacy.spec.ts`. |
| **Metric Filter Toggles** | Toggle Switches | `#comp-metric-points-btn`, `#comp-metric-calls-btn` | Click (Multi-select) | Chart Metrics | Allows displaying `Points`, `Calls`, or **both simultaneously** on the computation timeline SVG chart. | Add dedicated E2E assert for dual metric toggle. |
| **Computation Total Overlay Toggle** | Checkbox / Toggle Button | `#comp-toggle-total` | Click | Chart Overlay | Toggles cumulative total line across active computation series, mirroring the `#legacy-toggle-total` control in the Legacy tab. | Covered in computation specs. |
| **Chart Series Color Pickers** | Color Swatch Popovers | `.interactive-legend-item`, `.legend-swatch-option` | Click / Swatch select | Chart Series Symbology | Allows user to customize and select colors for individual series lines/bars and Total line in the timeline chart. Closes on click-outside. | Covered in `legacy.spec.ts`. |

---

### 6.3 Sub-View 2: Spatial (`SpatialMatrixChart.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Spatial KPI Cards** | KPI Summary Cards | `.spatial-overview-banner` | Read-Only | Spatial View | Displays spatial aggregations (Total Metric Value, Active Cells Count). **Filter synchronization**: Strictly synchronizes with **all applied filters** (time window, custom date range, sources). Formats numbers according to active locale. | Covered. |
| **Quadrant Share Cards** | Mathematical KPI Cards | `.quadrant-card` (NW, NE, SW, SE) | Read-Only | Regional Distribution | Displays metric value and percentage share across Romania's 4 geographic quadrants divided by projection center origin. | Covered. |
| **Spatial Metric Selector** | Mutually Exclusive Radio Selector | `#spatial-metric-points-btn`, `#spatial-metric-calls-btn` | Click (Single Choice) | Heatmap Density | Unlike computation view, this is a **strict radio selector**: selects either `Points` or `Calls` so coordinate density values do not combine into meaningless sums. | Add E2E check verifying mutual exclusivity. |
| **Spatial Scale Mode Toggle** | Segmented Switch | `#spatial-scale-log-btn`, `#spatial-scale-linear-btn` | Click | Color Scaling | Toggles between `Logarithmic` and `Linear` heatmap bin color scaling for dense clusters. | Covered. |
| **Spatial Grid Matrix Canvas** | Interactive Canvas | `#spatial-matrix-slot` | Pan / Hover | Romanian Coordinate Grid | Interactive SVG/Canvas rendering transformation density binned into 10km geodetic cells over Romania's bounds. Hover displays cell count tooltip. | Covered. |
| **5-Step Color Ramp Legend** | Graphic Legend | `.spatial-legend-bar` | Read-Only | Density Scale | Displays discrete 5-step color scale: Empty, Low, Medium, High, Peak (labeled with dynamic max cell hits). | Covered. |

---

### 6.4 Sub-View 3: Legacy Historical Data (`LegacyTimelineChart.ts`)

| ID / Element Name | Control Type | DOM Locator / Test ID | User Action | Target Scope | Expected Behavior & Side Effects | Test Coverage Gap / Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Legacy Aggregator Cards** | KPI Summary Cards | `.telemetry-kpi-banner .telemetry-kpi-card` | Read-Only | Legacy View | Displays historical transformation statistics from legacy Java TransDatOnline service (2011–2026): Total, App, Service, Stereo 70, Stereo 30 points. Reacts to legacy filter selections. Formats numbers according to active locale. | Covered in `legacy.spec.ts`. |
| **Legacy Time & Granularity Presets** | Segmented Controls | `[data-legacy-preset]`, `[data-legacy-granularity]` | Click | Legacy Filters | Preset time windows (`All`, `12M`, `YTD`, `30D`, `7D`) and temporal bucket grouping (`Daily`, `Weekly`, `Monthly`). Presets update date pickers based on DB maximum history date (`2026-10-14`), not system clock. | Covered in `legacy.spec.ts`. |
| **Legacy Date Range Inputs** | Date Pickers | `#legacy-date-start`, `#legacy-date-end` | Date change | Archive Date Filter | Fixed archive boundaries clamped between `2011-01-09` and `2026-10-14`. | Covered in `legacy.spec.ts`. |
| **Legacy Source Filter** | Toggle Buttons | `[data-legacy-source="app"]`, `[data-legacy-source="service"]` | Click (Multi-select) | Legacy Sources | Limited strictly to legacy traffic sources: `Web App` and `REST API (Service)`. At least one source must remain active. | Covered in `legacy.spec.ts`. |
| **Legacy Operation Filter** | Toggle Buttons | `[data-legacy-op="st70"]`, `[data-legacy-op="st30"]` | Click (Multi-select) | Legacy Operations | Filters historical operations: `Stereo 70 ➔ ETRS89` and `Stereo 30 ➔ ETRS89` (historical regional projection). At least one operation must remain active. | Covered in `legacy.spec.ts`. |
| **Legacy Total Overlay Toggle** | Checkbox / Toggle | `#legacy-toggle-total` | Click | Chart Overlay | Toggles cumulative total line across active legacy series. | Covered in `legacy.spec.ts`. |

---

## 7. Coverage Gap Analysis & Recommended Action Items

From our analysis of existing Playwright tests in `webapp/e2e/` against the comprehensive inventory, here are the specific gaps to address when generating or updating tests:

1. **Locale-Aware Number Formatting (Usage Dashboard vs. Transformation Purity)**:
   - *Requirement*: Usage / Telemetry dashboard statistics format numbers with commas (`,`) in Romanian and periods (`.`) in English. Transformation inputs, outputs, clipboard operations, and export files strictly remain standard dot (`.`) decimals to preserve CAD/GIS interoperability.
   - *Action Item*: Add assertions in `i18n.spec.ts` and `legacy.spec.ts` verifying that Usage statistics reflect locale-aware decimal formatting while coordinate boxes remain dot-decimal.
2. **Single Distance Measurement Reset & Esc Key**:
   - *Requirement*: Measure tool must measure strictly one line at a time, reset on subsequent clicks, and cancel cleanly on `Escape`.
   - *Action Item*: In `map.spec.ts`, verify that starting a third click clears the previous measurement line/tooltip, and verify `Escape` deactivates the tool.
3. **Marker Color Invariance on Map Selection & Continuous Pulsing**:
   - *Requirement*: Point markers on the map must permanently stay the assigned layer color from the Layers drawer (never turn red, flash, or change color on selection or warning/error). Single point marker exhibits a continuous soft pulsing animation around it at all times. Errors are communicated exclusively through the workbench reporting area.
   - *Action Item*: Verify in `map.spec.ts` that the single point marker retains its exact layer color hex code even when outside grid.
4. **Invalid Row Filtering in GeoJSON / DXF Exporters vs. CSV Exact Replica**:
   - *Requirement*: CSV export includes all lines (including comment and error rows), but GeoJSON and DXF exporters must cleanly filter out out-of-grid coordinates, invalid angles, and syntax errors to prevent corrupted vector geometries.
   - *Action Item*: Add tests in `webapp/src/utils/exporters.test.ts` and `batch.spec.ts` asserting DXF and GeoJSON entity counts exclude invalid records.
5. **Spatial Metric Mutual Exclusivity & Global Filter Sync**:
   - *Requirement*: In the Spatial section of Usage, Metric must be a single-choice selector (`Points` OR `Calls`), whereas Computation allows multi-toggle (`Points` AND `Calls`). Spatial KPI cards must sync with all applied filters (time, custom date, sources).
   - *Action Item*: Add test in `legacy.spec.ts` verifying clicking `Calls` unchecks `Points` in the Spatial view and KPI cards reflect active time filters.
6. **URL Deep-Linking Hydration**:
   - *Requirement*: Booting the SPA with query parameters (`?op=Stereo70ToETRS89&c1=500000&c2=500000`) populates inputs and computes output coordinates on load.
   - *Action Item*: Add a new test in `point.spec.ts` asserting that query params successfully hydrate and calculate without user interaction.
7. **Floating Panel Controls (Grid Metadata & Feedback)**:
   - *Requirement*: Dragging headers, resizing windows, search tree filtering, and anti-bot math challenges function reliably without DOM errors.
   - *Action Item*: Add targeted integration tests validating floating panel open, search, minimize, and close cycles.
8. **Brand Logo Fresh-Load Reset Action vs. Title Navigation**:
   - *Requirement*: Clicking the Brand Logo (`#header-brand-logo-wrap`) triggers a complete application state reset (clearing Point inputs/outputs, Batch text, Files queue, map layers, measurements, and closing floating panels) while preserving language and theme. Clicking Brand Title text performs soft navigation to Studio without state destruction.
   - *Action Item*: In `smoke.spec.ts`, assert clicking Brand Logo returns to Studio and purges active coordinate inputs.
9. **Global Offline Alert Banner & Non-Disruptive Health Dot**:
   - *Requirement*: The header health dot remains a fixed 8px status indicator (`.online`/`.offline`) without shifting adjacent navigation elements. A dedicated top alert banner appears when the REST API is offline with an explicit click-to-report link.
   - *Action Item*: In error mock tests, verify that when offline, the global banner renders with an actionable report trigger and opens the Feedback modal.
10. **Untransformable Map Click Coordinate Erase**:
    - *Requirement*: Clicking outside the official TransDatRO grid or on an area without grid data immediately clears both input and output coordinate fields (`""` and `—`) in both ETRS89 and Stereo 70 directions.
    - *Action Item*: In `map.spec.ts`, simulate clicking outside Romania's grid boundary and assert input coordinate boxes are cleared.
11. **Golden Rule Quality Reporting Architecture & Error Banner Elimination**:
    - *Requirement*: Eliminate all ad-hoc inline `.error-banner` elements. Feedback is handled exclusively through the Two-Section Quality Reporting Card (`Coordinates` and `Transformation` separated by a dashed divider) using unified classifications (`Valid`, `Invalid Input`, `Invalid Angle`, `Outside Grid`, `No Grid Data`). DMS angles explicitly show sexagesimal symbols (`°`, `'`, `"`).
    - *Action Item*: In `point.spec.ts` and `batch.spec.ts`, verify that invalid numeric and angular inputs produce structured quality badges without inline error banners.
12. **Batch Textarea Height Persistence & Clean Error Row Formatting**:
    - *Requirement*: Clicking Transform preserves user-adjusted vertical textarea heights. Non-computable input rows output only Point Name (if present) and quality comment (e.g. `104, # Invalid input` or `101, # Invalid angle`), omitting invalid coordinate numbers.
    - *Action Item*: In `batch.spec.ts`, verify custom textarea heights are maintained across calculations and invalid rows omit coordinates in output.
13. **Download Action Toast Elimination**:
    - *Requirement*: Clicking export for CSV, GeoJSON, or DXF triggers native browser save/download dialogs without firing premature synthetic "Downloaded" toast notifications.
    - *Action Item*: In `batch.spec.ts`, assert clicking download buttons does not emit success toasts.
14. **Computation Timeline Total Series Toggle**:
    - *Requirement*: Computation sub-view timeline chart includes an interactive toggle button (`#comp-toggle-total`) to show or hide the combined total series line(s).
    - *Action Item*: In `legacy.spec.ts` / `telemetry.spec.ts`, verify `#comp-toggle-total` toggles the total series line visibility.
