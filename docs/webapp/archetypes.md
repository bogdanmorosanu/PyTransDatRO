# PyTransDatRO Web Application — Component Archetypes & Consistency Catalog

This catalog defines the **14 multi-instance UI element archetypes** across the PyTransDatRO web application. Elements within the same archetype share the same geometry, typography, active/hover transitions, and styling rules.

> [!NOTE]
> Global color tokens and layout invariants are defined in [`design_system.md`](./design_system.md). Operational behaviors and automated test mappings are documented in [`actionable_inventory.md`](./actionable_inventory.md).

---

## 1. Navigation Tabs

### Archetype 1.1: Primary View Tabs
*Top-level navigation between the main application views inside the 52px header.*
- **Styling Standard**: Clean uppercase/sentence-case typography (`font-weight: 500`), transparent surface background, and an active 2px bottom underline indicator in GitHub Primer coral (`--tab-indicator-color`). Inactive tabs dim to `--text-muted` and brighten on hover.
- **Member Instances**:
  - `Transform Tab (#nav-btn-studio)` (Global Header)
  - `Usage Tab (#nav-btn-telemetry)` (Global Header)

### Archetype 1.2: Secondary View Tabs
*Sub-view navigation within a primary workspace.*
- **Styling Standard**: 
  - *Workbench Modes*: Folder bookmark tabs (`.tab-btn`) that physically dock and merge with the `.studio-card` below with zero background gap (`margin-bottom: -1px`, matching surface background and top border).
  - *Usage Sub-Pages*: Underlined sub-navigation tabs (`.telemetry-nav-tab-btn`) with smooth transition indicator.
- **Member Instances**:
  - `Point Mode Tab (#tab-point)` (Workbench Mode Tabs)
  - `Text Mode Tab (#tab-text)` (Workbench Mode Tabs)
  - `Files Mode Tab (#tab-files)` (Workbench Mode Tabs)
  - `Computation Sub-Tab (#tab-timeline-btn)` (Usage Page Sub-Nav)
  - `Spatial Sub-Tab (#tab-spatial-btn)` (Usage Page Sub-Nav)
  - `Legacy Sub-Tab (#tab-legacy-btn)` (Usage Page Sub-Nav)

---

## 2. Selection & Mode Switchers

### Archetype 2.1: Slide Selectors (Animated Glider / Switch)
*Mutually exclusive two-state toggles with a sliding indicator and brief switch animation.*
- **Styling Standard**: Recessed `--bg-surface-input` outer track (`border-radius: var(--radius-sm)`), an elevated neutral `--btn-primary-bg` glider with `1px solid var(--btn-primary-border)`, and a 150ms CSS cubic-bezier slide animation. Inactive options dim to `opacity: 0.65`.
- **Member Instances**:
  - `Transformation Direction Selector (.sliding-toggle)` (Shared Config Bar)  
    *Options: `ETRS89 ➔ Stereo 70` vs `Stereo 70 ➔ ETRS89`*
  - `Coordinate Order Selector (#shared-order-control)` (Shared Config Bar)  
    *Options: `NE(H)` vs `EN(H)`*

### Archetype 2.2: Continuous Segmented Pill Trays (Single-Choice Enclosed Trays)
*Enclosed pill trays containing mutually exclusive options on a shared recessed track.*
- **Styling Standard**: Enclosed container (`.telemetry-btn-toggle-group`) with recessed background, active pill in `--btn-primary-bg` with medium border, compact padding (`4px 8px`), font size `11px`.
- **Member Instances**:
  - `Time Window Presets ([data-preset], [data-spatial-preset], [data-legacy-preset])` (Computation, Spatial, Legacy)  
    *Options: `All Time`, `12M`, `YTD`, `30D`, `7D`*
  - `Granularity Selectors ([data-granularity], [data-legacy-granularity])` (Computation, Legacy)  
    *Options: `Daily`, `Weekly`, `Monthly`*
  - `Spatial Metric Selector (#spatial-metric-points-btn, #spatial-metric-calls-btn)` (Spatial Sub-View)  
    *Options: `Points` vs `Calls` (Single-choice radio tray)*
  - `Spatial Scaling Mode Toggle (#spatial-scale-log-btn, #spatial-scale-linear-btn)` (Spatial Sub-View)  
    *Options: `Logarithmic` vs `Linear`*

### Archetype 2.3: Multi-Select Toggle Buttons (Multi-Check with Status Dot)
*Independent toggle buttons with embedded circular status dot, where multiple options can be active simultaneously.*
- **Styling Standard**: Physically separated buttons (`gap: 8px`, independent borders, no outer enclosing tray), featuring an embedded `6px × 6px` circular status dot (`::before`):
  - **ON**: Solid accent dot (`--color-primary`), `--btn-primary-bg` surface fill, and `--btn-primary-border` border.
  - **OFF**: Hollow ring dot, recessed `--bg-surface-input` background, `opacity: 0.62`.
- **Member Instances**:
  - `Computation Source Filter ([data-source-name])` (Computation Sub-View) — *`Web App`, `REST API`, `Desktop / CLI`*
  - `Computation Metric Filter (#comp-metric-points-btn, #comp-metric-calls-btn)` (Computation Sub-View) — *`Points`, `Calls` (Multi-select)*
  - `Spatial Source Filter ([data-spatial-source-name])` (Spatial Sub-View) — *`Web App`, `REST API`, `Desktop / CLI`*
  - `Legacy Source Filter ([data-legacy-source])` (Legacy Sub-View) — *`Web App`, `REST API (Service)`*
  - `Legacy Operation Filter ([data-legacy-op])` (Legacy Sub-View) — *`Stereo 70`, `Stereo 30`*

---

## 3. Menus & Header Actions

### Archetype 3.1: Header Flyout / Disclosure Menus
*Clickable disclosure buttons in the header displaying a chevron and opening a dropdown card.*
- **Styling Standard**: HTML5 `<details><summary>` disclosure wrapper (`.nav-menu-trigger`), rotating chevron (`▾`), elevated dropdown card with subtle border, and keyboard navigation (`Escape` closes menu).
- **Member Instances**:
  - `Help Menu (#header-help-details)` (Global Header) — *Icon + "Help" + chevron `▾`*
  - `Language Selector Menu (#header-lang-details)` (Global Header) — *Globe icon + active language ("RO" / "EN") + chevron `▾` with checkmarked radio options*

### Archetype 3.2: Header Icon Action Buttons
*Single-click circular icon buttons without text labels.*
- **Styling Standard**: Compact circular button (`.btn-icon`, 32×32px, `border-radius: 50%`), transparent or subtle surface fill, brightening on hover.
- **Member Instances**:
  - `Theme Toggle Button (#header-theme-btn)` (Global Header) — *Swaps Sun / Moon icon to toggle Dark / Light mode*
  - `Mobile View Switcher Button (#header-mobile-view-btn)` (Global Header) — *Visible on mobile (<768px) to toggle Workbench vs Map*

---

## 4. Typography: Labels & Section Headers

### Archetype 4.1: Primary Labels
*Prominent field and section identifiers defining primary controls.*
- **Styling Standard**: Clean, medium-weight typography (`font-size: 12px–12.5px`, `font-weight: 500`, `color: var(--text-muted)`).
- **Member Instances**:
  - `Angular Unit Label (.config-field-label)` (Shared Config Bar) — *"Angular Unit"*
  - `Coordinate Order Label (.config-field-label)` (Shared Config Bar) — *"Coordinate Order"*
  - `Angle Format Unit Label (.format-example-label)` (Shared Config Bar) — *"Format Example:"*

### Archetype 4.2: Secondary Labels
*Subordinate labels for secondary metrics or preview readouts.*
- **Styling Standard**: Sentence-case, neutral label (`font-size: 11.5px`, `color: var(--text-muted)`).
- **Member Instances**:
  - `Quasigeoid Height Label (.form-label)` (Point Output Card) — *"Quasigeoid Height (ζ = h - H)"*
  - `Angle Format Example Label / Prefix (.format-example-label)` (Shared Config Bar)

### Archetype 4.3: Element & Component Labels
*Standard field labels directly identifying inputs, outputs, filter groups, and date bounds.*
- **Styling Standard**: Consistent uppercase or sentence-case label (`font-size: 11px–12px`, `color: var(--text-muted)`).
- **Member Instances**:
  - `Coordinate 1 Input Label (.form-label)` (Point Mode) — *"Northing (X)"* / *"Latitude (B)"*
  - `Coordinate 2 Input Label (.form-label)` (Point Mode) — *"Easting (Y)"* / *"Longitude (L)"*
  - `Altitude/Height Input Label (.form-label)` (Point Mode) — *"Altitude (H)"* / *"Ellipsoidal Height (h)"*
  - `Coordinate 1 Output Label (.form-label)` (Point Mode) — *"Latitude (B)"* / *"Northing (X)"*
  - `Coordinate 2 Output Label (.form-label)` (Point Mode) — *"Longitude (L)"* / *"Easting (Y)"*
  - `Altitude/Height Output Label (.form-label)` (Point Mode) — *"Ellipsoidal Height (h)"* / *"Altitude (H)"*
  - `Batch Input Label (.form-label)` (Text Mode) — *"Input coordinates (multi-line)"*
  - `Batch Output Label (.form-label)` (Text Mode) — *"Transformed coordinates"*
  - `Telemetry Facet Group Labels (.facet-label)` (Usage Dashboard) — *"SOURCES"*, *"METRIC"*, *"OPERATIONS"*
  - `Earliest Available Date Label (.limit-subtle)` (Usage Date Limits) — *"Earliest Available:"*
  - `Current / Today Date Label (.limit-subtle)` (Usage Date Limits) — *"Current Date:"* / *"Today:"*

### Archetype 4.4: Discrete / Subtle Section Labels
*Faded, low-prominence labels that divide sub-sections without drawing visual attention.*
- **Styling Standard**: Uppercase, small, muted typography (`font-size: 10px–10.5px`, `letter-spacing: 0.5px`, `color: var(--text-dim)`).
- **Member Instances**:
  - `Coordinates Quality Sub-Section Header (.quality-section-header)` (Point / Batch / Files Reporting) — *"Coordinates"* / *"Coordonate"*
  - `Transformation Quality Sub-Section Header (.quality-section-header)` (Point / Batch / Files Reporting) — *"Transformation"* / *"Transformare"*
  - `Download As: Export Header (.download-group-label)` (Text Mode) — *"Download As:"*

---

## 5. Input & Output Fields

### Archetype 5.1: Monospace Text Input Fields
*User-editable single-line input boxes accepting coordinates.*
- **Styling Standard**: Monospace font (`JetBrains Mono`), recessed `--bg-surface-input`, `--border-medium` border, `--border-focus` focus ring, auto-selection on focus (`input.select()`).
- **Member Instances**:
  - `Coordinate Input 1 (#pt-input-1)` (Point Mode)
  - `Coordinate Input 2 (#pt-input-2)` (Point Mode)
  - `Altitude / Height Input 3 (#pt-input-3)` (Point Mode)

### Archetype 5.2: Monospace Read-Only Output Fields
*Non-editable output boxes displaying computed coordinate values with a persistent em-dash `—` when empty.*
- **Styling Standard**: Monospace typography, non-editable but selectable, elevated surface background (`--bg-surface-elevated`), persistent em-dash `—` when uncalculated.
- **Member Instances**:
  - `Coordinate Output 1 (.output-grid-3 input:nth-child(1))` (Point Mode)
  - `Coordinate Output 2 (.output-grid-3 input:nth-child(2))` (Point Mode)
  - `Altitude / Height Output 3 (.output-grid-3 input:nth-child(3))` (Point Mode)

### Archetype 5.3: Monospace Multi-Line Text Areas
*Vertically resizable multi-line text editors for batch processing.*
- **Styling Standard**: Monospace font, line-height 1.4, vertical resize grip, height persistence across calculations, default arrow cursor over scrollbar gutter.
- **Member Instances**:
  - `Batch Input Textarea (#batch-input-area)` (Text Mode) — *Editable with direction placeholder*
  - `Batch Output Textarea (#batch-output-area)` (Text Mode) — *Read-only with elevated background*

### Archetype 5.4: Selectable Informational Readouts
*Informational text readouts that the user can highlight, select, and copy.*
- **Styling Standard**: Monospace font, high-contrast readable text, selectable cursor, optional copy action.
- **Member Instances**:
  - `Angle Format Example Value (#angle-format-example-badge)` (Shared Config Bar)  
    *Selectable / copyable coordinate preview (e.g., `45°50'28.44420"`).*
  - `Quasigeoid Height Value (.metric-field)` (Point Output Card)  
    *Selectable / copyable numerical scalar value in cyan metric palette (e.g., `+0.421`).*
  - `Earliest Available Date Value (.limit-date-val)` (Usage Dashboard)  
    *Selectable historical boundary date (e.g., `2023-11-01`).*
  - `Current / Today Date Value (.limit-date-val)` (Usage Dashboard)  
    *Selectable current boundary date (e.g., `2026-10-22`).*

### Archetype 5.5: Native Themed Date Pickers
*Themed calendar date inputs with consistent borders, calendar icons, and boundary clamping.*
- **Styling Standard**: Native `<input type="date">` styled with `color-scheme`-aware theme colors, `--border-medium` outline, and subtle border radius (`--radius-sm`).
- **Member Instances**:
  - `Computation Date Range Inputs (#comp-date-start, #comp-date-end)` (Computation Sub-View)
  - `Spatial Date Range Inputs (#spatial-date-start, #spatial-date-end)` (Spatial Sub-View)
  - `Legacy Date Range Inputs (#legacy-date-start, #legacy-date-end)` (Legacy Sub-View)

---

## 6. Action Buttons

### Archetype 6.1: Primary Action CTA Buttons
*The prominent full-width action button that executes the transformation calculation.*
- **Styling Standard**: Option #8 Subtle Surface + Border (`--btn-primary-bg`, `--btn-primary-border`, `--btn-primary-text`), full panel width, loading label/spinner, zero decorative emojis.
- **Member Instances**:
  - `Point Transform Button (#pt-calculate-btn)` (Point Mode) — *"Transform" / "Calculating..."*
  - `Batch Transform Button (#batch-transform-btn)` (Text Mode) — *"Transform" / "Calculating..."*
  - `Files Transform Button (#file-start-btn)` (Files Mode) — *"Transform Files" / "Streaming progress..."*

### Archetype 6.2: Code-Block Copy Action Buttons
*Discrete top-right icon buttons with instant checkmark feedback to copy computed values to clipboard.*
- **Styling Standard**: Small icon button (`.copy-code-btn`, `.btn-icon`), subtle border, instant checkmark feedback (`✓ Copied`), zero layout shift.
- **Member Instances**:
  - `Copy Point Coordinates Button (#pt-copy-code-btn)` (Point Mode)
  - `Copy Batch Results Button (#batch-copy-results-slot)` (Text Mode)

### Archetype 6.3: Export & File Action Buttons
*Secondary action buttons triggering file downloads or exports.*
- **Styling Standard**: Subtle secondary surface fill (`.btn-secondary`), 12px typography, padding `6px 12px`, direct invocation of native browser Save As dialog without premature toasts.
- **Member Instances**:
  - `Download CSV Button (#batch-export-csv)` (Text Mode)
  - `Download GeoJSON Button (#batch-export-geojson)` (Text Mode)
  - `Download DXF Button (#batch-export-dxf)` (Text Mode)
  - `Individual File Download Button (.file-download-btn)` (Files Mode)
  - `Download All Completed Button (#file-download-all-btn)` (Files Mode)

### Archetype 6.4: Queue & List Item Action Buttons
*Utility buttons managing items in a list or staging queue.*
- **Styling Standard**: Subtle outline button (`.btn-outline`) or compact `✕` remove icon (`.btn-icon`).
- **Member Instances**:
  - `Clear Queue Button (#file-clear-queue-btn)` (Files Mode) — *Outlined button "✕ Clear Queue"*
  - `Remove File Item Button (.file-queue-remove)` (Files Mode) — *Compact "✕" icon button per row*

### Archetype 6.5: Map Navigation Action Buttons
*Secondary utility buttons triggering Leaflet map panning and boundary fitting.*
- **Styling Standard**: Outlined button (`.btn-outline`), font size `11.5px`, padding `4px 10px`, accompanied by crosshair or search icon, guarded by Romania bounds.
- **Member Instances**:
  - `Center on Map (#pt-center-map-btn, #batch-center-map-btn, #file-center-map-btn)` (Point, Text, Files)
  - `Zoom to Extents (#pt-zoom-extents-btn, #batch-zoom-extents-btn, #file-zoom-extents-btn)` (Point, Text, Files)

### Archetype 6.6: Graph Series Total Toggles
*Standalone toggle buttons turning the cumulative Total series curve on/off on SVG time charts.*
- **Styling Standard**: Independent toggle button with status indicator, mirroring source toggles.
- **Member Instances**:
  - `Computation Total Overlay Toggle (#comp-toggle-total)` (Computation Sub-View)
  - `Legacy Total Overlay Toggle (#legacy-toggle-total)` (Legacy Sub-View)

---

## 7. Quality Reporting & Feedback Labels

### Archetype 7.1: Quality Info Labels
*Neutral, informational pills indicating point counts and dimensional traits.*
- **Styling Standard**: Elevated background (`--bg-surface-elevated`), neutral border (`--border-medium`), text color `--text-muted`, corner radius `--radius-xs` (4px).
- **Member Instances**:
  - `Total Points Count (.quality-report-badges .status-pill)` (Batch / Files)
  - `2D Coordinates Pill / Count (.quality-report-badges .status-pill)` (Point / Batch / Files)
  - `3D Coordinates Pill / Count (.quality-report-badges .status-pill)` (Point / Batch / Files)

### Archetype 7.2: Quality Validation Labels
*Outcome and diagnostic pills rendered with semantic status colors (Green, Amber, Red).*
- **Styling Standard**: High-contrast semantic pill (`border-radius: var(--radius-xs)`):
  - Green (`#10b981`): `Valid`
  - Amber (`#f59e0b`): `Outside Grid`, `No Grid Data`
  - Red (`#ef4444`): `Invalid Input`, `Invalid Angle`
- **Member Instances**:
  - `Valid Transformation (.status-pill.success)` (Point / Batch / Files) — *Green*
  - `Outside Grid (.status-pill.warning)` (Point / Batch / Files) — *Amber*
  - `No Grid Data (.status-pill.warning)` (Point / Batch / Files) — *Amber*
  - `Invalid Input (.status-pill.danger)` (Point / Batch / Files) — *Red*
  - `Invalid Angle (.status-pill.danger)` (Point / Batch / Files) — *Red*

### Archetype 7.3: File Processing Status Badges
*Status indicators for individual files queued in File Upload mode.*
- **Styling Standard**: Small rounded pill indicating processing lifecycle state.
- **Member Instances**:
  - `Pending Badge (.file-status-badge.pending)` (Files Mode)
  - `Processing Badge (.file-status-badge.processing)` (Files Mode)
  - `Completed Badge (.file-status-badge.completed)` (Files Mode)
  - `Error Badge (.file-status-badge.error)` (Files Mode)

---

## 8. Analytics & Telemetry KPI Cards

### Archetype 8.1: Main KPI Cards (Single Metric)
*Standard standalone KPI cards with one large numerical value and one descriptive label below it.*
- **Styling Standard**: Elevated card surface (`--bg-surface`), subtle 1px border, large monospace metric value (`font-size: 20px–22px`, `--text-main`), quiet uppercase label below (`--text-muted`, `font-size: 11px`).
- **Member Instances**:
  - `Total Points Card (.telemetry-kpi-card)` (Computation Sub-View)
  - `Total Calls Card (.telemetry-kpi-card)` (Computation Sub-View)
  - `Points / Call Card [Batch Intensity] (.telemetry-kpi-card)` (Computation Sub-View)
  - `Legacy Total Points Card (.telemetry-kpi-card)` (Legacy Sub-View)
  - `Legacy Web App Points Card (.telemetry-kpi-card)` (Legacy Sub-View)
  - `Legacy Service API Points Card (.telemetry-kpi-card)` (Legacy Sub-View)
  - `Legacy Stereo 70 Points Card (.telemetry-kpi-card)` (Legacy Sub-View)
  - `Legacy Stereo 30 Points Card (.telemetry-kpi-card)` (Legacy Sub-View)

### Archetype 8.2: Main KPI Cards (Merged Multi-Column Banner)
*A unified full-width card container divided into distinct metric columns.*
- **Styling Standard**: Full-width elevated banner with vertical separator line between metric columns.
- **Member Instances**:
  - `Spatial Overview Banner (.spatial-overview-banner)` (Spatial Sub-View)  
    *Columns: Total Points/Calls (Left) + Active Cells Count (Right)*

### Archetype 8.3: Secondary Regional KPI Cards (Quadrants)
*Compact geographic distribution cards displaying quadrant name above, numerical value below, and percentage share alongside.*
- **Styling Standard**: 4-column deck (`.spatial-quadrants-deck`), quadrant code above (`NW`, `NE`, `SW`, `SE`), monospace value + percentage share (`font-size: 12px`, `--text-muted`).
- **Member Instances**:
  - `Quadrant NW Card (.quadrant-card.nw)` (Spatial Sub-View)
  - `Quadrant NE Card (.quadrant-card.ne)` (Spatial Sub-View)
  - `Quadrant SW Card (.quadrant-card.sw)` (Spatial Sub-View)
  - `Quadrant SE Card (.quadrant-card.se)` (Spatial Sub-View)
