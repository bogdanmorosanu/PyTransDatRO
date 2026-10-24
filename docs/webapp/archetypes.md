# PyTransDatRO Web Application — Component Archetypes & Consistency Catalog

This catalog defines the **27 multi-instance UI element archetypes across 9 functional categories** in the PyTransDatRO web application. Elements within the same archetype share the same geometry, typography, active/hover transitions, and styling rules.

> [!NOTE]
> Global color tokens and layout invariants are defined in [`design_system.md`](./design_system.md). Operational behaviors and automated test mappings are documented in [`actionable_inventory.md`](./actionable_inventory.md).

---

## 1. Navigation Tabs

### Archetype 1.1: Primary View Tabs
*Top-level navigation between the main application views inside the 52px header.*
- **Styling Standard**: Prominent sentence-case typography (`font-size: 13.5px`, `font-weight: 500` inactive, `font-weight: 600` active), height `32px` with `border-radius: var(--radius-sm)` (6px). On hover, renders a discrete neutral gray pill background (`--tab-hover-bg`) with a clean 10px vertical gap above the bottom header border. On active, displays an active 2.5px bottom underline indicator in neutral gray (`--tab-indicator-color`: `#d0d7de` Dark / `#1f2328` Light) flush with the header boundary.
- **Member Instances**:
  - `Transform Tab (#nav-btn-studio)` (Global Header)
  - `Usage Tab (#nav-btn-telemetry)` (Global Header)

### Archetype 1.2: Secondary View Tabs
*Sub-view navigation within primary workspaces (Studio Workbench & Usage Dashboard).*
- **Styling Standard**: Unified underline navigation tabs across all secondary workspaces. Container height `40px`, button height `28px`, `font-size: 12.5px`, `padding: 0 10px`, `border-radius: var(--radius-sm)`. Inactive tabs show neutral hover pill (`--tab-hover-bg`), active tabs display bold text and an active 2px bottom underline indicator in `--tab-indicator-color` (`#d0d7de` Dark / `#1f2328` Light).
- **Member Instances**:
  - `Point Mode Tab (#tab-point)` (Studio Workbench Tabs)
  - `Text Mode Tab (#tab-text)` (Studio Workbench Tabs)
  - `Files Mode Tab (#tab-files)` (Studio Workbench Tabs)
  - `Computation Sub-Tab (#tab-timeline-btn)` (Usage Page Sub-Nav)
  - `Spatial Sub-Tab (#tab-spatial-btn)` (Usage Page Sub-Nav)
  - `Legacy Sub-Tab (#tab-legacy-btn)` (Usage Page Sub-Nav)

---

## 2. Selection & Mode Switchers

### Archetype 2.1: Slide Selectors (Animated Glider / Switch)
*Mutually exclusive two-state toggles with an animated sliding indicator and smooth 240ms cubic-bezier transition.*
- **Styling Standard**: Unified compact 28px height (`height: 28px; padding: 2px; border-radius: var(--radius-sm)`), recessed `--bg-surface-input` outer track, elevated neutral `--btn-primary-bg` glider with `1px solid var(--btn-primary-border)`, and a 240ms CSS cubic-bezier slide animation (`transform: translateX(...)`). Inactive options dim to `--text-dim` (`font-size: 11.5px; font-weight: 500`).
- **Member Instances**:
  - `Transformation Direction Selector (.sliding-toggle)` (Shared Config Card, Row 1)  
    *Options: `ETRS89 to Stereo70` vs `Stereo70 to ETRS89`*
  - `Coordinate Order Selector (#shared-order-control)` (Shared Config Card, Row 2)  
    *Options: `NE(H) / LatLon(h)` vs `EN(H) / LonLat(h)`*

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
  - `Computation Total Overlay Toggle (#comp-toggle-total)` (Computation Sub-View) — *Toggle button with dashed total line sample*
  - `Spatial Source Filter ([data-spatial-source-name])` (Spatial Sub-View) — *`Web App`, `REST API`, `Desktop / CLI`*
  - `Legacy Source Filter ([data-legacy-source])` (Legacy Sub-View) — *`Web App`, `REST API (Service)`*
  - `Legacy Operation Filter ([data-legacy-op])` (Legacy Sub-View) — *`Stereo 70`, `Stereo 30`*
  - `Legacy Total Overlay Toggle (#legacy-toggle-total)` (Legacy Sub-View) — *Toggle button with dashed total line sample*

### Archetype 2.4: Boolean Toggle Switch (Sliding Track & Thumb)
*Compact, accessible on/off toggle switch for binary options within input cards.*
- **Styling Standard**: Semantic checkbox hidden with `.toggle-switch-input`, styled via `.toggle-switch-track` (30×17px, border radius 9px, `--bg-surface-elevated`) and sliding `.toggle-switch-thumb` (12×12px circular nub, transitions 12px rightward to `--color-primary` when `:checked`), accompanied by `.toggle-switch-label` (11.5px, medium weight).
- **Member Instances**:
  - `Ignore Header / First Line Switch (#file-skip-header-toggle)` (Files Input Card) — *Toggles ignoring the header row in uploaded bulk files.*


---

## 3. Menus & Header Actions

### Archetype 3.1: Header Flyout / Disclosure Menus
*Clickable disclosure buttons in the header displaying a chevron and opening a dropdown card.*
- **Styling Standard**: HTML5 `<details><summary>` disclosure wrapper (`.nav-menu-trigger`), rotating chevron (`▾`), elevated dropdown card with subtle border, and keyboard navigation (`Escape` closes menu). All action and navigation items within the Help menu dropdown are purely typographical (icon-free) to preserve visual calm and consistency with external links.
- **Member Instances**:
  - `Help Menu (#header-help-details)` (Global Header) — *"Help"* + chevron `▾` (purely typographical item ecosystem)
  - `Language Selector Menu (#header-lang-details)` (Global Header) — *Globe icon + active language ("RO" / "EN") + chevron `▾` with checkmarked radio options*

### Archetype 3.2: Header Icon Action Buttons
*Single-click circular icon buttons without text labels.*
- **Styling Standard**: Compact circular button (`.btn-icon`, 32×32px, `border-radius: 50%`), transparent or subtle surface fill, brightening on hover.
- **Member Instances**:
  - `Theme Toggle Button (#header-theme-btn)` (Global Header) — *Swaps Sun / Moon icon to toggle Dark / Light mode*
  - `Mobile View Switcher Button (#header-mobile-view-btn)` (Global Header) — *Visible on mobile (<768px) to toggle Workbench vs Map*

---

## 4. Typography: Labels & Section Headers

### Archetype 4.1: Standard Field & Control Labels
*Standard field labels directly identifying inputs, outputs, selects, and filter group decks.*
- **Styling Standard**: Consistent sentence-case or title-case label (`font-size: 11.5px–12.5px`, `font-weight: 500`, `color: var(--text-muted)`, `margin-bottom: 5px`, `text-transform: none`). Trailing colons (`:`) are strictly omitted.
- **Member Instances**:
  - `Operation Direction Label (.form-label)` (Shared Config Card) — *"Operation Direction"* / *"Sensul transformării"*
  - `Angular Unit Label (.form-label)` (Shared Config Card) — *"Angular Unit"* / *"Unitate unghiulară"*
  - `Coordinate Order Label (.form-label)` (Shared Config Card) — *"Coordinate Order"* / *"Ordinea coordonatelor"*
  - `Coordinate 1 Input Label (.form-label)` (Point Mode) — *"Northing (X)"* / *"Latitude (B)"*
  - `Coordinate 2 Input Label (.form-label)` (Point Mode) — *"Easting (Y)"* / *"Longitude (L)"*
  - `Altitude/Height Input Label (.form-label)` (Point Mode) — *"Altitude (H)"* / *"Ellipsoidal Height (h)"*
  - `Coordinate 1 Output Label (.form-label)` (Point Mode) — *"Latitude (B)"* / *"Northing (X)"*
  - `Coordinate 2 Output Label (.form-label)` (Point Mode) — *"Longitude (L)"* / *"Easting (Y)"*
  - `Altitude/Height Output Label (.form-label)` (Point Mode) — *"Ellipsoidal Height (h)"* / *"Altitude (H)"*
  - `Batch Input Label (.form-label)` (Text Mode) — *"Input coordinates (multi-line)"*
  - `Batch Output Label (.form-label)` (Text Mode) — *"Transformed coordinates"*
  - `Telemetry Facet Group Labels (.form-label)` (Usage Dashboard) — *"Sources"*, *"Metrics"* (Computation) / *"Metric"* (Spatial), *"Operations"*
  - `Active Period Main Label (.form-label.telemetry-active-main-label)` (Usage Active Period) — *"Active Period"*
  - `Report Type Label (.form-label)` (Feedback Drawer) — *"Report Type"* / *"Tipul semnalării"*
  - `Grid Metadata Labels (.fgp-meta-label)` (Grid Metadata Drawer) — *"Size"*, *"SHA-256"*

### Archetype 4.2: Sub-Labels & Contextual Bounds
*Subordinate, discrete labels for preview readouts, boundary dates, and metric qualifiers.*
- **Styling Standard**: Sentence-case, neutral label (`font-size: 11.5px`, `font-weight: 500`, `color: var(--text-dim)` or `--text-muted`).
- **Member Instances**:
  - `Format Example Label (.info-label)` (Shared Config Card) — *"Format Example"* / *"Exemplu format"*
  - `Quasigeoid Height Label (.info-label)` (Point Output Card) — *"Quasigeoid Height (ζ = h - H)"*
  - `Earliest Available Date Label (.info-label)` (Usage Date Limits) — *"Earliest Available"*
  - `Current / Today Date Label (.info-label)` (Usage Date Limits) — *"Current Date"*
  - `Active Period Start & End Sub-Labels (.info-label)` (Usage Active Period) — *"Start"*, *"End"*

### Archetype 4.3: Discrete Section Landmarks
*Faded, low-prominence uppercase landmark labels that organize functional cards without drawing visual attention.*
- **Styling Standard**: Uppercase, small, muted typography (`font-size: 10.5px`, `font-weight: 600`, `letter-spacing: 0.6px`, `color: var(--text-dim)`, `margin-bottom: 8px`, `line-height: 1.2`).
- **Member Instances**:
  - `Settings Landmark (.section-micro-label)` (Global Config Card) — *"SETTINGS"* / *"SETĂRI"*
  - `Input Landmark (.section-micro-label)` (Point & Text Modes) — *"INPUT"* / *"INTRARE"*
  - `Output Landmark (.section-micro-label)` (Point & Text Modes) — *"OUTPUT"* / *"IEȘIRE"*
  - `Report Landmark (.section-micro-label)` (Point & Text Modes) — *"REPORT"* / *"RAPORT"*
  - `Coordinates Quality Header (.quality-label)` (Reporting Card) — *"Coordinates"* / *"Coordonate"* (Sentence case, discrete)
  - `Transformation Quality Header (.quality-label)` (Reporting Card) — *"Transformation"* / *"Transformare"* (Sentence case, discrete)
  - `Download as: Export Header (.download-group-label)` (Text Mode) — *"Download as:"* / *"Descarcă ca:"*

### Archetype 4.4: Content Topic Headings & Explanatory Blocks (Tier 2 Subheadings)
*Structured topic headings paired with descriptive body text inside documentation and information cards.*
- **Styling Standard**: Two-tier explanatory hierarchy. Subheading (`.content-card-title`, `.topic-heading`, `.fhp-topic-heading`): `font-size: 11.5px`, `font-weight: 600`, `color: var(--text-main)`, subordinate to the 13.5px Section Header. Successive topic items are divided by a subtle dashed line (`border-top: 1px dashed var(--border-dashed); padding-top: 8px`). Body (`.topic-desc`, `.fhp-topic-desc`): `font-size: 11.5px`, `line-height: 1.5`, `color: var(--text-muted)`. **The Golden Rule**: Trailing colons (`:`) are strictly omitted from all headings.
- **Member Instances**:
  - `Stereo 70 Topic Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"Stereo 70 (Dealul Piscului 1970 · EPSG:31700)"*
  - `ETRS89 Topic Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"ETRS89 (GRS80 · EPSG:4258)"*
  - `Altimetry Topic Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"Altimetry & Quasigeoid (Black Sea 1975)"*
  - `Angular Formats Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"Angular Formats (ETRS89)"*
  - `Coordinate Order Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"Coordinate Order"*
  - `Dimensionality Auto-Detection Heading (.fhp-topic-heading)` (Geodetic Help Panel) — *"Dimensionality Auto-Detection (2D / 3D)"*

### Archetype 4.5: Informational Contextual Tooltip Indicators
*Discrete, non-actionable hoverable icons providing contextual domain guidance or explanations without button semantics.*
- **Styling Standard**: GitHub-style hoverable SVG indicator (`.info-help-indicator`). Rendered as an accessible inline element (`<span class="info-help-indicator" tabindex="0">`) with `cursor: help`, transparent background, zero border, and muted color (`--text-dim` or `--text-muted`) that smoothly transitions to `--text-main` on hover/focus. Pairs with an embedded theme-reactive contextual tooltip bubble (`.info-tooltip-bubble`, `role="tooltip"`) that strictly honors Dark/Light mode (`background-color: var(--bg-surface)`, `color: var(--text-main)`, `border: 1px solid var(--border-medium)`, `box-shadow: var(--shadow-md)`), avoiding native OS black box tooltip popovers on light backgrounds.
- **Member Instances**:
  - `Active Period Explanation Indicator (.info-help-indicator)` (Usage Dashboard — Computation & Spatial) — *Contextual help icon explaining monthly and quarterly telemetry database storage snapping.*

---

## 5. Input & Output Fields

### Archetype 5.1: Monospace Text Input Fields
*User-editable single-line input boxes accepting coordinates.*
- **Styling Standard**: Monospace font (`JetBrains Mono`), compact 30px height (`height: 30px; padding: 4px 10px; font-size: 12.5px`), recessed `--bg-surface-input`, `--border-medium` border, `--border-focus` focus ring, auto-selection on focus (`input.select()`).
- **Member Instances**:
  - `Coordinate Input 1 (#pt-input-1)` (Point Mode)
  - `Coordinate Input 2 (#pt-input-2)` (Point Mode)
  - `Altitude / Height Input 3 (#pt-input-3)` (Point Mode)

### Archetype 5.2: Monospace Read-Only Output Fields
*Non-editable output boxes displaying computed coordinate values, remaining strictly empty (`""`) when uncalculated or untransformable.*
- **Styling Standard**: Monospace typography, compact 30px height (`height: 30px; padding: 4px 10px; font-size: 12.5px`), non-editable but selectable, elevated surface background (`--bg-surface-elevated`), remaining empty (`""`) when uncalculated or untransformable.
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
    *Selectable historical boundary date in UK format (e.g., `01/11/23`).*
  - `Current / Today Date Value (.limit-date-val)` (Usage Dashboard)  
    *Selectable current boundary date in UK format (e.g., `24/10/26`).*
  - `Active Period Start & End Values (.active-period-value)` (Usage Dashboard)  
    *Selectable effective snapped period boundaries in UK format (e.g., `Start [ 01/01/24 ]` and `End [ 31/12/24 ]`).*
  - `Benchmark Residual Metric Chips (.fhp-metric-chip)` (Geodetic Help Panel)  
    *Monospace read-only info chips displaying numerical planimetric arcseconds / mm and elevation mm without inequality clutter.*
  - `Angle Units Format Preview Chips (.fhp-format-chip)` (Geodetic Help Panel)  
    *Selectable coordinate format chips (e.g. `46°00'00.00000"`, `46.00000000`, `0.8028514559`).*
  - `Grid File Info Chips (.fgp-meta-chips-row .fhp-metric-chip)` (Grid Metadata Panel)  
    *Selectable / copyable file size and SHA-256 hash chips without colons.*

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

### Archetype 6.2: Secondary & Utility Action Buttons
*Secondary buttons triggering auxiliary actions, exports, list operations, map navigation, and modal confirmations.*
- **Styling Standard**: Compact secondary action buttons (`.btn.btn-secondary`, `.btn.btn-outline`, `.btn-secondary-action`, `.btn-sm`), height 26px–30px, typography 11.5px–12px, padding `0 9px` to `0 12px`, border-radius `--radius-sm` (6px).
- **Member Instances**:
  - `Batch Export Buttons (#batch-export-csv, #batch-export-geojson, #batch-export-dxf)` (Text Mode)
  - `File Download Buttons (.file-download-btn, #file-download-all-btn)` (Files Mode)
  - `Clear Queue & Output Buttons (#file-clear-queue-btn, #file-clear-output-btn)` (Files Mode)
  - `Map Navigation Dock Buttons (#pt-center-map-btn, #batch-center-map-btn, #file-center-map-btn, #pt-zoom-extents-btn, etc.)` (Point, Text, Files)
  - `Modal Footer Buttons (#modal-close-btn, #fb-submit-btn, #fb-cancel-btn)` (Modals & Drawers)
  - `Feedback Auxiliary Buttons (#fb-reset-template-btn, #fb-copy-btn, #fb-github-btn)` (Feedback Drawer)

### Archetype 6.3: Icon Action & Copy Buttons
*Compact icon-only triggers for copy-to-clipboard actions, drawer window controls, and list item removal.*
- **Styling Standard**: Discrete icon button (`.btn-icon`, `.copy-code-btn`, `.drawer-icon-btn`), height 24px–32px, transparent background brightening on hover, instant checkmark feedback (`✓ Copied`) for copy actions without layout shift.
- **Member Instances**:
  - `Copy Point Coordinates Button (#pt-copy-code-btn)` (Point Mode)
  - `Copy Batch Results Button (#batch-copy-results-slot)` (Text Mode)
  - `Drawer Minimize & Close Buttons (.drawer-icon-btn, #fhp-close-btn, #ffp-close-btn, #fgp-close-btn)` (Drawers)
  - `Remove File Item Button (.file-queue-remove)` (Files Mode)

---

## 7. Badges & Diagnostic Feedback

### Archetype 7.1: Diagnostic Status & Validation Pills
*Compact status pills indicating calculation validity, error notices, or neutral point counts.*
- **Styling Standard**: High-contrast semantic pill (`.status-pill`, `border-radius: var(--radius-sm)` or `var(--radius-xs)`):
  - Neutral (`.status-pill.neutral`, `.status-pill-neutral`): `--bg-surface-elevated`, `--text-muted`, `--border-subtle` (`Ready to transform`, `Total (15)`, `2D (10)`, `3D (5)`).
  - Success (`.status-pill.success`): Green (`#10b981`), `--color-success-bg` (`Valid`).
  - Warning (`.status-pill.warning`): Amber (`#f59e0b`), `--color-warning-bg` (`Outside Grid`, `No Grid Data`).
  - Danger (`.status-pill.danger`): Red (`#ef4444`), `--color-danger-bg` (`Invalid Input`, `Invalid Angle`).
- **Member Instances**:
  - `Ready to Transform Pill (.status-pill.neutral)` (Point / Batch / Files)
  - `Dimension Breakdown Pills (.quality-report-badges .status-pill)` (Point / Batch / Files)
  - `Valid Transformation Pill (.status-pill.success)` (Point / Batch / Files)
  - `Outside Grid / No Data Warning Pills (.status-pill.warning)` (Point / Batch / Files)
  - `Invalid Syntax / Angle Danger Pills (.status-pill.danger)` (Point / Batch / Files)

### Archetype 7.2: Queue Lifecycle Badges
*Status indicators representing asynchronous file processing states.*
- **Styling Standard**: Small rounded pill (`.file-status-badge`) indicating queue lifecycle: pending (neutral), processing (primary glow), completed (success green), error (danger red).
- **Member Instances**:
  - `Pending Badge (.file-status-badge.pending)` (Files Mode)
  - `Processing Badge (.file-status-badge.processing)` (Files Mode)
  - `Completed Badge (.file-status-badge.completed)` (Files Mode)
  - `Error Badge (.file-status-badge.error)` (Files Mode)

---

## 8. Analytics & Telemetry KPI Displays

### Archetype 8.1: KPI Metric Cards & Banners
*Elevated numerical readout containers displaying key volume indicators with large monospace figures and subtle descriptions.*
- **Styling Standard**: Elevated surface (`.telemetry-kpi-card`, `.spatial-overview-banner`, `--bg-surface`), 1px border (`--border-subtle`), prominent monospace metric value (`20px–22px`, `--text-main`), quiet sentence-case label (`11px`, `--text-muted`).
- **Member Instances**:
  - `Total Points / Calls / Intensity Cards (.telemetry-kpi-card)` (Computation & Legacy)
  - `Spatial Overview Banner (.spatial-overview-banner)` (Spatial Sub-View)

### Archetype 8.2: Regional Geographic Quadrant Cards
*Compact geographic distribution cards displaying quadrant code above, numerical volume below, and percentage share alongside.*
- **Styling Standard**: 4-column deck (`.spatial-quadrants-deck`), quadrant code above (`NW`, `NE`, `SW`, `SE`), monospace value + percentage share (`font-size: 12px`, `--text-muted`).
- **Member Instances**:
  - `Quadrant NW Card (.quadrant-card.nw)` (Spatial Sub-View)
  - `Quadrant NE Card (.quadrant-card.ne)` (Spatial Sub-View)
  - `Quadrant SW Card (.quadrant-card.sw)` (Spatial Sub-View)
  - `Quadrant SE Card (.quadrant-card.se)` (Spatial Sub-View)

---

## 9. Content Sections & Notice Blocks

### Archetype 9.1: Content Section Cards
*Universally shared elevated modular cards grouping related informational topics, benchmarks, documentation, and user reporting.*
- **Styling Standard**: Shared card container (`.content-card`, `.fhp-card`, `.ffp-card`), elevated background (`--bg-surface-elevated`), 1px subtle border (`--border-subtle`), rounded corners (`--radius-md`), standard internal padding (`14px 16px`), flex-column layout with gap `10px`. Features an authoritative Tier 1 title (`.content-card-title`, `13.5px`, `600` weight, `--text-main`) with a solid divider underline (`border-bottom: 1px solid var(--border-subtle); padding-bottom: 8px`), optional subtitle (`11.5px`, `--text-muted`), and structured topic rows (`.content-row-divided`).
- **Member Instances**:
  - `About & Disclaimer Card (.content-card)` (Geodetic Help Panel)
  - `Official Engine Replication Benchmark Card (.content-card)` (Geodetic Help Panel)
  - `Reference Systems & Datums Card (.content-card)` (Geodetic Help Panel)
  - `Formats & Input Conventions Card (.content-card)` (Geodetic Help Panel)
  - `Report & Feedback Card (.content-card)` (Feedback Drawer)
  - `Human Verification Card (.content-card)` (Feedback Drawer)
  - `Active Grid Specification Card (.content-card)` (Legacy Stats View)

### Archetype 9.2: Callout & Notice Boxes
*Universally shared recessed notice boxes for disclaimers, acceptance criteria, operational instructions, and origin notices.*
- **Styling Standard**: Shared callout container (`.callout-box`, `.fhp-callout-box`, `.fgp-callout-box`, `.ffp-template-callout`), recessed background (`--bg-surface-input`), 1px medium border (`--border-medium`), border radius `--radius-sm` (6px), padding `10px 12px`, flex-column layout with gap `4px`. Features a bold sentence-case title (`.callout-title`, `11.5px`, `600` weight, `--text-main`, strictly no trailing colon) and body text (`.callout-body`, `11.5px`, line-height `1.45`, `--text-muted`).
- **Member Instances**:
  - `About Section Disclaimer Box (.callout-box)` (Geodetic Help Panel) — *Title: "Disclaimer"*
  - `Benchmark Acceptance Criteria Box (.callout-box)` (Geodetic Help Panel) — *Title: "Acceptance Criteria"*
  - `Binary Grid Origin Callout Box (.callout-box)` (Grid Metadata Panel) — *Title: "Binary Grid Origin"*
  - `Feedback Description Template Callout Box (.callout-box)` (Feedback Drawer) — *Title: "Detailed Description"*
  - `Accuracy & Distortion Note (.callout-box.accent-primary)` (Legacy Stats View) — *Title: "Accuracy & Limitations"*
