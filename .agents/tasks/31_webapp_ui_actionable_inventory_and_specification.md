# Task 31: Web Application Actionable Elements Inventory & UI Living Specification

## Objectives
- [x] Create a comprehensive living UI Actionable Elements Inventory (`docs/ui_actionable_inventory.md`) detailing all actionable controls, reactive surfaces, and event mechanics across the web application.
- [x] Classify UI elements by interaction taxonomy (Segmented Switch, Dropdown Menu, Action CTA Button, Numeric/Text Input, Map Control, Floating Window, Aggregator Card).
- [x] Document expected behaviors, target scopes, state mutations, and side effects for Top Navigation, Transformation Workbench (Point, Batch Text, File Upload), GIS Map Canvas, Floating Diagnostic Panels, and Usage Telemetry Dashboard.
- [x] Conduct an exhaustive code-level scan across `webapp/src/` to identify unrecorded functional capabilities, behavioral divergences, and UX/testability recommendations.
- [x] Reconcile the inventory with user clarifications and finalize the authoritative specification for subsequent automated E2E test generation.

## Implementation Details

1. **Authoritative Living UI Inventory (`docs/ui_actionable_inventory.md`)**:
   - Structured into 7 functional sections:
     - **Section 1: Global Taxonomy & Interaction Classification**: Defines standard UX mechanics and UI factory patterns (`renderSegmentedControl`, `renderDropdownMenu`, `renderSelect`, `renderCopyBadge`, `renderCopyButton`).
     - **Section 2: Top Navigation Bar & Global State (`Header.ts`, `main.ts`)**: Documents brand actions, dynamic release version badge, REST API health checking & recovery polling, repository attribution, primary tab navigation, Help menu dropdown, locale switching, dark/light theme switching, mobile viewport switcher, URL query parameter deep-linking (`?op=`, `?c1=`, `?c2=`, `?c3=`, `?unit=`, `?order=`), and global `Escape` keyboard dismissal.
     - **Section 3: Transformation Workbench (`PointConverter.ts`, `BatchConverter.ts`, `FileUpload.ts`)**:
       - *Shared Configuration Bar*: Direction toggle, angular unit dropdown, coordinate order switch, copyable angle format example badge. Defines the **deferred processing rule**: direction/unit/order adjustments apply immediately to Point mode, while for Batch Text and File Uploads, existing input remains untouched until the user clicks `Transform`.
       - *Single Point Mode*: Auto-sanitized inputs (quotes, `m/metri`, compass letters stripped), `input.select()` focus behavior, `Enter` key calculation trigger, smart delimiter-split paste (strictly manual calculation trigger), read-only selectable outputs, cyan-isolated Quasigeoid height ($\zeta$), non-selectable visual reporting labels, comma-separated copy button (`val1, val2, val3`), and bounded Center/Zoom map actions.
       - *Batch Text Mode*: Resizable textarea, auto-detected delimiters (comma, semicolon, tab, multiple spaces), 1-to-1 line correspondence (comments `#` and blank rows preserved), mixed 2D/3D and named/unnamed coordinate parsing, input locking during calculation, quality report summary badges, CSV line-by-line replica export, clean vector filtering for GeoJSON FeatureCollections and AutoCAD DXF POINT entities.
       - *File Upload Queue*: Staged multi-file queue allowing setting adjustments between additions, cumulative map rendering across all processed files, volume warning guard (>1,000 pts), individual file downloads, staggered bulk download (`Download All`), and queue clearing.
     - **Section 4: Web GIS Map Canvas (`MapView.ts`)**:
       - Direction-aware smart map click (the only action that automatically computes output coordinates and updates the map pin), marker layer color permanence (point marker strictly preserves its layer color swatch without flashing red on selection or error), CAD marker label zoom visibility threshold ($\ge 15$), cluster spiderfying at max zoom, unified geodetic attribute popups with in-place copy feedback, top-left navigation dock (Zoom In/Out, Zoom to Romania, single-distance measurement tool with reset and `Escape` cancel, fullscreen toggle), expandable legend widget, segmented graphic scale bar, bottom-left layers drawer (OSM, Carto Dark/Positron coupled to theme, Satellite, customizable layer swatches), and live cursor coordinate badge.
     - **Section 5: Floating Windows & Diagnostic Panels (`GridMetadataPanel.ts`, `GeodeticHelpPanel.ts`, `FeedbackPanel.ts`)**:
       - Draggable headers, bottom-right corner resizing, minimize/restore mechanics, and close buttons.
       - Grid Metadata live search filter, `Expand All`, `Collapse All`, and raw JSON clipboard export.
       - Geodetic Help live benchmark comparison against ANCPI test vectors.
       - Feedback modal with contextual template pre-population, anti-bot math challenges, and direct `Open as GitHub Issue` generation.
     - **Section 6: Usage & Telemetry Dashboard (`TelemetryDashboard.ts`)**:
       - Sub-navigation tabs (Computation, Spatial, Legacy).
       - Computation timeline: KPI summary cards, time presets with date picker synchronization, granularity grouping, dynamic source toggles, dual-metric toggles (`Points` AND `Calls` simultaneously), and interactive legend color pickers.
       - Spatial grid matrix: 10km geodetic cell density canvas over Romania's bounds, quadrant share cards (NW, NE, SW, SE), **mutually exclusive** metric selector (`Points` OR `Calls`), and logarithmic/linear color scaling.
       - Legacy archive: Historical transformation volumes (2011–2026), fixed archive boundaries (`2011-01-09` to `2026-10-14`), legacy sources (`Web App`, `REST API`), legacy operations (`Stereo 70`, `Stereo 30`), and cumulative total overlay.
     - **Section 7: Coverage Gap Analysis & Recommended Action Items**: Detailed mapping of every element to Playwright tests in `webapp/e2e/`, highlighting test holes for automated suite generation.

2. **Codebase Functional Audit (`scratch/ui_inventory_audit_report.md`)**:
   - Performed deep inspection of all UI components, map controllers, and utilities in `webapp/src/`.
   - Identified and resolved 16 previously unrecorded features, 5 behavioral divergences, and 8 UX/testability recommendations.

## Verification & Testing
- **Specification Coherence**: Reconciled all user requirements with living code mechanics in `webapp/src/`.
- **TypeScript & Build Verification**: Verified that `npm --prefix webapp run build` builds cleanly with zero TypeScript errors or broken imports.
- **Test Accompaniment**: All documented controls are linked to existing Playwright specs in `webapp/e2e/` (`smoke.spec.ts`, `point.spec.ts`, `batch.spec.ts`, `file.spec.ts`, `map.spec.ts`, `legacy.spec.ts`, `i18n.spec.ts`, `errors.spec.ts`).
