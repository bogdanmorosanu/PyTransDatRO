# Web Application & Map Design Rules

These rules govern all UI, UX, styling, and map-related development in `webapp/`. All AI agents working on frontend components or web map interactions must strictly adhere to these invariants.

For the full design system specification and component archetypes, see [docs/webapp/design_system.md](../../docs/webapp/design_system.md), [docs/webapp/archetypes.md](../../docs/webapp/archetypes.md), and [docs/webapp/actionable_inventory.md](../../docs/webapp/actionable_inventory.md).

---

## 1. Core Purpose & Aesthetics
- **Singular Purpose**: The application exists strictly for geodetic coordinate transformation between Stereo 70 and ETRS89. Do not introduce extraneous non-essential validations or auxiliary widgets that detract from this purpose.
- **Zero Decorative Emojis**: Tab labels, action buttons, and segmented options must never contain emojis (no ⚡, ⬍, ⬄, 🧭, etc.). Clean typography communicates technical authority.
- **Micro-Headers**: Use subtle, sentence-case micro-headers (e.g. `Input coordinates (or click on map)`, `Computed coordinates`, `Download As:`).
- **Workbench Panel Width**: Strictly maintain `480px` for the left workbench panel (`.workbench-layout`).
- **No Ad-Hoc Styling**: Never style individual controls or inputs separately with bespoke CSS rules. All controls (segmented toggles, selects, copy badges, inputs) must use unified classes and derive from CSS variables in `variables.css`.
- **Modular Stylesheet Architecture**: Styles must be cleanly divided by category (`forms.css`, `buttons.css`, `badges.css`, `layout.css`, `telemetry.css`) rather than bundled into a single monolithic file.
- **UI Helper Factories**: Use vanilla TS helper functions in `webapp/src/ui/` for rendering standard controls to guarantee identical DOM structure and styling across tabs.

---

## 2. Navigation & Layout Continuity
- **Seamless Folder Bookmark Tabs**: The `.tabs-bar` must physically connect with the active card container below with zero background gap. The active tab merges seamlessly with the top border of the card surface.
- **Prefix Label**: Navigation tabs must include the subtle prefix `Transform:` (`Transformă:`) before concise tab buttons (`Point`, `Text`, `Files`).
- **Header Direct Links**: `Grid Info` and `Geodetic Help` belong in the top navbar (`Header.ts`), keeping input tabs focused strictly on transformation modes.

---

## 3. Shared Transformation Configuration Bar
- Positioned at the top of the workbench panel, globally controlling all tabs (`Point`, `Text`, `Files`):
  1. **Tier 1 (Direction Toggle)**: Sliding directional toggle with sleek 6px rounded-rectangle glider (`var(--radius-sm)`). Default direction is `ETRS89 ➔ Stereo 70` (first button).
  2. **Tier 2 (50/50 Split)**: `Angular unit: [select]` (Degrees, DMS, Radians) | `Coordinate order: [ NE(H) | EN(H) ]`. Must use sentence case.
  3. **Tier 3 (Angle Format Example)**: Copyable badge updating reactively with chosen unit (`45.84123450`, `45°50'28.44420"`, `0.7997813210`). Clicking badge copies strictly the numeric angle value without label text.

---

## 4. Buttons, Colors & Metric Isolation
- **Primary CTA Anchor (`.btn-primary`)**: Must use Option #8 (**Subtle Surface + Border**) theme-native styling:
  - Dark Mode: Subtle charcoal surface (`--btn-primary-bg: #1f2630`, hover: `#28313e`) with medium slate border (`--btn-primary-border: #48586a`) and crisp light text (`--btn-primary-text: #e6edf3`).
  - Light Mode: Subtle light grey surface (`--btn-primary-bg: #eaeef2`, hover: `#dde3ea`) with defined border (`--btn-primary-border: #8c959f`) and dark text (`--btn-primary-text: #1f2328`).
  - Label: Clean text `Transform` / `Transformă` (no emojis, no icons).
- **Metric Color Isolation**:
  - Quasigeoid Height ($\zeta = h - H$) is a physical measurement, not a status. It must **never** be colored green, red, or yellow.
  - Render Quasigeoid Height in its own dedicated box on its own row, styled with `--color-metric-text` (cyan-blue) and `--color-metric-bg`.
  - Green (`#10b981`), yellow (`#f59e0b`), and red (`#ef4444`) are strictly reserved for system health, warnings, and errors.

---

## 5. Coordinate Inputs, Outputs & Ergonomics
- **Asymmetrical 3-Column Layout**: Northing/Lat (`1.35fr`), Easting/Lon (`1.35fr`), and Altitude/Height (`0.8fr`) so long DMS strings never wrap or clip.
- **Dynamic Coordinate Order Swapping**: Toggling between `NE(H)` and `EN(H)` must dynamically swap Column 1 and Column 2 across inputs, outputs, and clipboard exports.
- **Pure Numeric Outputs**: Output text boxes must contain only pure numbers (strip `m`, `°`, `rad`).
- **Persistent Output Layout**: Output boxes must persistently display `—` when cleared or invalid, preventing vertical layout shift.
- **Focus Auto-Select**: All coordinate inputs must execute `input.select()` on focus.
- **Unspaced DMS Formatting**: Always format DMS strings without spaces after `°` and `'` (`DD°MM'SS.SSSSS"`, e.g. `46°00'00.00000"`).
- **Single-DMS Paste Guard**: Pasting a full DMS string into the first input box must not split across multiple boxes.
- **Auto 2D/3D Detection**: 2D when 3rd coordinate is empty; 3D when provided.

---

## 6. Batch Text & File Processing
- **Strict Numeric Validation**: Use strict decimal validation (`isValidDecimal`) to prevent typos like `5001a50.000` from being partially parsed into coordinates.
- **Strict 1-to-1 Line Correspondence**: Line $N$ in corresponds strictly to line $N$ out. Empty lines, comments, and invalid rows maintain their exact sequence index.
- **Human 1-Based Indexing**: Refer to unnamed points as `Point 1`, `Point 2`, etc.
- **Out-of-Grid Suppression**: Suppress approximate coordinates for out-of-grid points; emit only `[Point_ID], # Out of grid`.
- **Bulk File Queue**: Use interactive queue (add $\to$ remove $\to$ manual Transform $\to$ cumulative quality report $\to$ download). Never show raw 5-line text previews.

---

## 7. Web GIS Map & Viewport Invariants
- **3 Independent LayerGroups**: Split map markers into `pointLayer`, `textLayer`, and `filesLayer`, independently toggleable in Leaflet layer control.
- **Decoupled Viewports**: Coordinate transformations must **never** automatically pan or zoom the map.
- **Dual Navigation Actions**: Bottom-right action footer of each tab must provide:
  - `Center on Map`: Pans to centroid of valid in-grid coordinates, preserving current zoom level.
  - `Zoom to Extents`: Pans and fits bounds (`fitBounds`) to encompass all points, or zooms to level 14 for a single point.
- **Boundary-Safe Guard**: Centering/zooming is strictly bounded by Romania's envelope (`[43.5, 20.2]` to `[48.3, 29.8]` and `gridBounds.ts`). Never pan outside Romania on out-of-grid points.
- **Smart Map Click**:
  - `ETRS89ToStereo70`: Populates clicked lat/lon formatted in active angular units.
  - `Stereo 70ToETRS89`: Reverse-calculates clicked point into Stereo 70 projection meters.

---

## 8. Build Verification
Always run `npm --prefix webapp run build` before completing any frontend task to guarantee zero TypeScript or bundler errors.
