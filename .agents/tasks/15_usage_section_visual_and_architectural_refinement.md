# Usage Section Visual and Architectural Refinement

## Objectives
- [x] Remove the redundant "Back to Geodetic Studio" button and update page title to "Usage Data".
- [x] Style view switcher tabs with GitHub-style underline indicators ("Computation" and "Spatial") adjacent to the title with a subtle separator, matching the primary `Transform` and `Usage` header navigation.
- [x] Remove the obsolete "Grid Specs" tab completely.
- [x] Replace `noUiSlider` time slider with native themed `<input type="date">` calendar controls matching the Spatial section, including range limits pill (`Earliest Available` • `Current Date`).
- [x] Add global `color-scheme: dark | light` ensuring native browser calendar pickers and form controls render natively in dark or light theme without white glare in dark mode.
- [x] Implement dynamic source discovery via lookup dictionary (`rest_api` -> "REST API", `web_app`/`web_map` -> "Web App", `desktop_cli` -> "Desktop / CLI") with fallback to capitalized database names and filtering out `unknown`.
- [x] Ensure Total line strictly reflects the sum of active selected sources ($\text{Total} \ge \text{each active source}$).
- [x] Enhance Total line symbology with wider dash spacing (`stroke-dasharray="7 5"`), thicker stroke (`2.6px`), and visible 24px repeating dashed sample in the legend.
- [x] Add 3-way Metric selector: `Points`, `Calls`, or `Both`, dynamically adapting curves, left/right Y-axes with explicit title labels (`Calls` on left, `Points` on right), and interactive tooltips.
- [x] Standardize naming to "Points" and "Calls". Neutralize top KPI cards, removing saturated green accents.
- [x] Separate source toggle checkboxes (clean minimalist checkboxes without attached line symbols) from the dedicated graph symbology legend.
- [x] Implement dynamic graph symbology legend synchronized with the visible curves and active metric filter.
- [x] Harmonize Spatial view active period controls into a unified compact preset bar matching Computation.
- [x] Redesign Spatial overview banner into a full-width summary card spanning from Total Points/Calls (left) to Active Cells (right) with the average cell size indicator.
- [x] Optimize quadrant summary cards with vertical stacking (quadrant name above, volume below) and compact adjacent percentage (`1,250 (34.2%)` without redundant "of total" text).
- [x] Verify heatmap axis labels use `°N` and `°E` cartographic notation.
- [x] Relocate scaling mode selector (`Logarithmic` | `Linear`) to the bottom-right of the heatmap canvas viewport alongside the 5-step discrete color ramp.
- [x] Combine grid resolution, average cell area, and quadrant divider notes into a single cohesive footnote below the heatmap.

## Implementation Details

1. **Usage Header & Navigation (`webapp/src/components/TelemetryDashboard.ts`, `webapp/src/style/components.css`)**:
   - Replaced button group with GitHub-style underline tabs (`.telemetry-nav-tabs`, `.telemetry-nav-tab-btn`, `.telemetry-tab-indicator`).
   - Positioned tabs directly to the right of "Usage Data" (`.telemetry-page-title`), separated by a subtle vertical divider (`.header-brand-nav-divider`).
   - Removed "Back to Geodetic Studio" button and obsolete "Grid Specs" tab.

2. **Themed Date Controls & Native Dark Mode Support (`webapp/src/components/TelemetryDashboard.ts`, `webapp/src/style/components.css`, `webapp/src/style/variables.css`)**:
   - Replaced `noUiSlider` dependency in Computation with native `<input type="date" class="telemetry-date-input">` controls identical to Spatial.
   - Set `color-scheme: dark;` under `html[data-theme="dark"]` and `color-scheme: light;` under `html[data-theme="light"]` in `variables.css` and `.telemetry-date-input`, forcing native browser calendar dropdowns to render in dark theme when dark mode is active.
   - Added range limits pill displaying earliest available log date and today's date.

3. **Dynamic Source Resolution & Additive Total Calculation (`webapp/src/components/TelemetryDashboard.ts`, `webapp/src/components/TimelineChart.ts`)**:
   - Introduced `SOURCE_LABEL_LOOKUP` dictionary mapping database keys to canonical UI labels (`REST API`, `Web App`, `Desktop / CLI`), with automated sentence-case fallback for future sources. Completely filtered out `unknown` sources.
   - Refactored `TimelineChart` series and total calculation: the Total curve strictly sums values across `activeSourceKeys` at every date bucket, guaranteeing $\text{Total} \ge \text{each individual source}$.
   - Updated dash pattern to `stroke-dasharray="7 5"` and width to `2.6px`. Configured `TOTAL_SERIES_CONFIG.callsColor` with `var(--text-main, #f8fafc)` for theme responsiveness.

4. **Metric Filter Toggle & Y-Axis Labels (`webapp/src/components/TelemetryDashboard.ts`, `webapp/src/components/TimelineChart.ts`)**:
   - Added `computationMetricFilter` (`'both' | 'points' | 'calls'`).
   - Dynamic SVG rendering only displays left Y-axis and calls curves when calls are active, and right Y-axis and points curves when points are active.
   - Added explicit axis title labels (`Calls` on left, `Points` on right) anchored at the top of the axes.
   - Filtered SVG crosshair tooltip to only show metrics corresponding to the active filter.

5. **Decoupled Source Checkboxes and Synchronized Legend (`webapp/src/components/TelemetryDashboard.ts`)**:
   - Source toggles in the controls deck are rendered as clean, symbology-free checkboxes.
   - Built `renderTimelineLegendHtml()` which dynamically outputs swatches and labels for each active source curve and total curve currently visible on the chart.

6. **Spatial Visual Alignment (`webapp/src/components/TelemetryDashboard.ts`, `webapp/src/components/SpatialMatrixChart.ts`, `webapp/src/style/components.css`)**:
   - Unified active period preset buttons into a single row.
   - Replaced top stats with `.spatial-overview-banner` spanning full width: left column shows Total Points/Calls, right column shows Active Cells with average cell area note.
   - Restructured quadrant summary cards into compact cards displaying quadrant name above and volume below with adjacent percentage (e.g. `1,250 (34.2%)`).
   - Relocated scaling mode selector (`Logarithmic` | `Linear`) to the bottom-right of the heatmap canvas viewport alongside the discrete color ramp.
   - Consolidated footnotes into a single concise note below the heatmap canvas.
   - Verified vertical axis labels render as `°N` and horizontal as `°E`.

8. **Iteration 2 Refinements (`TelemetryDashboard.ts`, `TimelineChart.ts`, `SpatialMatrixChart.ts`, `components.css`)**:
   - **Exact Database Fallback**: Unmapped database sources display their exact string (`sourceName`) without transformation.
   - **Symbology Lightness Convention**: Standardized color lightness across curves and legend swatches: Calls use lighter tints (luminous), Points use darker shades of the respective hue.
   - **Filter Rows Harmonization**:
     - Computation Row 1: Presets aligned left, granularity toggles (`Daily`, `Weekly`, `Monthly`) aligned right (`justify-content: space-between`).
     - Merged Row: SOURCES on left, METRIC on right across both views with labels `SOURCES` and `METRIC`.
     - SOURCES (both views) and METRIC (Computation) are rendered as individual standalone toggle buttons with theme-native neutral borders and backgrounds (`--bg-surface-elevated`), generous 8px spacing, and zero continuous frame enclosure, completely removing blue borders and checkboxes.
     - Toggle button styling `.telemetry-toggle-btn.active` and `.telemetry-toggle-btn.individual.active` uses theme-native high contrast (`--btn-primary-bg` / `--btn-primary-text`), removing electric blue tints.
   - **Strict Date Clamping**: Event listeners clamp input values to database bounds (`startVal < minHistoryDate` clamped to `minHistoryDate`, `endVal > maxHistoryDate` clamped to `maxHistoryDate`), preventing selection of dates prior to `2023-11-01`.
   - **Neutral Computation Chart Axes**: Neutralized left and right Y-axes tick labels and axis titles to match the grey X-axis (`var(--text-muted)`).
   - **2-Line Grouped Legend**: Line 1 groups active sources by client name (`Web App: Calls, Points`), Line 2 displays `Total: Calls, Points` on its own line labeled simply "Total".
   - **Spatial Top Banner**: Retained only the clean numerical Active Cells value in the top overview banner.
   - **Stacked Spatial Heatmap Footnotes**: Split into two distinct stacked notes: simplified grid cell dimensions (`0.1° latitude by 0.1° longitude (~10 km, ~86 km² avg. cell area at 46°N)`) and quadrant divider explanation.
   - **Spatial Matrix Axis Alignment**: Relocated the Latitude axis to the left side with 42px left padding and right-aligned text (`textAlign = 'right'`), eliminating the clipped letter "N".

## Verification & Testing

- **TypeScript Compilation & Production Build**:
  - `npm --prefix webapp run build` (`tsc && vite build`) passed with 0 errors and generated optimized production bundles (`dist/index.html`, `dist/assets/*.css`, `dist/assets/*.js`).
- **Python Backend Test Suite**:
  - `python -m pytest tests/` passed 83/83 tests.
- **Server Health**:
  - FastAPI server (`http://127.0.0.1:8000`) and Vite development server (`http://localhost:5173`) active and responding healthy.
