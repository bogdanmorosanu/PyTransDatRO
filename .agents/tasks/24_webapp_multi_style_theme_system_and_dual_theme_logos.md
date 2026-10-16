# Task 24: Web Application Design System Overhaul — GitHub Primer Palette, Subtle Surface + Border Primary Button, Enhanced Toggle ON/OFF States, Interactive Computation Legend Color Picker & Dual-Theme Logo

- **Date**: 2026-04-16
- **Status**: Completed
- **Scope**: `webapp/`, `docs/web_app.md`

---

## 1. Objective & Context

Following the structural CSS and UI component refactoring of the web application, the user requested a comprehensive visual design overhaul and iterative refinement (`Iteration 1` exploration $\to$ `Iteration 2` logo/menu consolidation $\to$ `Iteration 3` GitHub Primer palette, button color selection, toggle ON/OFF contrast, and interactive Computation chart legend color picker):
1. **GitHub Primer Neutral Dual-Theme Palette & Subtle Surface + Border (`#8`) Primary Action Button**:
   - Adopt the GitHub Primer neutral Dark (`#0d1117` app canvas, `#161b22` card surface, `#21262d` elevated controls) and Light (`#f6f8fa` app canvas, `#ffffff` card surface, `#f3f4f6` elevated controls) theme tokens, eliminating dark navy-blue backgrounds (`#090d16` / `#0f172a` / `#1e293b`).
   - Lock in Option `#8` (**Subtle Surface + Border**) as the permanent primary action button style (`--btn-primary-bg: #1f2630`, `--btn-primary-hover: #28313e`, `--btn-primary-text: #e6edf3`, `--btn-primary-border: #48586a` in Dark mode; `--btn-primary-bg: #eaeef2`, `--btn-primary-hover: #dde3ea`, `--btn-primary-text: #1f2328`, `--btn-primary-border: #8c959f` in Light mode) and remove the temporary button color selector from the header.
2. **Enhanced ON/OFF State Visibility for Toggle Controls (Transform & Usage Pages)**:
   - Slightly increase ON/OFF visual contrast on the Transform page direction `.sliding-toggle` glider and `.segmented-control` coordinate order pills.
   - Add a `6px × 6px` circular ON/OFF status dot (`::before`), elevated `--btn-primary-bg` active surface, and crisp `--btn-primary-border` outline on standalone Usage filter toggles (`REST API`, `Web App`, `Desktop / CLI`, `Points`, `Calls`), while recessed/dimming (`opacity: 0.62`) inactive buttons.
3. **High-Contrast Dual-Theme Computation Chart Palette & Interactive Legend Line Color Picker**:
   - Replace the low-contrast pale yellow (`#fde68a` / `#f59e0b`) used for `Desktop / CLI` (`desktop_cli`) in the Computation timeline chart with Coral Rose (`#fb7185`) / Crimson Rose (`#e11d48`), legible on both Dark and Light backgrounds.
   - Add click-to-change 10-swatch color popovers (`.legend-color-popover`, matching the Web Map layer widget UX) on every line item in the Computation chart legend, persisting custom line colors in `localStorage['pytransdat_timeline_colors']`.
4. **2-Line Lowercase Repository Attribution & External Link Rule**:
   - Keep the GitHub-style left-zone placement right after the `PyTransDatRO` brand block with 2 stacked lines (`powered by` in lowercase on line 1, `PyTransDatRO` on line 2) and zero oblique `↗` arrow.
   - Restrict the oblique external link arrow (`↗`, `.ext-icon`) strictly to external links inside dropdown menus.
5. **Official Dual-Theme Brand Logo (`C2: Map Projection Graticule`)**:
   - Consolidate on **C2: Map Projection Graticule** as the official dual-theme brand logo in `webapp/src/assets/logos.ts`, `webapp/src/assets/logo.svg`, and `webapp/public/favicon.svg`, reactive to Dark and Light themes both inline (`getLogoInlineSvg(36)`) and in the browser tab favicon (`syncFaviconWithTheme`).
6. **Status-Displaying Language Selector Dropdown Menu & `RomGEO Web API` Link**:
   - Convert the header Language control into a disclosure **Dropdown Menu** (`renderDropdownMenu` with globe icon, `RO`/`EN` status label, and `Română (RO)` / `English (EN)` options with active checkmark `✓`), and add **`RomGEO Web API`** (`https://api.romgeo.ro/api/v1/`, `isExternal: true`) with standardized **`RomGEO`** capitalization.

---

## 2. Technical Architecture & Key Changes

### 2.1 GitHub Primer Palette & Subtle Surface + Border (`#8`) Button (`webapp/src/style/variables.css`)
- Configured `:root, :root[data-theme="dark"]` and `:root[data-theme="light"]` with the GitHub Primer neutral palette and Option `#8` (`Subtle Surface + Border`) primary button tokens.
- Removed all temporary `data-btn-color` rules and the `#header-btn-color-select` dropdown from `Header.ts`, `main.ts`, and `layout.css`.

### 2.2 Enhanced ON/OFF Toggle State Visibility (`forms.css`, `telemetry.css`)
- **Transform Page (`forms.css`)**:
  - Styled `.sliding-toggle` and `.segmented-control` with a recessed `--bg-surface-input` track and an elevated `--btn-primary-bg` active pill/glider with a `1px solid var(--btn-primary-border)` border and dimmed (`opacity: 0.65`) inactive text.
- **Usage Page (`telemetry.css`)**:
  - Added a `6px × 6px` circular ON/OFF status dot (`::before`) to `.telemetry-toggle-btn.individual` (solid `--color-primary` dot when `.active`, hollow ring dot and `opacity: 0.62` when inactive).

### 2.3 Dual-Theme Computation Chart Colors & Interactive Legend Color Picker (`TimelineChart.ts`, `TelemetryDashboard.ts`, `telemetry.css`)
- Replaced the pale yellow `desktop_cli` colors in `KNOWN_SOURCE_COLORS` with Coral Rose (`#fb7185`) / Crimson Rose (`#e11d48`) and removed the hardcoded `.chart-total-line.total-calls { stroke: #f8fafc; }` override in `telemetry.css`.
- Added `TIMELINE_COLOR_PRESETS` (10 swatches) and `localStorage['pytransdat_timeline_colors']` persistence (`getLineColor`, `setLineColor`, `getTotalSeriesColors`) in `TimelineChart.ts`.
- Updated `TelemetryDashboard.ts` to render clickable `.interactive-legend-item` elements with `.legend-color-popover` swatch grids in `renderTimelineLegendHtml()` and wired up swatch selection and outside-click dismissal in `attachTimelineEventListeners()`.

### 2.4 Official Brand Logo, Header Attribution, Language Dropdown Menu & `RomGEO Web API`
- Consolidated `webapp/src/assets/logos.ts` around **C2: Map Projection Graticule** (`getLogoInlineSvg(36)`, `getLogoStandaloneSvg('dark' | 'light')`, `syncFaviconWithTheme(theme)`).
- Rendered lowercase `powered by` / `PyTransDatRO` in `.header-repo-link.zone-left-attribution`, the Language Selector Dropdown Menu (`#header-lang-details`), and the `RomGEO Web API` Help menu link.

---

## 3. Verification & Quality Audit

- **Frontend Unit Tests (`npm --prefix webapp test -- --run`)**: **67 / 67 passed**.
- **Frontend Production Build (`npm --prefix webapp run build`)**: **0 TypeScript or Vite bundler errors** (`dist/assets/index-CulpTjGK.css` 93.11 kB, `dist/assets/index-B_MwwnZw.js` 413.69 kB).
- **Python Linting (`ruff check pytransdatro api tests`)**: **0 issues**.
- **Full Regression Suite (`pytest` + ANCPI Benchmark)**: **98 / 98 passed**, 0 numerical regressions ($< 10^{-9}$ tolerance).
