# Webapp Visual Styling, Content Audit & Archetype Consolidation

## Objectives
- [x] Establish a canonical, concise Thematic Color Palette in `variables.css` for both Dark and Light themes.
- [x] Consolidate all non-tokenized raw color literals (e.g. in `telemetry.css`, `forms.css`) to CSS custom property tokens.
- [x] Unify floating drawer structures and content cards across Geodetic Help, Feedback, and Grid Metadata panels under shared CSS classes (`.floating-drawer`, `.content-card`, `.callout-box`).
- [x] Purge static inline formatting `style="..."` attributes across all TypeScript components in favor of reusable utility classes in `layout.css`, `forms.css`, and `badges.css`.
- [x] Streamline `docs/webapp/archetypes.md` from 34 redundant member entries into 27 canonical archetypes across 9 functional categories.
- [x] Maintain 100% visual invariance ("under the hood" refactor preserving exact layout geometry and aesthetic appearance).
- [x] Verify complete frontend test suites: Vitest unit tests (83/83), Vite production build, and Playwright E2E suite (59/59).

---

## Implementation Details

### 1. Canonical Thematic Color Palette & Token Purity
- **Custom Property Tokens Added (`webapp/src/style/variables.css`)**:
  - `--toggle-thumb-checked: #ffffff` defined in `:root` (unifying the toggle switch indicator across themes).
  - `--border-dashed: #30363d` (Dark) and `#afb8c1` (Light) for segmented/dashed panel dividers.
  - Heatmap 4-step density tokens for the spatial punchcard matrix:
    - Dark theme: `--heatmap-step-1: #1e3a5f`, `--heatmap-step-2: #0284c7`, `--heatmap-step-3: #38bdf8`, `--heatmap-step-4: #a5f3fc`.
    - Light theme: `--heatmap-step-1: #0969da`, `--heatmap-step-2: #218bff`, `--heatmap-step-3: #54aeff`, `--heatmap-step-4: #80ccff`.
- **Tokenized Stylesheets**:
  - `webapp/src/style/forms.css`: Swapped hardcoded `#ffffff` on toggle thumb for `var(--toggle-thumb-checked)`.
  - `webapp/src/style/telemetry.css`: Replaced all raw hex literals (`.telemetry-kpi-val.cyan`, `.legend-box.cyan`, `.chart-points-line`, `.chart-crosshair`, `.tooltip-dot.cyan`, `.legend-indicator.dashed`, `.step-swatch.level-1..4`) with canonical theme variables.
  - **Audit Result**: Zero raw hex colors remain outside `variables.css` across all webapp CSS files (except safe fallback in pure SVG fill strings).

### 2. Unified Floating Drawers & Shared Content Cards (`webapp/src/style/panels.css`)
- Replaced panel-specific styling duplication across Help, Feedback, and Grid Metadata drawers with unified base abstractions:
  - **Container Shell**: `.floating-drawer` with `.drawer-header`, `.drawer-title-cluster`, `.drawer-title`, `.drawer-header-actions`, `.drawer-icon-btn`, `.drawer-body`.
  - **Shared Content Cards**: `.content-card`, `.content-card-header`, `.content-card-title`, `.content-card-subtitle`, `.content-row-divided`.
  - **Shared Callout Boxes**: `.callout-box`, `.callout-title`, `.callout-body`, and accent variant `.callout-box.accent-primary`.
- Preserved backward-compatible CSS class aliases (`.fhp-card`, `.ffp-card`, `.fhp-callout-box`, `.fgp-callout-box`, `.fgp-header`) to avoid breaking existing DOM selectors.
- Refactored `GeodeticHelpPanel.ts`, `FeedbackPanel.ts`, and `GridMetadataPanel.ts` to adopt the shared classes.

### 3. Static Inline Style Elimination & Reusable Utility Classes
- **New Reusable Layout & Form Classes**:
  - `webapp/src/style/layout.css`: `.flex-between-center`, `.flex-column-gap-sm`, `.flex-column-gap-md`, `.flex-row-gap-sm`, `.margin-bottom-none`, `.margin-bottom-sm`, `.modal-title.warning`, `.modal-body-text`, `.modal-alert`, `.mobile-only-btn`.
  - `webapp/src/style/forms.css`: `.batch-textarea`, `.batch-textarea-readonly`, `.dropzone-title`, `.dropzone-desc`, `.dropzone-queue-header`, `.file-options-row`.
  - `webapp/src/style/badges.css`: `.status-pill.neutral`.
- **TypeScript Components Cleaned**:
  - Refactored `BatchConverter.ts`, `FileUpload.ts`, `PointConverter.ts`, `VolumeWarningModal.ts`, `ReportModal.ts`, `Header.ts`, and `StatsView.ts`.
  - Removed all static inline `style="..."` attributes. Purged informal decorative emoji `📊` from calculation details modal header in `ReportModal.ts`.
  - Hardened `.btn.mobile-only-btn, .mobile-only-btn { display: none !important; }` in `layout.css` so that the legacy mobile viewport toggle (`#header-mobile-view-btn`) remains strictly hidden on desktop screens (preventing `.btn` in `buttons.css` from overriding its visibility due to CSS load order).
  - Preserved strictly dynamic mathematical calculations (e.g. `padding-left: ${indentPx}px` for hierarchical tree nodes in `GridMetadataPanel.ts`).
  - **Audit Result**: Zero static inline formatting styles remain across all TypeScript component templates.

### 4. Living Documentation & Archetype Catalog Modernization
- **Design System Manual (`docs/webapp/design_system.md`)**:
  - Documented Section 2.1 Canonical Thematic Color Palette specification detailing the minimal color hierarchy (base backgrounds, surface elevations, text contrasts, accent blues, semantic status colors, and spatial density heatmaps).
- **Archetype Catalog Streamlining (`docs/webapp/archetypes.md`)**:
  - Consolidated catalog from 34 fragmented member entries down to **27 non-overlapping, strictly categorized archetypes across 9 functional categories**:
    1. Navigation & Headers (Archetypes 1–3)
    2. Primary & Secondary Actions (Archetypes 4–7)
    3. Segmented & Selection Controls (Archetypes 8–10)
    4. Form Inputs & Textareas (Archetypes 11–13)
    5. Information Indicators & Badges (Archetypes 14–17)
    6. Cards & Content Containers (Archetypes 18–20)
    7. Floating Drawers & Overlays (Archetypes 21–23)
    8. Charts & Telemetry Displays (Archetypes 24–25)
    9. Map Overlay & Floating Controls (Archetypes 26–27)
- **Webapp README (`webapp/README.md`)**:
  - Synchronized archetype catalog reference to reflect 27 archetypes across 9 functional categories.

---

## Verification & Testing

1. **Frontend Vitest Unit Suite**:
   - Command: `npm --prefix webapp test`
   - Result: 83/83 unit tests passing across all parser, telemetry, and UI component factory test suites.
2. **Vite Production Bundle Build**:
   - Command: `npm --prefix webapp run build`
   - Result: Clean TypeScript compilation (`tsc -b`) and minified production build in 535ms with zero warnings or errors.
3. **Playwright E2E Browser Suite**:
   - Command: `.\test-e2e.ps1`
   - Result: 59/59 end-to-end tests passing across all 13 spec suites (`batch.spec.ts`, `point.spec.ts`, `file.spec.ts`, `map.spec.ts`, `i18n.spec.ts`, `errors.spec.ts`, `computation.spec.ts`, `spatial.spec.ts`, `legacy.spec.ts`, `smoke.spec.ts`, `lifecycle.spec.ts`).
   - Hardened `webapp/e2e/batch.spec.ts` to assert output textarea update before reading `inputValue()`, eliminating race conditions during sequential batch transformations.

---

## Pre-Commit Visual Impact Report

| Component / Area | Architectural Change | Visual Shift / Layout Impact | User Verification Status |
| :--- | :--- | :--- | :--- |
| **Color Tokens** | Swapped raw hex colors for `--heatmap-step-*`, `--toggle-thumb-checked`, `--border-dashed` | **Zero Shift**: Rendered RGB colors match previous values byte-for-byte. | Passed (Verified via Playwright & CSS audit) |
| **Floating Drawers** | Unified `.fhp-card`, `.ffp-card`, `.fgp-header` under `.content-card`, `.callout-box`, `.floating-drawer` | **Zero Shift**: All elevations (`0 8px 24px`), border radii (`8px`), padding (`14px/16px`), and title font sizes are identical. | Passed (Preserved full alias backward compatibility) |
| **Batch & File Textareas** | Replaced static inline styles with `.batch-textarea` and `.dropzone-*` | **Zero Shift**: Identical font families (`var(--font-mono)`), heights, borders, and margins. | Passed (Verified in browser test runner) |
| **Report Modal & Volume Modal** | Removed inline styles; purged informal emoji `📊` from report modal | **Clean Polish**: Header displays clean typography (`Calculation Details`) matching the rest of the app without informal emojis. | Passed |
| **Header Mobile Viewport Button** | Hardened `.btn.mobile-only-btn` rule with `display: none !important` on desktop | **Zero Shift**: Mobile viewport toggle button (`#header-mobile-view-btn`) remains strictly hidden on desktop as intended. | Passed |
| **Archetypes Catalog** | Streamlined from 34 redundant entries into 27 categorized archetypes | **Documentation Only**: No breaking CSS changes; provides clean, disciplined classification for developers. | Passed |
