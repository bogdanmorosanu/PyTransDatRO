# Active Period Styling and Date Standardization (Computation & Spatial Dashboards)

## Objectives
- [x] Standardize all read-only telemetry dashboard dates to the **UK / Romanian date format (`DD/MM/YY`)** with 2-digit years across the Computation and Spatial subpages to optimize horizontal space.
- [x] Delineate the user date selection cluster from the database recorded active period using a subtle vertical divider line (`.telemetry-deck-vdivider`).
- [x] Group the user selection inputs (calendar date pickers + `Earliest Available` and `Current Date` limit chips) to the left side of the controls deck.
- [x] Group the database effective period to the right side of the deck with the main label **`Active Period`** ([Archetype 4.3](../../docs/webapp/archetypes.md#archetype-43-element--component-labels)) styled consistently with `Sources` and `Metrics`.
- [x] Style the **`Start`** and **`End`** sub-labels of the Active Period with the exact same archetype and typography (`.info-label`, Archetype 4.3) as `Earliest Available` and `Current Date`.
- [x] Replace the button-like info trigger with an **industry-standard hoverable information indicator (GitHub-style)** (`.info-help-indicator`, [Archetype 4.6](../../docs/webapp/archetypes.md#archetype-46-informational-contextual-tooltip-indicators)) utilizing `cursor: help` and transparent background.
- [x] Implement a **theme-reactive contextual tooltip bubble (`.info-tooltip-bubble`)** that strictly honors the Dark/Light theme (`background-color: var(--bg-surface)`, `color: var(--text-main)`), preventing OS/browser native black box tooltips in light mode.
- [x] Rename the Computation page metric label to plural **`Metrics`** (`stats.metricsFilterLabel`) to reflect multi-selection capability while maintaining singular **`Metric`** (`stats.metricFilterLabel`) in the Spatial view.
- [x] Add localized explanatory tooltip copy in English and Romanian clarifying database monthly/quarterly temporal snapping.
- [x] Update Vitest unit tests, Playwright E2E suites (`computation.spec.ts`, `spatial.spec.ts`), and living documentation (`archetypes.md`, `actionable_inventory.md`).

---

## Implementation Details

1. **Date Formatter Utility (`webapp/src/utils/numberFormat.ts`)**:
   - Implemented `formatDisplayDate(isoDateStr: string): string` converting ISO-8601 `YYYY-MM-DD` and `YYYY-MM` formats to UK standard `DD/MM/YY`.
   - Returns input unchanged if format is missing or invalid.
   - Comprehensive unit test coverage added in `webapp/src/utils/numberFormat.test.ts` (15/15 tests passing).

2. **Internationalization Dictionaries (`webapp/src/i18n/en.ts`, `webapp/src/i18n/ro.ts`)**:
   - Added `stats.metricsFilterLabel`: `'Metrics'` (EN) / `'Metrici'` (RO) for the multi-select Computation view.
   - Retained `stats.metricFilterLabel`: `'Metric'` (EN) / `'Metrică'` (RO) for the single-select Spatial view.
   - Added `stats.periodStart`: `'Start'` (EN) / `'Start'` (RO).
   - Added `stats.periodEnd`: `'End'` (EN) / `'Sfârșit'` (RO).
   - Added `stats.activePeriodTooltip`: Localized tooltip string explaining monthly/quarterly telemetry snapping against database intervals.

3. **Modular CSS Styling (`webapp/src/style/telemetry.css`)**:
   - Added `.telemetry-deck-vdivider`: 1px subtle divider (`--border-subtle`), height 24px, hidden on narrow screens (`@media (max-width: 900px)`).
   - Added `.telemetry-date-active-cluster`: Flexible right-aligned container (`margin-left: auto`, `align-items: center`, `gap: 10px`).
   - Styled `.telemetry-active-main-label`: Inline `.form-label` (12px, weight 500, color `--text-muted`).
   - Created `.info-help-indicator`: GitHub-style hoverable SVG indicator (`cursor: help`, transparent background, zero border, color `--text-dim` transitioning to `--text-main` on hover/focus).
   - Created `.info-tooltip-bubble`: Theme-reactive tooltip bubble (`role="tooltip"`) rendered above the indicator with caret arrow (`::after`). Styled using semantic CSS variables (`background-color: var(--bg-surface)`, `color: var(--text-main)`, `border: 1px solid var(--border-medium)`, `box-shadow: var(--shadow-md)`), ensuring a crisp white card in light mode and dark carbon card in dark mode. Suppressed native `title` attribute to prevent OS black tooltips from colliding.

4. **Component Updates (`webapp/src/components/TelemetryDashboard.ts`)**:
   - Computation subpage Row 2 reorganized:
     - Left group: Custom date pickers + `Earliest Available` and `Current Date` limit chips formatted as `DD/MM/YY`.
     - Divider: `.telemetry-deck-vdivider`.
     - Right group: `.telemetry-date-active-cluster` featuring `Active Period` main label, `Start` and `End` sub-labels with readonly chips (`DD/MM/YY`), and `.info-help-indicator`.
   - Computation subpage Row 3: Updated label to `${t('stats.metricsFilterLabel')}` ("Metrics").
   - Spatial subpage Row 2: Reorganized identically with active period right cluster and info indicator. Row 3 retained singular "Metric".

5. **Living Documentation & Specifications**:
   - [`docs/webapp/archetypes.md`](../../docs/webapp/archetypes.md): Added **Archetype 4.6: Informational Contextual Tooltip Indicators** for `.info-help-indicator`; updated Archetype 4.3 (Element & Component Labels) and Archetype 5.4 (Selectable Informational Readouts).
   - [`docs/webapp/actionable_inventory.md`](../../docs/webapp/actionable_inventory.md): Updated Active Period row specification and cross-link.

---

## Verification & Testing
- **Vitest Unit Tests**: `npm --prefix webapp test` -> 6 test suites, 83 tests passing (including `numberFormat.test.ts`).
- **TypeScript & Vite Build**: `npm --prefix webapp run build` -> Clean build in 237ms (0 errors).
- **Python Quality Gate**: `ruff check` -> All checks passed; `pytest tests/test_trans_ro.py` -> 54/54 core geodetic tests passed ($< 10^{-9}$ tolerance).
- **Playwright E2E Browser Suite**:
  - `e2e/computation.spec.ts`: Passed (1/1) - validated UK `DD/MM/YY` regex, `Active Period` main label, `Start` / `End` chips, info indicator, and theme-reactive tooltip bubble.
  - `e2e/spatial.spec.ts`: Passed (1/1) - validated UK `DD/MM/YY` regex, `Active Period` main label, `Start` / `End` chips, info indicator, and theme-reactive tooltip bubble.
