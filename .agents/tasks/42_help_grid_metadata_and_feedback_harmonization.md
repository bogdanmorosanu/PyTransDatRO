# Harmonization of Geodetic Help, Grid Metadata, and Feedback Panels

## Objectives
- [x] Modernize the **Geodetic Help Panel** (`GeodeticHelpPanel.ts`) with Section Content Cards (Archetype 9.1), Two-Tier Headings, Content Topic Items (Archetype 4.5), Notice / Callout Boxes (Archetype 9.2), structured benchmark comparison table with separated planimetric and elevation columns, and angle format preview chips (Archetype 5.4).
- [x] Modernize the **Grid Metadata Panel** (`GridMetadataPanel.ts`) with Notice / Callout Box (Archetype 9.2) for the authoritative binary grid origin, and Selectable Readout chips (Archetype 5.4) for `Size` and `SHA-256`.
- [x] Modernize the **Report and Feedback Modal** (`FeedbackPanel.ts`) with Section Content Cards (Archetype 9.1), Two-Tier Headings, dedicated template guidance callout box, and a separate **Human Verification** section card with twin-field algebraic layout (`[ math ] = [ result ]`) and non-bulky 28px height.
- [x] Achieve complete **icon-free typography** across the help ecosystem: remove `iconSvg` from Header Help dropdown items (`Header.ts`) and eliminate decorative SVG icons from panel title bars.
- [x] Standardize window controls across all floating panels with minimize glyph `—`, restore glyph `□`, and reactive localized tooltips (`t('common.restore')` / `t('common.minimize')`).
- [x] Strictly eliminate trailing colons (`:`) across all informational and form field labels (Archetype 4.3) across English and Romanian dictionaries.
- [x] Eliminate 100% of inline style attributes in favor of modular CSS tokens in `webapp/src/style/panels.css`.
- [x] Synchronize living architecture documentation in `docs/webapp/archetypes.md` and `docs/webapp/actionable_inventory.md`.

---

## Implementation Details

1. **Application Header Dropdown (`webapp/src/components/Header.ts`)**:
   - Removed `iconSvg` properties from `#menu-help-app-btn`, `#menu-grid-info-btn`, and `#menu-feedback-btn`.
   - The Help dropdown menu items render cleanly as pure typography, aligning with external resource links (API Docs, ROMPOS, RomGEO).

2. **Geodetic Help Panel Modernization (`webapp/src/components/GeodeticHelpPanel.ts`)**:
   - Structured content into 4 discrete **Section Content Cards** ([Archetype 9.1](../../docs/webapp/archetypes.md#archetype-91-section-content-cards)):
     - *About & Disclaimer*: Includes non-official ANCPI disclaimer wrapped in `.fhp-callout-box` (Archetype 9.2).
     - *Official Engine Replication Benchmark*: Evaluated dynamically against national control stations; displays separated Planimetric (`"`) and Elevation (`mm`) residuals using `.fhp-metric-chip` readouts and status pills without inequality clutter.
     - *Reference Systems & Datums*: Stereo 70, ETRS89, and Altimetry topics using `.fhp-topic-item` (Archetype 4.5) with dashed dividers.
     - *Formats & Input Conventions*: Angular format preview chips (`46°00'00.00000"`, `46.00000000`, `0.8028514559`), coordinate ordering, and 2D/3D auto-detection.
   - Enforced the **Two-Tier Heading Structure**: 13.5px Section Header with solid border-bottom underline, and 11.5px Topic Subheadings with dashed dividers. Trailing colons strictly eliminated.
   - Removed decorative SVG icon from header title cluster.

3. **Grid Metadata Panel Refinements (`webapp/src/components/GridMetadataPanel.ts`)**:
   - Converted the authoritative file origin notice into `.fhp-callout-box.fgp-callout-box` (Archetype 9.2).
   - Displayed `Size` and `SHA-256` as copyable, selectable monospace chips (`.fhp-metric-chip.mono`) without trailing colons.
   - Removed decorative `<svg class="fgp-header-icon">` from normal and minimized headers.
   - Localized minimize/restore button tooltips with `t('common.restore')` and `t('common.minimize')`.

4. **Report & Feedback Modal Refinements (`webapp/src/components/FeedbackPanel.ts`)**:
   - Structured the modal body into **two independent Section Content Cards**:
     - *Card 1 (Feedback Card)*: Report Type selector, template guidance in `.fhp-callout-box`, multiline monospace description textarea, and auxiliary action buttons (Reset Template, Copy Report, Open GitHub Issue).
     - *Card 2 (Human Verification Card)*: Dedicated card with its own Two-Tier Heading (`Human Verification` / `Verificare Umană`).
   - Algebraic Equation Twin-Field Layout:
     - Math problem box (`.ffp-math-expr`) styled with identical height (`28px`), font size (`12px`), regular font weight (`500`), subtle border, and radius as the result input field.
     - Hardcoded standalone `=` symbol (`.ffp-math-equals`) placed between the math display box and the answer input field.
     - Result field (`.ffp-math-input`) reduced to a non-bulky, compact `28px` height (matching the auxiliary buttons) and `68px` width with centered text.
   - Removed decorative `<svg class="fgp-header-icon">` and standardized minimize character to `—`.
   - Eliminated all inline style attributes, introducing `.ffp-card`, `.ffp-verify-card`, `.ffp-field-group`, `.ffp-message-group`, `.ffp-template-callout`, and `.ffp-honeypot`.

5. **Design System & Modular CSS (`webapp/src/style/panels.css`)**:
   - Added classes: `.fgp-callout-box`, `.fgp-meta-chips-row`, `.fgp-meta-chip-item`, `.fgp-meta-label`.
   - Added classes: `.ffp-card`, `.ffp-verify-card`, `.ffp-field-group`, `.ffp-message-group`, `.ffp-template-callout`, `.ffp-verify-content`, `.ffp-math-cluster`, `.ffp-math-expr`, `.ffp-math-equals`, `.ffp-math-input`, `.ffp-honeypot`.
   - Tuned `.ffp-body` with `gap: 12px` for multi-card spacing.

6. **Internationalization (`webapp/src/i18n/en.ts`, `webapp/src/i18n/ro.ts`)**:
   - Added `common.minimize` and `common.restore`.
   - In `gridMetadata`: added `originTitle`, removed colons from `sizeLabel` and `sha256Label`.
   - In `feedback`: added `verificationSubtitle`, removed colon from `mathLabel`.
   - In `geodeticHelp`: added full English and Romanian dictionaries for benchmark table, acceptance criteria, datum notes, and angle chips.

7. **Living Documentation Synchronization**:
   - [`docs/webapp/archetypes.md`](../../docs/webapp/archetypes.md): Documented Archetype 3.1 icon-free menu items, Archetype 4.3 colon-free labels, Archetype 5.4 selectable info chips, Archetype 9.1 section cards (including Report & Feedback and Human Verification), and Archetype 9.2 notice boxes.
   - [`docs/webapp/actionable_inventory.md`](../../docs/webapp/actionable_inventory.md): Updated Help dropdown and floating modal entries.

---

## Verification & Testing

- **Vitest Unit Suite (`webapp/src/`)**:
  - `npm --prefix webapp test`: **80 passed (80)** across 6 test files (`pointParser`, `dmsParser`, `gridBounds`, `exporters`, `ui`, `numberFormat`).
- **Production Bundle & Type Checking**:
  - `npm --prefix webapp run build`: `tsc && vite build` completed in **445ms with 0 errors**.
- **Python Linting**:
  - `python -m ruff check --no-cache pytransdatro api tests tools`: **All checks passed**.
- **Core Geodetic Precision Suite**:
  - `python -m pytest tests/test_trans_ro.py`: **54 passed (54)** with byte-for-byte ANCPI benchmark compliance ($< 10^{-9}$ tolerance).
- **Parallel Testing Agent Boundary**:
  - Confirmed `webapp/e2e/**` and `docs/testing.md` remained 100% untouched throughout the session.
