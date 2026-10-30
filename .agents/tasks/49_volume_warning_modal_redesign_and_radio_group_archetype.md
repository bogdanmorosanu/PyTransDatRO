# Task 49: High-Volume Point Warning Modal Redesign & Archetype 2.5 Radio Group

## Objectives
- [x] Standardize the High-Volume Point Warning modal window header to match the neutral styling of Help and Grid Metadata panels.
- [x] Remove the orange warning styling from the modal title (`.modal-title`), keeping warning colors restricted strictly to callout/message components.
- [x] Eliminate the 3 multi-button layout with long button labels, replacing it with an industry-standard radio selection group with a single confirm button.
- [x] Create a new reusable UI factory control `renderRadioGroup` in `webapp/src/ui/radioGroup.ts`.
- [x] Formalize and catalog **Archetype 2.5: Radio Box / Selection Card Group (Single-Choice Radio Cards)** in living documentation.
- [x] Update Playwright E2E browser tests and add comprehensive Vitest unit tests for the new UI factory.
- [x] Synchronize living documentation in `docs/webapp/` and `webapp/README.md`.

## Implementation Details

1. **Reusable Radio Card / Selection Group UI Factory (`webapp/src/ui/radioGroup.ts`, `webapp/src/ui/index.ts`)**:
   - Implemented `renderRadioGroup` supporting accessible, keyboard-navigable single-choice radio cards.
   - Built custom `.radio-circle` indicator with active center dot (`--color-primary-light`), title, optional monospace count badge (`.radio-badge`), and descriptive subtitle text (`.radio-desc`).
   - Integrated full keyboard navigation (`ArrowDown`/`ArrowRight` to increment, `ArrowUp`/`ArrowLeft` to decrement, cycling through choices) with native ARIA `role="radiogroup"` and `role="radio"`.
   - Exported `renderRadioGroup`, `RadioOption`, and `RadioGroupProps` through `webapp/src/ui/index.ts`.

2. **Modular CSS Styling (`webapp/src/style/forms.css`, `webapp/src/style/layout.css`, `webapp/src/style/panels.css`)**:
   - Added styles for `.radio-group`, `.radio-card`, `.radio-input`, `.radio-circle`, `.radio-content`, and `.radio-badge` in `forms.css`.
   - Standardized `.modal-header` across all dialogs to use `background-color: var(--bg-surface-elevated)`, uniform padding (`14px 20px`), and subtle divider border (`border-bottom: 1px solid var(--border-subtle)`), matching `.drawer-header`.
   - Neutralized `.modal-title` font size (`14.5px`, `700` weight) and removed orange color inheritance from `.modal-title.warning`.
   - Added `.callout-box.accent-warning` with amber left border indicator (`--color-warning`) in `panels.css`.
   - Refined `.modal-dialog.modal-dialog-compact` to `max-width: 460px` for balanced card ergonomics.

3. **High-Volume Point Warning Modal Redesign (`webapp/src/components/VolumeWarningModal.ts`)**:
   - Replaced raw action buttons with `renderRadioGroup` pre-selected to `'all'` (Render all points).
   - Rendered Archetype 9.2 `.callout-box.accent-warning` with *"Performance Notice"* (*"Notă de performanță"*) callout header and volume description.
   - Reduced modal footer to a single primary action button: `#vw-confirm-btn` (*"Continue"* / *"Continuă"*).
   - Standardized dismissive actions (clicking header `✕`, clicking backdrop, or pressing `Escape`) to default to skip map rendering.

4. **Localization (`webapp/src/i18n/en.ts`, `webapp/src/i18n/ro.ts`)**:
   - Added `calloutTitle`, `promptLabel`, `optAll`, `optAllDesc`, `optSample`, `optSampleDesc`, `optSkip`, `optSkipDesc`, and `btnConfirm` in English and Romanian.
   - Preserved legacy translation keys for backwards-compatibility.

5. **Living Design System & Archetypes Catalog (`docs/webapp/archetypes.md`, `docs/webapp/actionable_inventory.md`, `webapp/README.md`)**:
   - Added **Archetype 2.5: Radio Box / Selection Card Group (Single-Choice Radio Cards)** under Category 2 (Selection & Mode Switchers).
   - Updated archetype count to 28 across 9 functional categories.
   - Updated Archetype 6.2 member instances from legacy volume action buttons to `Volume Warning Confirmation Button (#vw-confirm-btn)`.
   - Synchronized actionable inventory with updated radio card triggers and `README.md` tree structure.

## Verification & Testing

- **Vitest Unit Tests (`webapp/src/ui/ui.test.ts`)**:
  - Added test suite for `renderRadioGroup` covering card rendering, initial selection, change event dispatching, and badge formatting.
  - Test suite: 6 test files, 84 tests passing.
- **Playwright E2E Tests (`webapp/e2e/file.spec.ts`, full suite)**:
  - Updated `file.spec.ts` line 316-320 to select `#vw-opt-sample` and confirm via `#vw-confirm-btn`.
  - All 8 tests in `file.spec.ts` passed in 15.0s.
  - Verified batch test suite (`e2e/batch.spec.ts`: 9 tests passed).
- **Core Geodetic Math (`tests/test_trans_ro.py`)**:
  - Full suite passed: 54 tests passing in 3.31s with zero numerical regressions.
- **Frontend Production Build**:
  - `npm --prefix webapp run build` compiled cleanly (0 TypeScript/Vite errors).
