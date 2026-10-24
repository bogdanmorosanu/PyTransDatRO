# Mobile Compatibility & Dockable Workbench Architecture

## Objectives
- [x] Implement desktop dockable workbench with vertical divider toggle (`#workbench-dock-toggle`) and keyboard shortcuts (`Ctrl+\` / `Alt+W`).
- [x] Implement Google Maps-style mobile bottom sheet (`.utility-panel`) with 3 states (`sheet-peek`, `sheet-half`, `sheet-full`), touch swipe gestures, and input focus auto-elevation.
- [x] Preserve standard `<input type="text">` across `PointConverter.ts` to support non-numeric DMS notation, degree symbols (`°`, `'`, `"`), spaces, and direction letters without mobile keyboard restriction.
- [x] Adapt coordinate inputs into an ergonomic 2-row layout at $\le 540\text{px}$ to prevent text clipping of 17-character DMS strings.
- [x] Eliminate unprompted iOS Safari viewport auto-zoom on focus by enforcing `font-size: 16px !important;` on mobile inputs and selects.
- [x] Harmonize Usage / Telemetry dashboard on mobile viewports (12px padding, stacked Active Period deck, touch-scrubbing SVG timeline chart).
- [x] Render floating dialogs (Help, Feedback, Grid Info) as full-width bottom modal sheets on mobile viewports ($\le 768\text{px}$).
- [x] Expand Playwright test suite with dedicated mobile specs (`webapp/e2e/mobile.spec.ts`) verifying all 11 mobile test scenarios.
- [x] Verify full regression suites: Vitest unit tests (83/83), Vite production build, Playwright full suite (70/70), and Pytest geodetic benchmarks (54/54).

---

## Implementation Details

### 1. Desktop Dockable Workbench (`webapp/src/main.ts`, `webapp/src/style/layout.css`)
- Added vertical divider collapse toggle button (`#workbench-dock-toggle`, `.workbench-dock-toggle`) pinned between `.utility-panel` and `.map-panel`.
- Clicking the toggle smoothly animates `.utility-panel` shut to the left (`margin-left: -480px; transition: margin-left var(--transition-normal);`).
- Invokes `mapView.invalidateSize()` both on transition completion and on initial application load when starting in the collapsed state from `localStorage['pytransdat_workbench_collapsed']`.
- Features directional chevron flip (`◀` when expanded, `▶` when collapsed), accessible tooltip, dynamic `aria-label`, and `aria-expanded` attributes.
- Added global keyboard shortcut: `Ctrl+\` or `Alt+W`.

### 2. Mobile Bottom Sheet (Google Maps Style)
- Under `@media (max-width: 768px)`, `.workbench-layout` switches from desktop flex row to a full-screen Leaflet map (`position: absolute; inset: 0;`), with `.utility-panel` positioned as an anchored bottom sheet.
- Added `.bottom-sheet-handle-bar` with a grab pill indicator.
- Supported 3 ergonomic states:
  - `sheet-peek`: ~72px height, keeping 90% of the map visible for location exploring and tracking.
  - `sheet-half`: 45% screen height (default on load), allowing simultaneous map tapping and coordinate viewing.
  - `sheet-full`: 85vh height for editing large text batches, viewing file queues, and configuring settings.
- Added touch gesture engine (`touchstart`, `touchmove`, `touchend`) with velocity/offset tracking, synthetic click suppression, and automatic elevation to `sheet-full` whenever any coordinate input receives focus.

### 3. Coordinate Inputs & Non-Numeric DMS Preservation
- Maintained `<input type="text">` across `PointConverter.ts` to allow typing degree, minute, second symbols (`°`, `'`, `"`), spaces, and cardinal direction prefixes (`N`, `E`).
- At `@media (max-width: 540px)`, `.coord-grid-3` switches from 3 columns to 2 equal columns for Coordinate 1 & 2 (`1fr 1fr`), with Coordinate 3 (Height / Z) spanning full width underneath.
- Enforced `font-size: 16px !important;` on `.form-input` and `.form-select` on mobile ($\le 768\text{px}$) to prevent iOS Safari from zooming into the page on focus.

### 4. Responsive Header Navbar
- At `@media (max-width: 768px)`, hides `.header-repo-link` ("powered by ANCPI/Romgeo") to save horizontal space.
- Condenses main tab labels to compact text (`Transformare` and `Statistici`).
- Retired `#header-mobile-view-btn` cleanly since the Google Maps-style Bottom Sheet natively provides intuitive map/tools access.

### 5. Telemetry / Usage Mobile Harmonization
- Reduced `.telemetry-page-container` padding on mobile viewports to `12px 14px`.
- Wrapped Active Period deck rows and secondary tabs cleanly into vertical flex groups.
- Added `overflow-x: auto;` to `.spatial-matrix-container` for smooth punchcard inspection on small viewports.
- Added touch event support (`touchstart`, `touchmove`, `touchend`) with `touch-action: pan-y;` on the SVG timeline chart, placing tooltips ergonomically above the user's finger with boundary flipping.

### 6. Mobile Bottom Modal Drawers
- Converted floating dialogs (Help, Feedback, Grid Info) to full-width bottom modal sheets on mobile viewports ($\le 768\text{px}$), hiding desktop minimize buttons and disabling mouse drag conflicts.

---

## Verification & Testing

1. **Frontend Vitest Unit Suite**:
   - Command: `npm --prefix webapp test`
   - Result: 83/83 unit tests passing across all 6 test files (`303ms`).
2. **Vite Production Bundle Build**:
   - Command: `npm --prefix webapp run build`
   - Result: Clean TypeScript compilation (`tsc -b`) and minified production build in `443ms` with zero warnings or errors.
3. **Playwright E2E Suite (70/70 Passed)**:
   - Command: `powershell -ExecutionPolicy Bypass -File .\test-e2e.ps1`
   - Result: All 70 tests passed across 12 spec files (`44.3s`).
   - Covered 11 new mobile scenarios in `webapp/e2e/mobile.spec.ts`:
     1. Mobile bottom sheet defaults to `sheet-half` and cycles through `sheet-full`, `sheet-peek`, and `sheet-half` on handle tap.
     2. Map remains interactive on mobile and receives clicks in half and peek modes.
     3. Coordinate inputs adapt to 2-row layout at $\le 540\text{px}$ and preserve non-numeric DMS text (`45°50'28.44420"`).
     4. Responsive header navbar hides repo link and retires mobile view toggle button.
     5. Usage/telemetry dashboard harmonizes for mobile viewports without horizontal clipping.
     6. Floating help and feedback panels render as bottom overlays on mobile.
     7. Desktop dockable workbench collapses and expands smoothly via toggle button and keyboard shortcut.
     8. Desktop dockable workbench starts collapsed on initial load when persisted in `localStorage`.
     9. Bottom sheet touch swipe gestures smoothly transition between states.
     10. Focusing coordinate input on mobile automatically expands bottom sheet to `sheet-full`.
     11. Timeline chart supports touch scrubbing displaying crosshair and tooltip on mobile.
4. **Core Geodetic Suite**:
   - Command: `pytest tests/test_trans_ro.py`
   - Result: 54/54 passed with $< 10^{-9}$ tolerance.
5. **Ruff Code Hygiene**: Clean with zero errors.
