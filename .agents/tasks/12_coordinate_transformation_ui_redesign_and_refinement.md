# Coordinate Transformation Configuration Redesign & Refinement

## Objectives
- [x] Streamline left workbench panel layout (`480px`) and folder bookmark tabs with seamless connection to active cards.
- [x] Implement theme-native high-contrast primary CTA button (`.btn-primary`) without generic saturated accent colors.
- [x] Create 3-tier Shared Configuration Bar with sliding direction toggle, balanced 50/50 angular unit and coordinate order selectors, and dynamic copyable angle format example badge.
- [x] Symmetrical 3-column coordinate layout with pure numeric outputs and persistent structure preventing vertical layout shift.
- [x] Promote Quasigeoid Height ($\zeta = h - H$) to a dedicated metric field on its own row styled in technical cyan-blue.
- [x] Integrate unified Quality Report cards across Single Point, Batch Text, and Bulk Files tabs.
- [x] Decouple Leaflet map auto-panning from coordinate transformation; implement decoupled, boundary-safe `Center on Map` buttons guarded strictly within Romania (`[43.5, 20.2]` to `[48.3, 29.8]`).
- [x] Separate Leaflet map layers into three distinct LayerGroups (`pointLayer`, `textLayer`, `filesLayer`) for independent visibility and rendering.
- [x] Implement smart map click capturing coordinates directly or reverse-calculating Stereo 70 projection meters based on active transformation direction.
- [x] Implement strict numeric validation (`Number(token)` / strict regex) preventing non-numeric typos (e.g. `5001a50.000`) from being parsed into false coordinates.
- [x] Enforce strict 1-to-1 row correspondence in Batch Text conversion with human 1-based indexing for unnamed coordinates.
- [x] Streamline Bulk Files processing: remove raw text previews, implement interactive queue management with per-file removal, and provide post-transformation cumulative quality reports.
- [x] Synchronize bilingual localization (`en` and `ro`) with standardized geodetic terminology.

## Implementation Details

1. **Design Tokens & Theme-Native Styling (`webapp/src/style/variables.css`, `webapp/src/style/components.css`)**:
   - Replaced generic cobalt blue with high-contrast anchor `.btn-primary` tokens:
     - Dark theme: Porcelain surface (`#f8fafc`) with deep charcoal text (`#090d16`).
     - Light theme: Slate surface (`#0f172a`) with crisp white text (`#f8fafc`).
   - Added geodetic metric color tokens: `--color-metric-text` (`#38bdf8` / `#0284c7`), `--color-metric-border`, and `--color-metric-bg`.
   - Eliminated the 16px background gap between `.tabs-bar` and active panel cards, creating physical folder bookmark tabs.
   - Styled `.metric-field`, `.quality-report-card`, `.format-example-badge`, and `.action-footer-right`.

2. **Application Orchestration & Shared Configuration (`webapp/src/main.ts`)**:
   - Lifted transformation direction (`ETRS89ToStereo70` default / first button in toggle) and coordinate order (`NE(H)` vs `EN(H)`) to the shared configuration header.
   - Configured dynamic angle format badge with click-to-copy toast, copying strictly the formatted numeric angle string without label text.
   - Wired smart map click handler:
     - `ETRS89ToStereo70`: formats clicked latitude and longitude into active angular units.
     - `Stereo70ToETRS89`: calls backend `apiClient.transformPoint` to reverse-calculate the clicked point into Stereo 70 meters.
     - Out-of-grid clicks display warning toasts and do not populate invalid values.

3. **Single Point Transformation Workbench (`webapp/src/components/PointConverter.ts`)**:
   - Initialized default coordinates to Romania's projection center origin ($\phi = 46.0^\circ\text{ N}, \lambda = 25.0^\circ\text{ E}, h = 100.0\text{ m}$).
   - Micro-headers: `Input coordinates (or click on map)` and `Computed coordinates`.
   - Added top-right code-block style copy icon button above computed coordinates with temporary `✓ Copied` visual feedback.
   - Configured focus auto-selection (`input.select()`) across all coordinate inputs.
   - Output boxes contain pure numeric values (stripped of `m`, `°`, `rad`), persistently displaying `—` when cleared or invalid to prevent vertical layout shifts.
   - Dedicated Quasigeoid Height box rendered on its own row below coordinates.
   - Unified Quality Report card placed below outputs summarizing valid coordinates and grid envelope status.
   - Decoupled `Center on Map` button placed in the bottom-right action footer.

4. **Batch Text Transformation Precision (`webapp/src/components/BatchConverter.ts`, `webapp/src/utils/pointParser.ts`)**:
   - Replaced permissive `parseFloat` parsing with `isValidDecimal` strict regex checking (`/^[+-]?(?:\d+\.?\d*|\.\d+)$/`) and `Number(token)`.
   - Guaranteed strict 1-to-1 line correspondence: input row $N$ maps directly to output row $N$. Empty lines, comments, and invalid rows maintain their original line index.
   - Human 1-based indexing (`Point 1`, `Point 2`, ...) applied for unnamed coordinate points.
   - Added top-right `Copy Results` button to computed results card header.
   - Integrated input pre-flight diagnostic badge and output Quality Report card.

5. **Bulk Files Queue & Quality Reporting (`webapp/src/components/FileUpload.ts`)**:
   - Removed raw 5-line text preview in favor of clean queue metadata.
   - Added individual `✕` remove buttons per queued file item.
   - Replaced automatic processing with explicit `Transform Files` action.
   - Rendered per-file and cumulative Quality Report cards catching valid points, `# Outside grid`, `# Wrong dimension`, and `# Invalid coordinate data`.
   - Integrated `Download All (ZIP)` and individual download links with decoupled bottom-right `Center on Map` action.

6. **Web GIS Map Layer Separation & Dual Navigation Actions (`webapp/src/components/MapView.ts`, `webapp/src/utils/gridBounds.ts`)**:
   - Created three independent Leaflet LayerGroups (`pointLayer`, `textLayer`, `filesLayer`) added to Leaflet layer controls.
   - Decoupled coordinate transformation from automatic map zooming/panning.
   - Implemented dual map navigation controls across all tabs:
     - `Center on Map`: Pans to centroid of valid in-grid coordinates while strictly preserving the current map zoom.
     - `Zoom to Extents`: Pans and calls `fitBounds` to encompass all plotted points, or zooms to level 14 for a single point.
   - Bounded by coarse SPG envelope bounds (`isWithinStereo70Grid`, `isWithinETRS89Grid`, `isInsideRomaniaGrid`), preventing erratic viewport jumps.

7. **File Streaming Robustness & Test Datasets (`api/routes/file.py`, `api/schemas/point.py`, `tests/`)**:
   - Enhanced backend file transformation pipeline:
     - Auto-detects Point ID column (distinguishes 3 or 4 columns with optional ID).
     - Parses DMS coordinates using `geoutils.sexa_dms_chars_to_rad`.
     - Supports `unit="dms"` returning clean unspaced DMS strings (`DD°MM'SS.SSSSS"`).
     - Decodes legacy files with Windows-1252 / Latin-1 degree symbol fallback.
   - Extracted 200-point benchmark datasets (`tests/data/st70_to_etrs89_input_200.csv`, `tests/data/etrs89_to_st70_input_200.csv`) with full UTF-8 formatting.
   - Added automated tests in `tests/test_api_modern.py`.

8. **Localization Alignment & Quality Warnings (`webapp/src/i18n/en.ts`, `webapp/src/i18n/ro.ts`)**:
   - Synchronized all new and modified strings across English and Romanian locales.
   - Standardized concise Quality Report warnings: `In Grid`, `Outside Grid`, `No Grid Data`, and `Syntax Error` (RO: `În grid`, `În afara gridului`, `Fără date în grid`, `Eroare sintaxă`).
   - Standardized geodetic terminology: `Quasigeoid Height (ζ = h - H)` (RO: `Înălțimea cvasigeoidului (ζ = h - H)`).

## Verification & Testing

- **Frontend Production Build**:
  - Command: `npm --prefix webapp run build`
  - Output: `tsc && vite build` completed in `350ms` with `0` errors.
- **Backend Test Suite**:
  - Command: `D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/`
  - Result: `73/73 passed` in 5.09s (including new 200-point streaming tests, legacy API, modern API, rate limiting, telemetry, and pure Python core).
- **Runtime Services**:
  - FastAPI Backend: `http://127.0.0.1:8000` (Status: `200 OK`).
  - Vite Dev Server: `http://localhost:5173` (Status: `200 OK`).
