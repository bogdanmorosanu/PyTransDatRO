# PyTransDatRO Web Application (`webapp/`)

Modern Web GIS Single Page Application (SPA) serving as the browser interface for **PyTransDatRO** (Romania's official geodetic coordinate transformation platform).

Built with **Vite**, **TypeScript**, **Vanilla CSS**, and **Leaflet**.

---

## 1. Quick Start

### Prerequisites
- **Node.js**: v18+ (LTS recommended)
- **FastAPI Backend**: Running at `http://127.0.0.1:8000` (see repo root instructions)

### Development Workflow

```bash
# Navigate to the frontend directory
cd webapp

# Install dependencies
npm install

# Start development server (with hot module replacement)
npm run dev
```

The application will be accessible at:
👉 **`http://localhost:5173/`**

> **Note on API Proxying**: Vite is configured (`vite.config.ts`) to transparently proxy all `/api`, `/health`, and `/transdatonline` calls to `http://127.0.0.1:8000`. Keep the Python backend running during frontend development.

### Production Build & Integration

```bash
# Type-check TypeScript and build production bundle
npm run build
```

This compiles optimized, minified assets into `webapp/dist/`.

When the FastAPI server starts (`python -m uvicorn api.main:app`), it automatically detects `webapp/dist/` and mounts it on the root URL (`http://127.0.0.1:8000/`) using FastAPI's `StaticFiles(html=True)`.

---

## 2. Directory Structure

```
webapp/
├── index.html              # HTML entrypoint with preloaded fonts and layout shells
├── package.json            # Node dependencies and build scripts
├── tsconfig.json           # Strict TypeScript configuration
├── vite.config.ts          # Vite build options and dev server API proxy
└── src/
    ├── main.ts             # Application orchestrator & component lifecycle manager
    ├── components/         # Interactive UI workbench components
    │   ├── Header.ts       # Branding, tab navigation, language & theme switcher
    │   ├── PointConverter.ts  # Single-point transformation workbench
    │   ├── BatchConverter.ts  # Multi-line text conversion & table preview
    │   ├── FileUpload.ts      # Drag-and-drop file processor with streaming progress
    │   ├── MapView.ts         # Leaflet map, basemap selector, clustering & grid extent
    │   ├── StatsView.ts       # Operational telemetry and active grid metadata
    │   ├── VolumeWarningModal.ts # Guardrail prompt for transforms > 2,000 points
    │   ├── ReportModal.ts     # Diagnostic calculation details modal
    │   └── FeedbackModal.ts   # Discrepancy reporting and user feedback modal
    ├── i18n/               # Internationalization system
    │   ├── index.ts        # Reactive translation store & persistent state
    │   ├── ro.ts           # Romanian translation dictionary (default)
    │   └── en.ts           # English translation dictionary
    ├── style/              # Design system & styles (Pure Vanilla CSS)
    │   ├── variables.css   # Color palette tokens, dark/light themes, typography
    │   ├── base.css        # CSS reset, responsive desktop split-screen grid
    │   ├── components.css  # Buttons, pills, form inputs, tables, modals, toasts
    │   └── map.css         # Leaflet container, custom popups & cluster styling
    └── utils/              # Client-side geodetic & formatting helpers
        ├── apiClient.ts    # Type-safe REST client (passes X-Client-Source: web_app)
        ├── dmsParser.ts    # High-precision DMS <-> Decimal Degrees converter
        ├── pointParser.ts  # Auto-delimiter detector & coordinate cleaner
        ├── exporters.ts    # Client-side CSV, GeoJSON, and AutoCAD DXF generator
        └── toast.ts        # Non-intrusive floating toast notifications
```

---

## 3. Key Architecture & Design Invariants

1. **Decoupled REST API Integration**:
   - The frontend communicates with the transformation engine exclusively over HTTP via `src/utils/apiClient.ts`.
   - All requests send the HTTP header `X-Client-Source: web_app`, allowing the backend `SqliteUsageLogger` to attribute telemetry appropriately.

2. **Vanilla CSS Design System**:
   - Built with pure CSS custom properties (`src/style/variables.css`).
   - Supports **Dark Mode**, **Light Mode**, and **System Preference** with no heavy UI framework dependencies.
   - Clean, modern geodetic aesthetic: Slate neutrals, Cyan and Emerald precision accents, monospace coordinate tables (`JetBrains Mono`).

3. **Bilingual Internationalization (`i18n`)**:
   - Romanian (`ro`) is the primary default language; English (`en`) is fully supported.
   - User language preference is stored in browser `localStorage`.
   - UI re-renders reactively upon language change without full page reloads.

4. **Map Performance & Guardrails**:
   - Built with **Leaflet** and **Leaflet.markercluster**.
   - **Large volume guardrail**: When batch transformations exceed **2,000 points**, the user is prompted via `VolumeWarningModal` to choose between full rendering, a 1,000-point sample preview, or skipping map rendering to preserve browser responsiveness.
   - Base maps: OpenStreetMap, CartoDB Positron, CartoDB Dark Matter, and ESRI World Imagery.
   - Raster shift overlays: Placeholders prepared for forthcoming $\Delta E$, $\Delta N$, and quasigeoid $\zeta$ map services.

5. **Client-Side Exporting**:
   - Results can be downloaded directly in the browser as:
     - **CSV** (Comma, Semicolon, or Tab delimited)
     - **GeoJSON** (`Point` features with coordinate and height properties)
     - **AutoCAD DXF** (R12 ASCII format with `POINT` entities for immediate CAD import)

---

## 4. Verification & Testing

To verify code quality and build integrity:

```bash
# Type check and verify production bundle compiles cleanly
npm run build
```
