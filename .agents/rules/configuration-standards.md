# Configuration & Environment Standards

This document establishes the configuration hierarchy and rules for externalizing settings in PyTransDatRO.

---

## 1. Principles

1. **Separation of Concerns**: Geodetic physics invariants remain in code constants; operational deployment knobs live in external configuration.
2. **Fail-Fast Startup**: Any missing required configuration or invalid data type must halt startup with an actionable human-readable message.
3. **12-Factor App Compatibility**: Configuration flows from environment variables > `.env` file > code defaults.

---

## 2. Parameter Classification

| Type | Examples | Location | Configurable via `.env`? |
| :--- | :--- | :--- | :--- |
| **Geodetic Invariants** | WGS84 $a$, $f$, Stereo 70 origin $\phi_0, \lambda_0$, false northing/easting | `pytransdatro/constants.py` | **NO** |
| **Grid Parameters** | Helmert 2D translations, rotation, scale, grid dimensions | `.spg` file binary header | **NO** (Dynamic from SPG) |
| **Server Operations** | Host, port, log level, DB path, CORS origins | `api/config.py` (Pydantic) | **YES** |
| **Safety & Throttling** | Rate limits, max batch rows, max file upload size | `api/config.py` (Pydantic) | **YES** |
| **Web Client Settings** | Default map center, zoom, client batch warning cap | `webapp/src/config.ts` | **YES** (`VITE_*`) |

---

## 3. Implementation Rules

### Backend (`api/config.py`)
- Built with `pydantic-settings` (`BaseSettings`).
- Provides sensible local development defaults so the application can run out-of-the-box without requiring a `.env` file.
- Reads environment variables with a clean prefix where appropriate (e.g. `PYTRANSDAT_PORT`, `PYTRANSDAT_LOG_LEVEL`).
- Must be registered in FastAPI via dependency injection (`get_settings`).

### Frontend (`webapp/src/config.ts`)
- Accesses variables via `import.meta.env.VITE_*`.
- Provides explicit fallbacks for local development:
  ```typescript
  export const APP_CONFIG = {
    apiBaseUrl: import.meta.env.VITE_API_BASE_URL || '',
    defaultMapCenter: [45.9432, 24.9668] as [number, number],
    defaultMapZoom: 7,
    batchWarningThreshold: 5000,
  };
  ```

### Documentation (`.env.example`)
- Every configurable variable must be listed in `.env.example` at repo root with:
  - An explanatory comment.
  - The default value.
- Real `.env` files and local database files (`*.db`, `*.sqlite`) must remain in `.gitignore`.
