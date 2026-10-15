# Code Quality & Architectural Standards

These standards govern code modifications across all layers of the PyTransDatRO project. All AI agents and developers working on this repository must adhere to these invariants.

---

## 1. General Repository Invariants

- **Zero Regression Rule**: No code modification may alter mathematical transformation results by more than numerical machine precision ($10^{-9}$ units). Always verify against ANCPI comparison reports.
- **Strict Layer Decoupling**:
  - `pytransdatro/` (Core): Standard library Python only (+ `numpy` strictly for deserializing binary `.spg` grids). Must NEVER import from `api/`, `fastapi`, or web dependencies.
  - `api/`: Built using FastAPI and Pydantic. Consumes `pytransdatro` as a core library. Must NEVER contain client-side business logic.
  - `webapp/`: Static SPA client using TypeScript. Consumes `api/` via HTTP endpoints. Must remain cleanly decoupled from backend internals.
- **Fail-Safe Exception Handling**:
  - Never use bare `except:` or `except Exception: pass`.
  - Always catch specific exceptions (e.g. `(ValueError, KeyError, FileNotFoundError)`).
  - Log suppressed or handled exceptions at `logger.debug` or `logger.warning` with contextual diagnostic details.

---

## 2. Python Core & API Standards

### Type Annotations
- All functions, methods, and class attributes must include complete Python type annotations (e.g. `tuple[float, float]`, `Optional[float]`, `dict[str, Any]`).
- Internal helper methods must be cleanly typed and prefixed with a single underscore `_`.

### Configuration & Constants
- Never hardcode arbitrary operational numbers (ports, timeouts, batch caps, file paths) in business logic.
- Operational settings reside in `api/config.py` backed by `pydantic-settings` and `.env`.
- Mathematical and geodetic invariants (ellipsoid axes, projection origins, false coordinates) reside strictly in `pytransdatro/constants.py`.

### Docstrings
- Follow reST / Sphinx field list formatting:
  - Document `:param <name>:` and `:type <name>:`.
  - Document `:return:` and `:rtype:`.
  - Document `:raises <Error>:` for all explicitly raised custom exceptions.
- Provide clear mathematical formulas in algorithmic modules (`proj_stereo.py`, `trans_grid.py`).

### Code Formatting
- Code is formatted using `ruff` with a maximum line length of 100 characters.
- Imports must be sorted using standard PEP 8 grouping: standard library, third-party libraries, local package imports.

---

## 3. Frontend (Web Application) Standards

### Component & Styling Architecture
- Never style individual elements ad-hoc with inline styles or arbitrary CSS classes.
- All styles must derive from CSS custom properties in `variables.css`.
- Use modular CSS stylesheets divided by role (`forms.css`, `buttons.css`, `badges.css`, `layout.css`, `telemetry.css`).
- Use reusable vanilla TS UI factory functions (`webapp/src/ui/`) for common controls (segmented toggles, select pickers, copy badges).

### Coordinate Parsers & Exporters
- Any modification to `dmsParser.ts`, `pointParser.ts`, or `exporters.ts` must be covered by unit tests in `vitest`.
- Parsers must be resilient to empty lines, comments, and whitespace variations without throwing unhandled exceptions.

### Accessibility (a11y)
- All custom toggles and buttons must support keyboard navigation (`Tab`, `Space`, `Enter`).
- Custom toggles must include `role="switch"` or `role="button"` with `aria-checked` or `aria-label` attributes.
- High-contrast focus outlines (`var(--focus-border)`) must be maintained on all interactive inputs.
