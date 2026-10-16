# Post-Refactor AI Agent Continuous Development Readiness

## Objectives
- [x] Transition AI agent system instructions from the completed phased refactor roadmap to a permanent continuous development lifecycle.
- [x] Seal `refactor_roadmap.md` with a permanent completion notice while preserving it as an architectural baseline.
- [x] Diagnose and document the Windows sandbox testing hang issue caused by Starlette `TestClient` / `anyio` loopback socket blocking.
- [x] Codify permanent refactoring quality standards across `GEMINI.md`, `.agents/rules/code-quality-standards.md`, and `.agents/skills/`.
- [x] Integrate Vitest frontend test execution (`npm --prefix webapp test`) across code auditing and regression verification workflows.

---

## Implementation Details

1. **Agent Workflow & System Instructions (`GEMINI.md`)**:
   - Replaced "Refactoring & Code Quality Mode" (which bound agent focus to the phased roadmap) with "Continuous Development & Maintenance Mode".
   - Declared `refactor_roadmap.md` as 100% completed, establishing it as the immutable baseline for clean architecture and code standards.
   - Codified permanent refactoring invariants for all subsequent work: pure Python core library purity, strict 3-tier layer isolation (`pytransdatro/` $\to$ `api/` $\to$ `webapp/`), zero functional regression ($< 10^{-9}$ tolerance), complete typing, Sphinx/reST docstrings with mathematical formulas, centralized configuration, and modular CSS/UI factories.
   - Codified Windows Sandbox Test Execution Rules: documented that Starlette's `TestClient` uses `anyio` to bridge synchronous test calls to asynchronous FastAPI endpoints via worker threads and local loopback sockets; under Windows sandbox isolation (`BypassSandbox: false`), local loopback socket communication is blocked, causing full `pytest` or API tests to hang indefinitely. Documented the requirement to run `pytest` and API tests with `BypassSandbox: true`.

2. **Roadmap Sealing (`refactor_roadmap.md`)**:
   - Added a top-level alert banner declaring status as `COMPLETED (Phases 1–6 Fully Executed & Sealed)`.
   - Directed all active and future development to `GEMINI.md` and `.agents/rules/`.

3. **Code Quality Standards (`.agents/rules/code-quality-standards.md`)**:
   - Integrated frontend unit testing expectations using Vitest (`npm --prefix webapp test`) for all parsers (`dmsParser.ts`, `pointParser.ts`), boundary envelopes (`gridBounds.ts`), exporters (`exporters.ts`), and UI factory components (`webapp/src/ui/`).
   - Reinforced production build verification (`npm --prefix webapp run build`).

4. **Agent Tooling Skills Updates (`.agents/skills/`)**:
   - **`code-audit`**: Removed obsolete Phase 1 note regarding Ruff installation. Added Vitest unit test execution (`npm --prefix webapp test`) to verify frontend logic alongside `npm --prefix webapp run build`.
   - **`regression-verify`**: Added Windows sandbox execution notice explaining why Step 2 (and full `pytest`) must use `BypassSandbox: true`. Added Step 4: Web Application Utility & UI Regression Suite (`npm --prefix webapp test`, 65 tests).
   - **`webapp-record`**: Standardized command execution using `npm --prefix webapp test` and `npm --prefix webapp run build` without directory changes.

---

## Verification & Testing
- **Core Geodetic Unit Suite (Sandboxed)**:
  `D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_trans_ro.py`
  Result: 46/46 passed in 2.45s.
- **Frontend Vitest Suite (Sandboxed)**:
  `npm --prefix webapp test`
  Result: 65/65 passed across 5 test suites in 322ms.
- **Frontend Production Bundle (Sandboxed)**:
  `npm --prefix webapp run build`
  Result: Compiled cleanly in 373ms with zero TypeScript or bundler errors.
- **Full Python Test Suite (BypassSandbox: true)**:
  `D:\anaconda3\envs\pytransdat\python.exe -m pytest`
  Result: 98/98 passed in 6.22s.
- **Ruff Linter & Formatter (Sandboxed)**:
  `ruff check` and `ruff format --check` passed clean across all 50 files.
- **ANCPI Reference Benchmark Comparison (Sandboxed)**:
  `tools/generate_ancpi_comparison_reports.py` executed with 100% byte-for-byte fidelity across 9,046 points ($< 10^{-9}$ tolerance).
