---
name: code-audit
description: Audits code quality across pytransdatro, api, and webapp layers using Ruff, TypeScript compiler, and pattern checks. Use when reviewing code quality, checking for regressions or code smells, or before finishing a refactoring chat.
---

# Code Audit Skill

## Purpose
This skill executes a complete, cross-layer static analysis and code hygiene audit across the PyTransDatRO repository. It verifies that Python formatting, type conventions, layer isolation boundaries, and TypeScript builds comply with project invariants.

---

## Audit Workflow

### 1. Python Code Hygiene & Ruff Lint
Run the Ruff linter and formatter check from the designated environment:

```powershell
D:\anaconda3\envs\pytransdat\python.exe -m ruff check pytransdatro api tests tools
D:\anaconda3\envs\pytransdat\python.exe -m ruff format --check pytransdatro api tests tools
```

*Note: If Ruff is not yet installed in Phase 1, verify syntax and basic compilation using `python -m py_compile`.*

### 2. Forbidden Anti-Pattern Scan
Check for prohibited anti-patterns across the codebase:

1. **Bare `except: pass` or `except Exception: pass`**:
   Ensure no silent error suppression exists:
   ```powershell
   Get-ChildItem -Path pytransdatro, api -Recurse -Filter "*.py" | Select-String -Pattern "except.*:\s*pass"
   ```

2. **Core Layer Boundary Violation**:
   Ensure `pytransdatro/` has ZERO imports from `api/`, `fastapi`, `pydantic`, `uvicorn`, `httpx`, or other web frameworks:
   ```powershell
   Get-ChildItem -Path pytransdatro -Recurse -Filter "*.py" | Select-String -Pattern "from api|import api|fastapi|uvicorn|pydantic"
   ```

3. **Hardcoded Localhost / Magic Ports**:
   Ensure application logic does not contain hardcoded URLs/ports that bypass configuration:
   ```powershell
   Get-ChildItem -Path api -Recurse -Filter "*.py" | Select-String -Pattern "localhost:8000|127.0.0.1:8000"
   ```

### 3. Frontend TypeScript & Bundler Check
Verify that the web application compiles cleanly without type errors:

```powershell
npm --prefix webapp run build
```

---

## Output Report
Summarize findings categorized by:
- **Clean**: Checks that passed with zero issues.
- **Warnings**: Potential code smells, missing type hints, or unorganized imports.
- **Errors / Invariants**: Boundary breaches, silent error swallowing, or TypeScript compilation failures.
