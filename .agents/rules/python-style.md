# Python Code Style & Conventions

## Naming Conventions
- **Modules**: `snake_case` (e.g., `trans_ro.py`, `trans_grid.py`, `proj_stereo.py`).
- **Classes**: `PascalCase` (e.g., `TransRO`, `Grid2D`, `StereoProj`). Internal helper classes use a single leading underscore (e.g., `_Ellipsoid`, `_ConfSph`). Exceptions use `PascalCase` ending in `Err` (e.g., `OutOfGridErr`, `NoDataGridErr`).
- **Functions & Methods**: `snake_case` (e.g., `st70_to_etrs89`, `coo_at_idx`, `sexa_to_rad`). Protected/internal methods use a leading underscore (e.g., `_grid1d_sel`, `_sgrid_idxs`, `_init_interp`).
- **Variables**: Concise, domain-specific `snake_case`. Coordinate variables standardly use `n` (Northing/latitude axis), `e` (Easting/longitude axis), and `z` / `h` (elevation/height).
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `DEG_TO_RAD_FACTOR`, `WGS84_ELL_SMAJ_AXIS`, `FALSE_N`).

## Docstrings & Documentation
- **Format**: Follow Sphinx / reStructuredText (reST) field list conventions:
  - Document parameters with `:param <name>:` and `:type <name>:`.
  - Document return values with `:return:` and `:rtype:`.
  - Document exceptions with `:raises <Error>:`.
  - Use Sphinx directives when applicable (e.g., `.. warning::`).
- **Module Docstrings**: Every module starts with a triple-quoted docstring outlining its purpose, function/class listings, notes on coordinate units, and basic usage examples.

## Commenting Patterns
- **Inline Annotations**: Provide brief inline comments clarifying mathematical formulas, units, or algorithm parameters (e.g., `# DMS value`, `# translation on north`, `# scale (ppm)`).
- **Block Comments**: Use block or section comments to demarcate major algorithmic steps or conditional branches (e.g., `# IF z value provided, apply 1D grid transformation (height correction)`).
- **Purity**: Rely exclusively on standard library modules (e.g., `math`, `struct`, `functools`, `abc`, `os`).
