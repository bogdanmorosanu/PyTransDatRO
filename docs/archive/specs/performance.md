# PyTransDatRO Performance Improvements Report

This report outlines potential performance improvements for the `pytransdatro` pure-Python coordinate transformation library. The goal is to optimize execution speed to support hundreds of API requests per day, with some requests containing thousands to hundreds of thousands of coordinates, while strictly adhering to the "Pure Python" project invariants.

The improvements are categorized by the scope and effort of the changes.

---

## 1. Minor Improvements (Mathematical and Syntax Optimizations)

These are small changes consisting of a few lines of code, primarily focused on reducing the overhead of Python's execution inside tight mathematical loops.

*   **Replace `math.pow` with direct multiplication:**
    *   **Description:** In `proj_stereo.py` (e.g., inside the iterative loop of `to_geo`), `math.pow(math.sin(r_lat), 2)` is used. Calling `math.pow` has higher overhead than direct multiplication `x * x`.
    *   **Implementation:** Replace it with `sin_r_lat = math.sin(r_lat)` followed by `sin_r_lat * sin_r_lat`.
    *   **Impact Estimate:** Low to Moderate. It saves a few nanoseconds per point, but inside an iterative `while` loop over 100,000 points, it adds up.

*   **Localizing Global Functions in Loops:**
    *   **Description:** Python has a slight overhead when resolving global module functions (like `math.sin`, `math.cos`, `math.tan`).
    *   **Implementation:** When executing calculations on multiple points, assign math functions to local variables before the loop (e.g., `_sin = math.sin`).
    *   **Impact Estimate:** Low. Provides a marginal (~5-10%) speedup for pure mathematical blocks.

*   **Inline Degree/Radian Conversions:**
    *   **Description:** The 1D height grid interpolations use `math.degrees(r_n)` and `math.radians`.
    *   **Implementation:** Inline these as `r_n * (180.0 / math.pi)` or define them as direct lambda/constants if function call overhead is high in bulk processing.
    *   **Impact Estimate:** Low.

*   **Optimize Subgrid Index Generation:**
    *   **Description:** In `trans_grid.py`, `_sgrid_idxs` initializes a list `r = [None] * 16` and uses a double `for` loop to populate it before returning `tuple(r)`.
    *   **Implementation:** Use a tuple comprehension or a flattened unrolled expression to construct the tuple directly in memory, avoiding the overhead of list mutations and Python `for` loops.
    *   **Impact Estimate:** Low to Moderate. This runs for every single point transformed.

---

## 2. Medium Improvements (Function and Concurrency Optimizations)

These changes involve modifying entire functions or adding parallel processing capabilities.

*   **Multiprocessing / Concurrent Processing:**
    *   **Description:** Coordinate transformation is a heavily CPU-bound task. Processing 100k points in a single thread leaves other server cores idle.
    *   **Implementation:** For requests with high point counts (e.g., >10,000), use Python's built-in `concurrent.futures.ProcessPoolExecutor` or `multiprocessing.Pool` to chunk the input arrays and process them in parallel.
    *   **Impact Estimate:** High. Will provide a near-linear speedup proportional to the number of available CPU cores on the server (e.g., ~4x faster on a 4-core machine).

*   **Using `__slots__` or NamedTuples for Points:**
    *   **Description:** If coordinates are ever wrapped in classes or objects as they pass through the API, standard Python classes use dictionaries for attribute storage.
    *   **Implementation:** Use `__slots__ = ['n', 'e', 'z']` in coordinate classes, or use `collections.namedtuple`.
    *   **Impact Estimate:** Moderate. Drastically reduces memory footprint and slightly improves attribute access times when holding hundreds of thousands of objects in memory.

*   **Optimize iterative loop in `StereoProj.to_geo`:**
    *   **Description:** The latitude approximation iterative `while` loop has a tolerance condition.
    *   **Implementation:** Evaluate if `math.isclose()` or a fixed number of iterations (e.g., unrolling the loop for a guaranteed 3-4 iterations based on known numerical stability) can replace the `while` loop condition overhead.
    *   **Impact Estimate:** Moderate.

---

## 3. Major Improvements (Architecture and Structural Changes)

These changes require refactoring how the data flows through the application or changing the execution environment.

*   **Implement Bulk Processing Methods (`_many` suffix):**
    *   **Description:** Currently, `st70_to_etrs89` processes a single point. If 100k points arrive via API, calling this function 100k times from a `for` loop introduces massive Python function call overhead.
    *   **Implementation:** Create `st70_to_etrs89_many(n_array, e_array, z_array=None)` that accepts standard library `array.array` or plain lists. The method internally loops over the arrays, meaning the function is called once, and only native structures and built-in math operate in the loop. 
    *   **Impact Estimate:** High. Can reduce total execution time by 20-30% simply by avoiding function call overhead (Python `def` invocations are relatively slow).

*   **Spatial Grouping for Grid Interpolation Caching:**
    *   **Description:** In `Grid2D`, bicubic interpolation requires fetching 16 coefficients for a 4x4 subgrid (`BiInterp` instantiation). `functools.lru_cache` helps, but we can do better for bulk data.
    *   **Implementation:** When processing bulk points, first group them by their grid cell ID. For all points falling inside the same grid cell, fetch the 16 grid values and instantiate `BiInterp` exactly **once**. Then, evaluate the polynomial for all points in that cell.
    *   **Impact Estimate:** Very High. Will drastically reduce memory lookups, tuple hashing (from `lru_cache`), and class instantiation overhead.

*   **Deploy using PyPy (JIT Compiler):**
    *   **Description:** Given the strict invariant of "Pure Python Implementation" with no C extensions (like Cython or NumPy), the math-heavy execution remains bottlenecked by the standard CPython interpreter.
    *   **Implementation:** Run the API webserver and `pytransdatro` package using **PyPy** instead of CPython. PyPy's Just-In-Time (JIT) compiler optimizes pure-Python `for` loops and mathematical operations at runtime.
    *   **Impact Estimate:** Extreme. PyPy typically executes tight numerical mathematical loops **5x to 10x faster** than CPython without changing a single line of code. This is the most significant performance gain possible while strictly respecting the project's invariants.
