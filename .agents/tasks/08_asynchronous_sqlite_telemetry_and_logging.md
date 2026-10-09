# Asynchronous SQLite Usage Telemetry & Logging

## Objectives
- [x] Design an optional, high-performance telemetry and logging system to track coordinate transformation activity over time (operational usage trends and spatial density heatmaps).
- [x] Maintain zero external dependencies, adhering strictly to the pure Python standard library (`sqlite3`, `queue`, `threading`, `datetime`, `math`, `pathlib`, `atexit`).
- [x] Guarantee zero latency and zero disk I/O overhead on geodetic transformation threads by using an in-memory queue and background daemon worker thread.
- [x] Ensure strict failure isolation: geodetic math never fails or throws exceptions due to telemetry errors (disk full, lock contention, permission errors).
- [x] Protect cadastral privacy and minimize storage footprint by anonymizing coordinates into ~1 km geographic cells with integer degree/grad splitting.
- [x] Implement in-place aggregation (atomic UPSERTs) in SQLite using Write-Ahead Logging (WAL) and `NORMAL` synchronous mode for multi-worker concurrency.
- [x] Track client sources (`unknown`, `web_map`, `rest_api`, `desktop_cli`) and daily call/point volume metrics.
- [x] Integrate `SqliteUsageLogger` into `TransRO` via a pluggable observer parameter (`telemetry_logger`) and `source` ID parameter.
- [x] Create comprehensive unit tests verifying queue batching, binning, atomic UPSERTs, daily metrics, and failure isolation.
- [x] Provide diagnostic and performance tooling (`tools/inspect_telemetry_db.py`, `tools/run_perf_telemetry_test.py`).
- [x] Document the architecture and operational usage in `docs/logging.md` and design specifications in `scratch/telemetry_specification.md`.

## Implementation Details

1. **Architecture & Invariants**:
   - **Zero Overhead by Default**: When `TransRO` is initialized without a logger, logging is completely inactive, performing no disk I/O and zero allocations.
   - **Non-Blocking Queue Push**: Calling threads push lightweight tuple pointers into an in-memory `queue.Queue(maxsize=10000)` via `put_nowait()`, completing in $< 0.001\text{ ms}$. If the queue fills under extreme I/O stall, payloads are dropped silently to preserve memory stability.
   - **Failure Isolation Guarantee**: Calls to the logger inside `TransRO` are wrapped in `try/except Exception: pass`, ensuring transformation calculations never fail due to logging errors.
   - **Cadastral Privacy & Anonymization**: Coordinates are converted from radians to degrees, rounded to 2 decimal places ($\sim 1.1\text{ km}$ latitude, $\sim 0.75\text{ km}$ longitude), and aggregated without persisting exact client coordinates.

2. **Storage Architecture & Database Schema (`logs/usage.db`)**:
   - Database connection opened with `timeout=5.0`, `PRAGMA journal_mode = WAL`, and `PRAGMA synchronous = NORMAL` for robust multi-process concurrency on web servers.
   - **`sources` Dimension Table**: Prepopulated with client IDs: `0: unknown`, `1: web_map`, `2: rest_api`, `3: desktop_cli`.
   - **`daily_stats` Volume Table**: Tracks total transformation calls (`count_calls`) and processed points (`count_points`) per calendar day and client source using atomic UPSERTs (`ON CONFLICT(date, source) DO UPDATE ...`).
   - **`usage` Spatial Heatmap Table**: Stores monthly ~1 km cell counts keyed by `(year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source)`.
   - **Degree & Grad Coordinate Splitting**: Coordinates are split into integer degrees (`lat_deg`, `lon_deg`) and hundredths (`lat_grad`, `lon_grad`), fitting into compact 1-byte SQLite variable-length integers instead of 8-byte floating point numbers.
   - **Storage Efficiency**: Entire country-wide monthly usage across all active cells and directions requires only $\approx 3.2\text{ MB}$ due to in-place incrementing (`count = count + excluded.count`).

3. **Asynchronous Daemon Engine (`pytransdatro/telemetry.py`)**:
   - Implemented `SqliteUsageLogger`:
     - Initializes bounded thread-safe `queue.Queue`.
     - Starts a background daemon thread (`_worker_thread`).
     - Maintains thread affinity by opening and managing the `sqlite3` connection strictly inside the worker loop.
     - Drains batches up to 1,000 items from the queue, aggregating spatial hits and daily call/point counts in in-memory dictionaries before executing bulk `executemany()` transactions.
     - Registers an `atexit` hook and `close()` method that signals stop and joins the worker thread with a 2-second timeout to flush all pending queue items upon process exit.

4. **Pipeline Modernization & Integration (`pytransdatro/trans_ro.py`)**:
   - Added optional `telemetry_logger=None` parameter to `TransRO.__init__()`.
   - Added optional `source=0` parameter to `st70_to_etrs89()` and `etrs89_to_st70()`.
   - Intercepts transformation points polymorphically (scalars and sequences) and forwards them to `self._collector.record()` with transformation direction flag (`is_st70_to_etrs89`), dimension flag (`is_3d`), and source ID.

5. **Tooling & Verification**:
   - **Unit Tests (`tests/test_telemetry.py`)**:
     - `test_telemetry_batching_and_binning`: Validates queue flushing, ~1 km coordinate binning, atomic count increments in `usage`, daily call/point counters in `daily_stats`, and foreign key resolution in `sources`.
     - `test_telemetry_failure_isolation`: Verifies that a broken or throwing logger does not interrupt or fail geodetic transformations.
     - Test suite pass rate: 48/48 passed.
   - **Database Inspection Tool (`tools/inspect_telemetry_db.py`)**:
     - Formats and displays database file size, sources table, daily stats, spatial usage totals, directional breakdown, and top 10 geographical hotspots in clean ASCII tables.
   - **Stress & Benchmark Harness (`tools/run_perf_telemetry_test.py`)**:
     - Exercises scenarios S1 through S4 (up to 25,000 points) in forward and reverse transformations, writing realistic multi-source telemetry data to `./logs/usage.db`.
