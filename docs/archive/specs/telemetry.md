# PyTransDatRO: Usage Telemetry & Logging Specification

## 1. Overview & Objectives
This document specifies the architecture, data model, and operational mechanics for the **optional usage telemetry and logging system** in `PyTransDatRO`.

The primary goal is to assess service usage over time (e.g., generating annual density heatmaps and traffic reports) without degrading performance, risking data corruption, or exposing private cadastral coordinates.

### Dual-Use Scenarios:
1. **Server / Web Service Deployment (Primary Logging Scenario):**
   When deployed as an API backend (e.g., via FastAPI/Flask with multiple worker processes) handling coordinate transformations. The service can log aggregated usage metrics safely across concurrent workers.
2. **Embedded Library / Application Integration (Default Scenario):**
   When imported as a library in desktop applications or local scripts. Logging is **disabled by default**, incurring **zero overhead** and writing no files. If a third-party developer wishes to enable logging, they can opt in via the same standard interface.

---

## 2. Invariants & Constraints
- **Pure Python Standard Library:** Compliant with `GEMINI.md`, all telemetry code relies strictly on Python standard library modules (`sqlite3`, `queue`, `threading`, `datetime`, `math`, `pathlib`, `json`, `atexit`). Zero third-party packages.
- **Separation of Geodetic Math from I/O:** Geodetic classes (`TransRO`, `Grid`, `Helmert2D`, `StereoProj`) remain purely mathematical. Telemetry logic is loosely coupled via a pluggable observer interface.
- **Failure Isolation (Non-Fatal Logging):** Under no circumstances can a logging error (e.g., disk full, file locked, permission error) cause a coordinate transformation to fail or raise an unhandled exception.
- **True Zero Latency Impact:** The user's calculation thread must not perform file I/O, database writes, or coordinate rounding.
- **Privacy & Anonymization:** Raw cadastral coordinates are never persisted. Coordinates are binned into ~1 km geographic cells.

---

## 3. Data Specification & Aggregation

### 3.1 Captured Attributes
| Field | Type | Description |
| :--- | :--- | :--- |
| **`year`** | `INTEGER` | Calendar year (e.g., `2026`). Ensures records are self-contained for multi-year analysis. |
| **`week`** | `INTEGER` | ISO week number of the year (1–53) via `datetime.date.isocalendar()`. |
| **`lat`** | `REAL` | Latitude in ETRS89/WGS84 degrees, rounded to **2 decimal places** (~1.1 km). |
| **`lon`** | `REAL` | Longitude in ETRS89/WGS84 degrees, rounded to **2 decimal places** (~0.75 km). |
| **`direction`** | `TEXT` | Transformation direction (`"st70_to_etrs89"` or `"etrs89_to_st70"`). |
| **`is_3d`** | `INTEGER` | `1` if elevation $Z/h$ was supplied, `0` for 2D transformations. |
| **`source`** | `TEXT` | Client origin (e.g., `"web_api"`, `"desktop_app"`, `"unknown"`). |
| **`count`** | `INTEGER` | Accumulated hit count for this exact combination in the given week. |
| **`metadata`** | `TEXT (JSON)` | Serialized JSON payload for forward-compatible, arbitrary client data. |

### 3.2 Coordinate Precision & Sequence Handling
- **Scalar vs. List Support:** `TransRO` accepts both single numbers and sequences (lists/tuples). The queue payload accepts either. The background worker checks if input coordinates are scalar or iterable:
  - If scalar: binned as a single point count.
  - If sequence: iterates over the points, converting and binning all coordinates into the in-memory batch.
- **Units & Rounding:** Coordinates in `TransRO` are in radians; the background worker converts them to degrees (`math.degrees()`) before binning.
- **2 Decimal Places Resolution:**
  - $0.01^\circ$ Latitude $\approx 1.1\text{ km}$.
  - $0.01^\circ$ Longitude (at Romania's latitude $\sim 45^\circ$) $\approx 0.75\text{ km}$.
- **Elevation Handling:** Elevation ($Z$) is ignored during spatial binning so points from the same area map to the same 2D geographic cell.

---

## 4. Storage Architecture: SQLite Database (`usage_<YYYY>.db`)

### 4.1 Storage Rationale & In-Place Aggregation (UPSERT)
Rather than an append-only log that duplicates identical locations, the system performs **in-place aggregation (UPSERT)**:
- Multiple transformation hits in the same 1 km cell during the same week update an existing record (`count = count + N`).
- Stored as a single binary file per year (e.g., `logs/usage_2026.db`).
- An entire year of active service across thousands of unique cells takes under **1 MB** of disk space.
- Adding `year` directly into the table ensures that merged databases or multi-year analytics queries remain unambiguous.

### 4.2 Database Schema & Atomic UPSERT
```sql
CREATE TABLE IF NOT EXISTS usage (
    year INTEGER NOT NULL,
    week INTEGER NOT NULL,
    lat REAL NOT NULL,
    lon REAL NOT NULL,
    direction TEXT NOT NULL,
    is_3d INTEGER NOT NULL,
    source TEXT NOT NULL,
    count INTEGER NOT NULL DEFAULT 1,
    metadata TEXT,
    PRIMARY KEY (year, week, lat, lon, direction, is_3d, source)
);

CREATE INDEX IF NOT EXISTS idx_usage_year_week ON usage(year, week);
```

**Atomic Upsert Query:**
```sql
INSERT INTO usage (year, week, lat, lon, direction, is_3d, source, count, metadata)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT(year, week, lat, lon, direction, is_3d, source)
DO UPDATE SET count = count + excluded.count;
```

### 4.3 Multi-Worker Server Concurrency & Thread Affinity
1. **Thread Affinity:** In Python, `sqlite3` connections cannot be shared across threads. The SQLite connection **must be opened inside the worker thread's execution loop**, not in `__init__`.
2. **Write-Ahead Logging (WAL):** Enables concurrent readers/writers across multiple server worker processes:
   ```sql
   PRAGMA journal_mode = WAL;
   PRAGMA synchronous = NORMAL;
   ```
3. **Busy Timeout:** SQLite connection initialized with `timeout=5.0` to wait if another worker process is actively writing.

---

## 5. Architectural Design: Asynchronous Queue & Observer Pattern

### 5.1 Architecture Diagram

```
[ User Calling Thread (TransRO) ]
      │
      ▼
1. Compute coordinates via geodetic math (Fast)
      │
2. If telemetry_logger:
   try:
       queue.put_nowait((direction, is_3d, lats, lons, metadata))  [< 0.001 ms pointer push]
   except queue.Full:
       pass  (Memory protection: drop telemetry if buffer full)
      │
3. Return coordinates to caller immediately! 🚀 (Zero latency, zero disk I/O)
      │
      └─────────────────────────────────────────────────────────┐
                                                                ▼
                                      [ Background Daemon Worker Thread ]
                                                                │
                                      1. Drain queue batch (pull available items)
                                      2. Check current record year; switch active DB if year changes
                                      3. Convert radians -> degrees & round(x, 2)
                                      4. Safely serialize metadata with json.dumps()
                                      5. Aggregate counts into in-memory dictionary
                                      6. Execute batched executemany() UPSERT into SQLite
```

### 5.2 Failure Isolation Guarantee
Inside `TransRO`:
```python
if self._collector is not None:
    try:
        self._collector.record(direction, lats, lons, is_3d=is_3d, metadata=metadata)
    except Exception:
        # Silently suppress any logging error so coordinate math NEVER fails
        pass
```

### 5.3 Batching & Performance Optimization
- To maximize throughput under load, the worker does not write single rows one by one.
- The worker drains all items currently available in the queue (or up to a batch size of 1,000), aggregates their counts in an in-memory dictionary `dict[key, count]`, and commits them in a single `executemany(...)` transaction.

### 5.4 Database Year Rotation
- The worker dynamically inspects the record timestamp/year.
- If a request belongs to a new calendar year (e.g., transition from 2026 to 2027), the worker closes the connection to `usage_2026.db` and opens `usage_2027.db`, ensuring continuous operation without restarting the server.

### 5.5 Metadata JSON Serialization
- The worker safely serializes `metadata` using `json.dumps(metadata)`.
- If `metadata` contains non-serializable objects, the worker catches `TypeError` and falls back to string representation or an empty dictionary, ensuring the thread never crashes on bad user input.

### 5.6 Memory Overflow Protection
- The in-memory `queue.Queue` is bounded by `maxsize=10000`.
- If the disk stalls and the queue fills up, `queue.put_nowait()` catches `queue.Full` and drops the telemetry payload silently, guaranteeing server RAM remains bounded and stable.

### 5.7 Lifecycle & Clean Shutdown
- The logger registers a cleanup function using Python's `atexit` standard library.
- Upon process termination, `close()` signals the worker thread and allows up to **2.0 seconds** to drain pending records to the SQLite database before terminating.

---

## 6. Developer Interface & API Ergonomics

### 6.1 Creating and Configuring TransRO
```python
from pytransdatro import TransRO
from pytransdatro.telemetry import SqliteUsageLogger

# 1. Default (no logging, zero overhead):
tr = TransRO()

# 2. Server with SQLite usage logging:
logger = SqliteUsageLogger(log_dir="./logs")
tr = TransRO(telemetry_logger=logger)
```

### 6.2 Executing Transformations with Metadata
```python
# Optional metadata parameter passed to transformation methods
lat, lon = tr.st70_to_etrs89(
    500000.0, 300000.0,
    metadata={"source": "web_api"}
)
```

---

## 7. Verification & Testing Strategy
All components can be thoroughly tested locally with pure Python `pytest`:
1. **Mathematical Invariant Test**: Run full test suite with `telemetry_logger=None` to verify zero regression in geodetic results.
2. **Coordinate Binning & Radian Conversion Test**: Verify known inputs (radians, both scalar and lists) are correctly converted to degrees and rounded to 2 decimal places.
3. **In-Place UPSERT Verification Test**: Perform 10 transformations in the same 1 km cell and assert that the database contains exactly **1 row** with `count = 10` and valid `year` and `week`.
4. **Failure Isolation Test**: Pass a mock logger that raises `RuntimeError("Simulated disk error")` and assert that `st70_to_etrs89` returns valid coordinates without raising.
5. **Drain on Close Test**: Push multiple batches, call `logger.close()`, and verify all records were committed to SQLite.
6. **Thread Affinity Test**: Verify that the SQLite connection is established and closed cleanly within the daemon thread without `ProgrammingError`.
