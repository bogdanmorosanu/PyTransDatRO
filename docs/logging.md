# PyTransDatRO: Usage Telemetry & Logging

## 1. Overview & Objectives
`PyTransDatRO` includes an optional, high-performance telemetry and logging system designed to track coordinate transformation activity over time (e.g. operational usage trends and spatial density heatmaps).

### Key Architectural Principles
- **Disabled by Default (Zero Overhead)**: When used as an embedded desktop library or in local scripts, logging is completely inactive by default and produces zero disk I/O and zero latency overhead.
- **Asynchronous & Non-Blocking**: Coordinate calculation threads push pointer references to an in-memory queue (`< 0.001 ms`). All disk I/O, binning, and database writes occur in an independent background daemon thread.
- **Failure Isolation (Non-Fatal)**: Coordinate mathematics will **never** fail or raise exceptions due to telemetry errors (e.g. disk full, permission denied, database locks).
- **Privacy & Anonymization**: Exact cadastral coordinates are never stored on disk; they are binned into ~1 km geographic cells.
- **Pure Python Standard Library**: Relies solely on Python standard library modules (`sqlite3`, `queue`, `threading`, `datetime`, `calendar`, `math`, `pathlib`, `atexit`). Zero third-party dependencies.

---

## 2. Storage Architecture (`logs/usage.db`)

All logging data is written to a single SQLite database: `logs/usage.db` configured with **Write-Ahead Logging (WAL)** and busy timeouts to ensure multi-process concurrency on web servers.

The database contains five relational tables:

### 2.1 Sources Dimension Table (`sources`)
Maps client source IDs to human-readable names.
```sql
CREATE TABLE IF NOT EXISTS sources (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT
);
```
Predefined defaults:
- `0`: `unknown` (Default / unspecified client or standalone script)
- `1`: `rest_api` (Direct REST API / backend microservice)
- `2`: `web_map` (Web Map UI client)
- `3`: `desktop_cli` (Desktop GIS application / CLI plugin)

### 2.2 Operations Dimension Table (`operations`)
Maps transformation direction IDs to standardized operation identifiers:
```sql
CREATE TABLE IF NOT EXISTS operations (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT
);
```
Predefined defaults:
- `1`: `st70_to_etrs89` (Stereo 70 to ETRS89)
- `2`: `etrs89_to_st70` (ETRS89 to Stereo 70)

### 2.3 Window Configuration Registry Table (`telemetry_config`)
Maintains dynamic time window definitions and their validity periods (SCD Type 2 interval registry):
```sql
CREATE TABLE IF NOT EXISTS telemetry_config (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    metric_type TEXT NOT NULL,
    window_days INTEGER NOT NULL,
    anchor_date TEXT NOT NULL,
    valid_from TEXT NOT NULL,
    valid_to TEXT,
    description TEXT
);
```
Default entries:
- `computation`: `window_days = 1` (daily volume metrics).
- `spatial`: `window_days = 30` (monthly approximated spatial heatmaps, or configurable to 7, 5, etc.).

### 2.4 Computation Volume Metrics Table (`computation_usage`)
Tracks operational volume (invocations and transformed coordinates) per time window, operation, dimensionality, and source:
```sql
CREATE TABLE IF NOT EXISTS computation_usage (
    period_start TEXT NOT NULL,
    period_end TEXT NOT NULL,
    operation INTEGER NOT NULL,
    is_3d INTEGER NOT NULL,
    source INTEGER NOT NULL DEFAULT 0,
    count_calls INTEGER NOT NULL DEFAULT 1,
    count_points INTEGER NOT NULL DEFAULT 1,
    PRIMARY KEY (period_start, operation, is_3d, source)
);
CREATE INDEX IF NOT EXISTS idx_comp_usage_period ON computation_usage(period_start);
```
- **`period_start` / `period_end`**: Inclusive ISO date boundaries (`YYYY-MM-DD`). For daily granularity, start and end dates are identical.
- **`operation`**: Foreign key to `operations` (`1` or `2`).
- **`is_3d`**: `0` (planimetric 2D) or `1` (3D elevation/quasigeoid height).
- **`count_calls`**: Number of transformation function/API calls made during the window.
- **`count_points`**: Total number of individual coordinate points processed during the window.

### 2.5 Spatial Usage Table (`spatial_usage`)
Stores geographic heatmap data grouped by dynamic period window, ~1 km cell, operation, dimension, and source:
```sql
CREATE TABLE IF NOT EXISTS spatial_usage (
    period_start TEXT NOT NULL,
    period_end TEXT NOT NULL,
    lat_deg INTEGER NOT NULL,
    lat_grad INTEGER NOT NULL,
    lon_deg INTEGER NOT NULL,
    lon_grad INTEGER NOT NULL,
    operation INTEGER NOT NULL,
    is_3d INTEGER NOT NULL,
    source INTEGER NOT NULL DEFAULT 0,
    count INTEGER NOT NULL DEFAULT 1,
    PRIMARY KEY (period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
);
CREATE INDEX IF NOT EXISTS idx_spatial_usage_period ON spatial_usage(period_start);
```

#### Coordinate Splitting (Degrees & Grads)
To minimize disk space, coordinates are converted to degrees, rounded to 2 decimal places (~1 km resolution), and split into integer degrees and hundredths (grads):
- Example: Latitude `45.79°` $\rightarrow$ `lat_deg = 45`, `lat_grad = 79`
- Reconstructing decimal degrees: `lat = lat_deg + (lat_grad / 100.0)`
- In SQLite, small integers ($\le 127$) consume only **1 byte** per value, compared to 8 bytes for standard IEEE floating-point numbers.

#### Macroscopic $0.1^\circ$ Aggregation & Visual Access
The REST API endpoint `GET /api/v1/telemetry/spatial-matrix` dynamically aggregates the ~1 km stored bins (`FLOOR(lat_grad / 10)`, `FLOOR(lon_grad / 10)`) into macroscopic $0.1^\circ$ ($\approx 10\text{ km} \times 7.5\text{ km}$) tiles for punchcard heatmaps and calculates mathematical quadrant volumes (NW, NE, SW, SE).

---

## 3. Storage Footprint & Dynamic Windows

> [!NOTE]
> **Self-Describing Windows**:
> By storing both `period_start` and `period_end` directly in the record, each row remains self-contained. Even if window granularity is adjusted from 30 days to 7 days or 5 days, historical rows preserve their exact duration and can be queried or normalized without ambiguity.

### In-Place Aggregation (UPSERT)
Any additional transformation occurring in an already visited 1 km cell during the same active window increments `count` in-place:
```sql
ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
DO UPDATE SET count = count + excluded.count;
```

---

## 4. Developer Interface & Usage

### 4.1 Enabling Telemetry in `TransRO`
```python
from pytransdatro import TransRO
from pytransdatro.telemetry import SqliteUsageLogger

# 1. Default (no logging, zero overhead):
tr = TransRO()

# 2. Server deployment with telemetry enabled:
logger = SqliteUsageLogger(log_dir="./logs")
tr = TransRO(telemetry_logger=logger)
```

### 4.2 Executing Transformations with Source IDs
Pass the optional `source` integer parameter:
```python
# Call from Web Map (source=1)
lat, lon = tr.st70_to_etrs89(500000.0, 300000.0, source=1)

# Call from REST API (source=2)
n, e = tr.etrs89_to_st70(lat, lon, source=2)
```

### 4.3 Clean Shutdown
`SqliteUsageLogger` registers an automatic `atexit` hook. When the host process exits, pending items in the queue are flushed to SQLite with a 2-second timeout. You can also trigger an explicit flush and close:
```python
logger.close()
```

---

## 5. Inspecting the Database

### Useful SQL Queries

#### Period Traffic by Client Source and Operation
```sql
SELECT 
    c.period_start,
    c.period_end,
    o.name AS operation_name,
    s.name AS source_name,
    c.is_3d,
    c.count_calls,
    c.count_points
FROM computation_usage c
LEFT JOIN operations o ON c.operation = o.id
LEFT JOIN sources s ON c.source = s.id
ORDER BY c.period_start DESC, c.source;
```

#### Top Geographic Hotspots (~1 km Cells)
```sql
SELECT 
    printf('%02d.%02d', lat_deg, lat_grad) AS lat,
    printf('%02d.%02d', lon_deg, lon_grad) AS lon,
    o.name AS operation_name,
    u.period_start,
    u.period_end,
    count
FROM spatial_usage u
LEFT JOIN operations o ON u.operation = o.id
ORDER BY count DESC
LIMIT 10;
```
