# PyTransDatRO: Usage Telemetry & Logging

## 1. Overview & Objectives
`PyTransDatRO` includes an optional, high-performance telemetry and logging system designed to track coordinate transformation activity over time (e.g. operational usage trends and spatial density heatmaps).

### Key Architectural Principles
- **Disabled by Default (Zero Overhead)**: When used as an embedded desktop library or in local scripts, logging is completely inactive by default and produces zero disk I/O and zero latency overhead.
- **Asynchronous & Non-Blocking**: Coordinate calculation threads push pointer references to an in-memory queue (`< 0.001 ms`). All disk I/O, binning, and database writes occur in an independent background daemon thread.
- **Failure Isolation (Non-Fatal)**: Coordinate mathematics will **never** fail or raise exceptions due to telemetry errors (e.g. disk full, permission denied, database locks).
- **Privacy & Anonymization**: Exact cadastral coordinates are never stored on disk; they are binned into ~1 km geographic cells.
- **Pure Python Standard Library**: Relies solely on Python standard library modules (`sqlite3`, `queue`, `threading`, `datetime`, `math`, `pathlib`, `atexit`). Zero third-party dependencies.

---

## 2. Storage Architecture (`logs/usage.db`)

All logging data is written to a single SQLite database: `logs/usage.db` configured with **Write-Ahead Logging (WAL)** and busy timeouts to ensure multi-process concurrency on web servers.

The database contains three relational tables:

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
- `1`: `web_map` (Web Map UI client)
- `2`: `rest_api` (Direct REST API / backend microservice)
- `3`: `desktop_cli` (Desktop GIS application / CLI plugin)

### 2.2 Daily Volume Metrics Table (`daily_stats`)
Tracks operational volume per calendar day grouped by client source:
```sql
CREATE TABLE IF NOT EXISTS daily_stats (
    date TEXT NOT NULL,
    source INTEGER NOT NULL,
    count_calls INTEGER NOT NULL DEFAULT 1,
    count_points INTEGER NOT NULL DEFAULT 1,
    PRIMARY KEY (date, source)
);
CREATE INDEX IF NOT EXISTS idx_daily_stats_date ON daily_stats(date);
```
- **`count_calls`**: Number of transformation function/API calls made that day.
- **`count_points`**: Total number of individual coordinate points processed that day.

### 2.3 Spatial Usage Table (`usage`)
Stores geographic heatmap data grouped by year, month, ~1 km cell, direction, dimension, and source:
```sql
CREATE TABLE IF NOT EXISTS usage (
    year INTEGER NOT NULL,
    month INTEGER NOT NULL,
    lat_deg INTEGER NOT NULL,
    lat_grad INTEGER NOT NULL,
    lon_deg INTEGER NOT NULL,
    lon_grad INTEGER NOT NULL,
    is_st70_to_etrs89 INTEGER NOT NULL,
    is_3d INTEGER NOT NULL,
    source INTEGER NOT NULL DEFAULT 0,
    count INTEGER NOT NULL DEFAULT 1,
    PRIMARY KEY (year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source)
);
CREATE INDEX IF NOT EXISTS idx_usage_year_month ON usage(year, month);
```

#### Coordinate Splitting (Degrees & Grads)
To minimize disk space, coordinates are converted to degrees, rounded to 2 decimal places (~1 km resolution), and split into integer degrees and hundredths (grads):
- Example: Latitude `45.79°` $\rightarrow$ `lat_deg = 45`, `lat_grad = 79`
- Reconstructing decimal degrees: `lat = lat_deg + (lat_grad / 100.0)`
- In SQLite, small integers ($\le 127$) consume only **1 byte** per value, compared to 8 bytes for standard IEEE floating-point numbers.

---

## 3. Storage Footprint & Disk Space Implications

> [!NOTE]
> **Storage Footprint Remark**:
> A full month's worth of active usage data covering the entire territory of Romania across all transformation directions and operations requires approximately **3.2 MB** of disk space.

### Why the Footprint is Minimal
1. **In-Place Aggregation (UPSERT)**: Any additional transformation occurring in an already visited 1 km cell during the same month simply increments the `count` column in-place:
   ```sql
   ON CONFLICT(...) DO UPDATE SET count = count + excluded.count;
   ```
   It does **not** create a new row or consume additional disk space.
2. **Compact All-Integer Schema**: Columns such as `is_st70_to_etrs89`, `is_3d`, `source`, `month`, and the split coordinate fields are stored as 0-byte or 1-byte variable-length integers.
3. **No Heavy String Metadata**: Client identities are tracked via low-cardinality integer IDs rather than repeated JSON strings.

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

A built-in inspection tool is provided in `tools/inspect_telemetry_db.py`:
```bash
python tools/inspect_telemetry_db.py
```

### Useful SQL Queries

#### Monthly Traffic by Client Source
```sql
SELECT s.name AS source_name, SUM(u.count) AS total_coordinates
FROM usage u
LEFT JOIN sources s ON u.source = s.id
WHERE u.year = 2026 AND u.month = 10
GROUP BY u.source;
```

#### Daily Call & Point Volumes
```sql
SELECT date, s.name AS source_name, count_calls, count_points
FROM daily_stats d
LEFT JOIN sources s ON d.source = s.id
ORDER BY date DESC, source;
```

#### Top Geographic Hotspots (~1 km Cells)
```sql
SELECT 
    printf('%02d.%02d', lat_deg, lat_grad) AS lat,
    printf('%02d.%02d', lon_deg, lon_grad) AS lon,
    count
FROM usage
ORDER BY count DESC
LIMIT 10;
```
