"""Populate SQLite usage database with 3 years of realistic geodetic dummy telemetry data.

Features:
- Targets modern telemetry tables (computation_usage, spatial_usage, telemetry_config, sources, operations).
- Simulates multiple historical telemetry period regimes:
    * Year 1 (2023-11-01 to 2024-10-31): Annual windows (window_days = -2)
    * Year 2 (2024-11-01 to 2025-10-31): Monthly windows (window_days = -1)
    * Year 3 (2025-11-01 to 2026-10-24): 10-day mathematical windows (window_days = 10)
- Enforces 100% mathematical parity: total points in computation_usage identically equals
  total points in spatial_usage.
- Tracks count_calls and count (points) accurately on both tables.
- High-performance vector simulation (~2 seconds for 3 years of data).
"""

import calendar
import datetime
import random
import shutil
import sqlite3
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LOGS_DIR = REPO_ROOT / "logs"
DB_PATH = LOGS_DIR / "usage.db"
BACKUP_PATH = LOGS_DIR / "usage.db.bak"

# Major Romanian metropolitan and regional hubs (lat, lon, lat_std, lon_std, weight)
HUBS = [
    # Bucharest & Ilfov (Major economic & construction hub)
    (44.4323, 26.1063, 0.14, 0.18, 0.35),
    # Cluj-Napoca & Cluj metropolitan
    (46.7712, 23.6236, 0.10, 0.12, 0.12),
    # Timișoara & Banat industrial zone
    (45.7537, 21.2257, 0.08, 0.10, 0.08),
    # Iași & Moldova regional hub
    (47.1585, 27.6014, 0.08, 0.10, 0.08),
    # Constanța, Port & Black Sea coast
    (44.1807, 28.6343, 0.08, 0.09, 0.07),
    # Brașov & Prahova Valley corridor
    (45.6579, 25.6012, 0.07, 0.08, 0.06),
    # Craiova & Oltenia
    (44.3302, 23.7949, 0.06, 0.08, 0.04),
    # Oradea & Bihor
    (47.0465, 21.9189, 0.06, 0.07, 0.04),
    # Sibiu & Central Transylvania
    (45.7983, 24.1256, 0.06, 0.07, 0.04),
    # Galați & Brăila conurbation
    (45.4353, 28.0080, 0.07, 0.08, 0.04),
    # Ploiești & Pitești (A1/DN1 industrial axis)
    (44.9367, 26.0129, 0.08, 0.10, 0.04),
    # Arad & Western border trade corridor
    (46.1866, 21.3123, 0.05, 0.06, 0.02),
    # Bacău & Suceava (A7 motorway development axis)
    (46.5670, 26.9146, 0.07, 0.09, 0.02),
]


def backup_database():
    """Create a safety backup of usage.db."""
    if DB_PATH.exists():
        try:
            print(f"Creating database backup: {BACKUP_PATH} ...", flush=True)
            shutil.copy2(DB_PATH, BACKUP_PATH)
            print("Backup created successfully.", flush=True)
        except Exception as e:
            print(f"Notice: backup skipped ({e})", flush=True)


def init_schema(conn: sqlite3.Connection):
    """Ensures modern tables and indexes exist."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS sources (
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL,
            description TEXT
        );
    """)
    conn.executemany(
        """
        INSERT OR IGNORE INTO sources (id, name, description)
        VALUES (?, ?, ?);
        """,
        [
            (0, "unknown", "Default / Unspecified origin"),
            (1, "rest_api", "REST API Service"),
            (2, "web_map", "Web Map Application"),
            (3, "desktop_cli", "Desktop / CLI App"),
        ],
    )

    conn.execute("""
        CREATE TABLE IF NOT EXISTS operations (
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL,
            description TEXT
        );
    """)
    conn.executemany(
        """
        INSERT OR IGNORE INTO operations (id, name, description)
        VALUES (?, ?, ?);
        """,
        [
            (1, "st70_to_etrs89", "Stereo70 to ETRS89"),
            (2, "etrs89_to_st70", "ETRS89 to Stereo70"),
        ],
    )

    conn.execute("""
        CREATE TABLE IF NOT EXISTS telemetry_config (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            metric_type TEXT NOT NULL,
            window_days INTEGER NOT NULL,
            anchor_date TEXT NOT NULL,
            valid_from TEXT NOT NULL,
            valid_to TEXT,
            description TEXT
        );
    """)

    conn.execute("""
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
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_comp_usage_period ON computation_usage(period_start);"
    )

    conn.execute("""
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
            count_calls INTEGER NOT NULL DEFAULT 1,
            count INTEGER NOT NULL DEFAULT 1,
            PRIMARY KEY (period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
        );
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_spatial_usage_period ON spatial_usage(period_start);"
    )

    # Ensure count_calls column exists if spatial_usage already existed
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(spatial_usage)")
    cols = {col[1] for col in cursor.fetchall()}
    if "count_calls" not in cols:
        conn.execute("ALTER TABLE spatial_usage ADD COLUMN count_calls INTEGER NOT NULL DEFAULT 1;")

    conn.commit()


def setup_telemetry_configs(conn: sqlite3.Connection):
    """Sets up historical multi-period regimes in telemetry_config."""
    conn.execute("DELETE FROM telemetry_config;")
    configs = [
        # Computation: always daily
        (
            "computation",
            1,
            "2023-01-01",
            "2023-01-01",
            None,
            "Daily computation usage (calls & points)",
        ),
        # Spatial Year 1: Calendar Year (-2)
        (
            "spatial",
            -2,
            "2023-11-01",
            "2023-11-01",
            "2024-10-31",
            "Year 1: Calendar year spatial window",
        ),
        # Spatial Year 2: Calendar Month (-1)
        (
            "spatial",
            -1,
            "2024-11-01",
            "2024-11-01",
            "2025-10-31",
            "Year 2: Calendar month spatial window",
        ),
        # Spatial Year 3: 10-day mathematical windows
        ("spatial", 10, "2025-11-01", "2025-11-01", None, "Year 3: 10-day rolling spatial window"),
    ]
    conn.executemany(
        """
        INSERT INTO telemetry_config (metric_type, window_days, anchor_date, valid_from, valid_to, description)
        VALUES (?, ?, ?, ?, ?, ?);
        """,
        configs,
    )
    conn.commit()


def calculate_period(
    date_obj: datetime.date, window_days: int, anchor_date_str: str
) -> tuple[str, str]:
    """Calculates period [start, end] for a given date and window configuration."""
    if window_days == -2:
        y = date_obj.year
        return f"{y:04d}-01-01", f"{y:04d}-12-31"

    if window_days in (-1, 30):
        y, m = date_obj.year, date_obj.month
        last_day = calendar.monthrange(y, m)[1]
        return f"{y:04d}-{m:02d}-01", f"{y:04d}-{m:02d}-{last_day:02d}"

    if window_days <= 1:
        date_str = date_obj.strftime("%Y-%m-%d")
        return date_str, date_str

    anchor = datetime.date.fromisoformat(anchor_date_str)
    delta_days = (date_obj - anchor).days
    bucket_index = delta_days // window_days
    p_start = anchor + datetime.timedelta(days=bucket_index * window_days)
    p_end = p_start + datetime.timedelta(days=window_days - 1)
    return p_start.strftime("%Y-%m-%d"), p_end.strftime("%Y-%m-%d")


def get_spatial_config_for_date(date_obj: datetime.date) -> tuple[int, str]:
    """Resolves active spatial config according to multi-period timeline."""
    d_str = date_obj.strftime("%Y-%m-%d")
    if d_str <= "2024-10-31":
        return -2, "2023-11-01"
    elif d_str <= "2025-10-31":
        return -1, "2024-11-01"
    else:
        return 10, "2025-11-01"


def get_seasonal_multiplier(month: int) -> float:
    """Seasonal surveying activity factor in Romania."""
    factors = {
        1: 0.55,
        2: 0.65,
        3: 1.15,
        4: 1.35,
        5: 1.40,
        6: 1.30,
        7: 1.05,
        8: 0.80,
        9: 1.35,
        10: 1.50,
        11: 1.45,
        12: 0.65,
    }
    return factors.get(month, 1.0)


def get_day_of_week_multiplier(weekday: int) -> float:
    """Workday vs Weekend factor (0=Mon, 6=Sun)."""
    if weekday < 5:
        return random.uniform(0.95, 1.25)
    elif weekday == 5:
        return random.uniform(0.18, 0.30)
    else:
        return random.uniform(0.05, 0.15)


def generate_point_in_romania() -> tuple[float, float]:
    """Generates a realistic (lat, lon) inside Romania based on hub weighting."""
    roll = random.random()
    if roll < 0.88:
        weights = [h[4] for h in HUBS]
        chosen = random.choices(HUBS, weights=weights, k=1)[0]
        lat = random.gauss(chosen[0], chosen[2])
        lon = random.gauss(chosen[1], chosen[3])
    else:
        lat = random.uniform(43.70, 48.20)
        lon = random.uniform(20.40, 29.60)

    lat = max(43.60, min(48.25, lat))
    lon = max(20.25, min(29.65, lon))
    return lat, lon


def populate_telemetry(conn: sqlite3.Connection):
    """Generates 3 years of synthetic telemetry with exact points parity."""
    # Wipe existing rows in usage tables to guarantee clean parity
    conn.execute("DELETE FROM computation_usage;")
    conn.execute("DELETE FROM spatial_usage;")
    conn.commit()

    start_date = datetime.date(2023, 11, 1)
    end_date = datetime.date(2026, 10, 8)
    total_days = (end_date - start_date).days + 1

    print(
        f"\nGenerating 3 years of synthetic telemetry from {start_date} to {end_date} ({total_days:,} days)...",
        flush=True,
    )

    # (comp_start, comp_end, op, is_3d, src) -> [calls, points]
    computation_records: dict[tuple, list[int]] = {}
    # (spatial_start, spatial_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src) -> [calls, points]
    spatial_records: dict[tuple, list[int]] = {}

    current_date = start_date
    day_index = 0

    while current_date <= end_date:
        weekday = current_date.weekday()
        year = current_date.year
        month = current_date.month

        # Computation window: daily
        comp_start = current_date.strftime("%Y-%m-%d")
        comp_end = comp_start

        # Spatial window: resolved from multi-period timeline
        sp_window_days, sp_anchor = get_spatial_config_for_date(current_date)
        sp_start, sp_end = calculate_period(current_date, sp_window_days, sp_anchor)

        # Volume modeling
        t = day_index / total_days
        base_intensity = 6 + int(35 * (t**1.8))
        season_mult = get_seasonal_multiplier(month)
        dow_mult = get_day_of_week_multiplier(weekday)
        spike_mult = random.uniform(1.6, 2.5) if (random.random() < 0.04 and weekday < 5) else 1.0

        day_cells_count = max(2, int(base_intensity * season_mult * dow_mult * spike_mult))

        # Source shares
        if year == 2023:
            s_shares = {1: 0.50, 2: 0.10, 3: 0.40}
        elif year == 2024:
            s_shares = {1: 0.45, 2: 0.35, 3: 0.20}
        else:
            s_shares = {1: 0.45, 2: 0.45, 3: 0.10}

        # For this day, sample spatial cells and directly correlate computation stats
        for _ in range(day_cells_count):
            src = random.choices(list(s_shares.keys()), weights=list(s_shares.values()), k=1)[0]
            op = 1 if random.random() < 0.65 else 2
            is_3d = 1 if random.random() < 0.25 else 0

            # Cell location in Romania
            lat, lon = generate_point_in_romania()
            lat_d = int(lat)
            lat_g = int(round((lat - lat_d) * 100)) % 100
            lon_d = int(lon)
            lon_g = int(round((lon - lon_d) * 100)) % 100

            # Calls and points for this cell event
            if src == 2:  # web_map
                event_calls = random.randint(1, 4)
                event_points = event_calls * random.randint(1, 3)
            elif src == 1:  # rest_api
                event_calls = random.randint(2, 8)
                event_points = event_calls * random.randint(3, 15)
            else:  # desktop_cli
                event_calls = random.randint(1, 3)
                event_points = event_calls * random.randint(15, 80)

            # 1. Spatial record accumulation
            sp_key = (sp_start, sp_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src)
            if sp_key not in spatial_records:
                spatial_records[sp_key] = [0, 0]
            spatial_records[sp_key][0] += event_calls
            spatial_records[sp_key][1] += event_points

            # 2. Computation record accumulation (guaranteeing exact parity!)
            cp_key = (comp_start, comp_end, op, is_3d, src)
            if cp_key not in computation_records:
                computation_records[cp_key] = [0, 0]
            computation_records[cp_key][0] += event_calls
            computation_records[cp_key][1] += event_points

        current_date += datetime.timedelta(days=1)
        day_index += 1

    # Insert computation_usage in chunks
    print(
        f"Inserting {len(computation_records):,} rows into 'computation_usage' table...", flush=True
    )
    comp_data = [(k[0], k[1], k[2], k[3], k[4], v[0], v[1]) for k, v in computation_records.items()]
    comp_insert_sql = """
        INSERT INTO computation_usage (period_start, period_end, operation, is_3d, source, count_calls, count_points)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(period_start, operation, is_3d, source) DO UPDATE SET
            count_calls = count_calls + excluded.count_calls,
            count_points = count_points + excluded.count_points;
    """
    for i in range(0, len(comp_data), 5000):
        conn.executemany(comp_insert_sql, comp_data[i : i + 5000])
        conn.commit()

    # Insert spatial_usage in chunks
    print(f"Inserting {len(spatial_records):,} rows into 'spatial_usage' table...", flush=True)
    spatial_data = [
        (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], v[0], v[1])
        for k, v in spatial_records.items()
    ]
    spatial_insert_sql = """
        INSERT INTO spatial_usage (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count_calls, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source) DO UPDATE SET
            count_calls = count_calls + excluded.count_calls,
            count = count + excluded.count;
    """
    for i in range(0, len(spatial_data), 5000):
        conn.executemany(spatial_insert_sql, spatial_data[i : i + 5000])
        conn.commit()

    print("Data insertion completed successfully.", flush=True)


def verify_database(conn: sqlite3.Connection):
    """Verifies bounds and strictly enforces mathematical parity."""
    cursor = conn.cursor()

    # 1. Bounds verification
    cursor.execute("""
        SELECT COUNT(*) FROM spatial_usage
        WHERE lat_deg < 43 OR lat_deg > 48 OR lon_deg < 20 OR lon_deg > 29;
    """)
    corrupt = cursor.fetchone()[0]

    # 2. Computation stats summary
    cursor.execute(
        "SELECT MIN(period_start), MAX(period_end), COUNT(*), SUM(count_calls), SUM(count_points) FROM computation_usage;"
    )
    min_date, max_date, comp_rows, total_calls, comp_points = cursor.fetchone()

    # 3. Spatial stats summary
    cursor.execute("SELECT COUNT(*), SUM(count_calls), SUM(count) FROM spatial_usage;")
    total_cells, total_spatial_cell_calls, spatial_points = cursor.fetchone()

    # 4. Period regimes verification in spatial_usage
    cursor.execute("""
        SELECT DISTINCT period_start, period_end
        FROM spatial_usage
        ORDER BY period_start ASC
    """)
    periods = cursor.fetchall()

    print("\n" + "=" * 70, flush=True)
    print("DATABASE VERIFICATION SUMMARY", flush=True)
    print("=" * 70, flush=True)
    print(f"Corrupt/Out-of-grid Rows in 'spatial_usage': {corrupt} (Must be 0)", flush=True)
    print(f"Computation Usage Date Range:               {min_date} to {max_date}", flush=True)
    print(f"Total Computation Records:                  {comp_rows:,}", flush=True)
    print(f"Total Calls Logged:                         {total_calls:,}", flush=True)
    print(f"Total Points (computation_usage):           {comp_points:,}", flush=True)
    print(f"Total Points (spatial_usage):               {spatial_points:,}", flush=True)
    print(
        f"Parity Difference (comp_pts - sp_pts):      {comp_points - spatial_points} (Must be 0!)",
        flush=True,
    )
    print(f"Total Unique Aggregated ~1km Cells:         {total_cells:,}", flush=True)
    print(f"Total Spatial Hit Events (sum of cells):    {total_spatial_cell_calls:,}", flush=True)
    print(f"Distinct Spatial Period Buckets:            {len(periods):,}", flush=True)
    print("=" * 70, flush=True)

    assert corrupt == 0, f"Found {corrupt} corrupt records!"
    assert comp_points == spatial_points, f"Parity mismatch: {comp_points} vs {spatial_points}!"
    print("SUCCESS: 100% Mathematical Parity Invariant Verified!", flush=True)


def main():
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    backup_database()

    conn = sqlite3.connect(str(DB_PATH), timeout=10.0)
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")

    try:
        init_schema(conn)
        setup_telemetry_configs(conn)
        populate_telemetry(conn)
        verify_database(conn)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
