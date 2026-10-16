"""Populate SQLite usage database with 3 years of realistic geodetic dummy telemetry data.
Cleans existing out-of-grid records, preserves existing valid data, and simulates authentic
usage trends from November 2023 to October 2026.
"""

import datetime
import math
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
        print(f"Creating database backup: {BACKUP_PATH} ...")
        shutil.copy2(DB_PATH, BACKUP_PATH)
        print("Backup created successfully.")


def clean_corrupt_records(conn):
    """Remove records where coordinates fall outside Romania's geographic bounds."""
    cursor = conn.cursor()
    cursor.execute("""
        SELECT COUNT(*) FROM usage
        WHERE lat_deg < 43 OR lat_deg > 48 OR lon_deg < 20 OR lon_deg > 29;
    """)
    corrupt_count = cursor.fetchone()[0]
    print(f"\nFound {corrupt_count} out-of-grid / corrupted rows in 'usage' table.")

    if corrupt_count > 0:
        cursor.execute("""
            DELETE FROM usage
            WHERE lat_deg < 43 OR lat_deg > 48 OR lon_deg < 20 OR lon_deg > 29;
        """)
        conn.commit()
        print(f"Cleaned {corrupt_count} invalid rows from 'usage' table.")


def get_seasonal_multiplier(month):
    """Seasonal surveying activity factor in Romania."""
    # Low activity in winter freeze, peak in spring/autumn construction season
    factors = {
        1: 0.55,  # Jan: snow/frozen ground, minimal field surveying
        2: 0.65,  # Feb
        3: 1.15,  # Mar: spring construction starts
        4: 1.35,  # Apr: peak surveying
        5: 1.40,  # May: high activity
        6: 1.30,  # Jun
        7: 1.05,  # Jul
        8: 0.80,  # Aug: vacation season
        9: 1.35,  # Sep: autumn surveying surge
        10: 1.50,  # Oct: high cadastral submission volume
        11: 1.45,  # Nov: year-end deadline rush
        12: 0.65,  # Dec: holidays after Dec 15
    }
    return factors.get(month, 1.0)


def get_day_of_week_multiplier(weekday):
    """Workday vs Weekend factor (0=Mon, 6=Sun)."""
    if weekday < 5:
        # Monday to Friday: surveying office business hours
        return random.uniform(0.95, 1.25)
    elif weekday == 5:
        # Saturday: partial field work
        return random.uniform(0.18, 0.30)
    else:
        # Sunday: minimal / hobbyist
        return random.uniform(0.05, 0.15)


def generate_point_in_romania():
    """Generates a realistic (lat, lon) inside Romania based on population/hub weighting."""
    roll = random.random()
    if roll < 0.88:
        # 88% chance: pick from weighted metropolitan hubs
        weights = [h[4] for h in HUBS]
        chosen = random.choices(HUBS, weights=weights, k=1)[0]
        lat = random.gauss(chosen[0], chosen[2])
        lon = random.gauss(chosen[1], chosen[3])
    else:
        # 12% chance: dispersed rural cadastral work across national territory
        lat = random.uniform(43.70, 48.20)
        lon = random.uniform(20.40, 29.60)

    # Hard clamp to Romania's bounding box
    lat = max(43.60, min(48.25, lat))
    lon = max(20.25, min(29.65, lon))
    return lat, lon


def populate_telemetry(conn):
    """Generates 3 years of realistic daily stats and geographic usage."""
    cursor = conn.cursor()

    # Timeframe: Nov 1, 2023 to Oct 8, 2026 (prior to existing Oct 9-11 records)
    start_date = datetime.date(2023, 11, 1)
    end_date = datetime.date(2026, 10, 8)
    total_days = (end_date - start_date).days + 1

    print(
        f"\nGenerating 3 years of synthetic telemetry from {start_date} to {end_date} ({total_days:,} days)..."
    )

    daily_stats_data = []
    # (year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source) -> count
    monthly_usage_cells = {}

    current_date = start_date
    day_index = 0

    while current_date <= end_date:
        year = current_date.year
        month = current_date.month
        weekday = current_date.weekday()
        date_str = current_date.strftime("%Y-%m-%d")

        # 1. Macro Adoption Growth Factor (smooth S-curve adoption over 3 years)
        # Fraction of progress from day 0 to end (0.0 to 1.0)
        t = day_index / total_days
        # Base daily calls grow from ~80 (late 2023) to ~4,500 (late 2026)
        base_calls = 80 + (4400 * (t**1.8))

        # Apply seasonal and weekly multipliers
        season_mult = get_seasonal_multiplier(month)
        dow_mult = get_day_of_week_multiplier(weekday)

        # Occasional random deadline spikes (e.g. quarterly cadastral submissions)
        spike_mult = 1.0
        if random.random() < 0.04 and weekday < 5:
            spike_mult = random.uniform(1.8, 3.2)

        day_total_calls = int(base_calls * season_mult * dow_mult * spike_mult)
        day_total_calls = max(5, day_total_calls)

        # Source volume distributions over adoption lifecycle:
        # 2023: desktop_cli (40%), rest_api (50%), web_map (10% - early beta)
        # 2024: web_map rises to 35%, rest_api 45%, desktop_cli 20%
        # 2025-2026: web_map (45%), rest_api (45%), desktop_cli (10%)
        if year == 2023:
            s1_share, s2_share, s3_share = 0.10, 0.50, 0.40
        elif year == 2024:
            s1_share, s2_share, s3_share = 0.35, 0.45, 0.20
        else:
            s1_share, s2_share, s3_share = 0.45, 0.45, 0.10

        s1_calls = max(1, int(day_total_calls * s1_share))
        s2_calls = max(1, int(day_total_calls * s2_share))
        s3_calls = max(1, int(day_total_calls * s3_share))

        # Points per call profiles:
        # Source 1 (web_map): interactive single points or small lists (1 - 8 pts/call)
        s1_points = int(s1_calls * random.uniform(1.2, 4.5))
        # Source 2 (rest_api): moderate API batches (5 - 35 pts/call)
        s2_points = int(s2_calls * random.uniform(8.0, 28.0))
        # Source 3 (desktop_cli): heavy cadastral file conversions (100 - 3,500 pts/call)
        # Occasional bulk file spike (e.g. 10,000 to 40,000 points in a single file)
        bulk_bonus = 0
        if random.random() < 0.10:
            bulk_bonus = random.randint(5000, 35000)
        s3_points = int(s3_calls * random.uniform(40.0, 300.0)) + bulk_bonus

        # Append to daily_stats
        daily_stats_data.append((date_str, 1, s1_calls, s1_points))
        daily_stats_data.append((date_str, 2, s2_calls, s2_points))
        daily_stats_data.append((date_str, 3, s3_calls, s3_points))

        # Sample points to generate spatial usage cells for this day
        # To keep DB size compact (~3.5-5 MB), sample representative points
        # Sample ~10 to 40 distinct cells per day, weighted by day volume
        sample_pts = min(45, max(8, int(math.sqrt(day_total_calls))))
        for _ in range(sample_pts):
            lat, lon = generate_point_in_romania()
            lat_deg = int(lat)
            lat_grad = int(round((lat - lat_deg) * 100)) % 100
            lon_deg = int(lon)
            lon_grad = int(round((lon - lon_deg) * 100)) % 100

            # 65% Stereo70 -> ETRS89, 35% ETRS89 -> Stereo70
            direction = 1 if random.random() < 0.65 else 0
            # 25% 3D, 75% 2D
            is_3d = 1 if random.random() < 0.25 else 0
            # Distribute sources
            src = random.choices([1, 2, 3], weights=[s1_share, s2_share, s3_share], k=1)[0]
            # Hits for this spatial cell
            hit_weight = random.randint(1, 15)

            key = (year, month, lat_deg, lat_grad, lon_deg, lon_grad, direction, is_3d, src)
            monthly_usage_cells[key] = monthly_usage_cells.get(key, 0) + hit_weight

        current_date += datetime.timedelta(days=1)
        day_index += 1

    # Insert daily_stats in batch
    print(f"Inserting {len(daily_stats_data):,} rows into 'daily_stats' table...")
    daily_insert_sql = """
        INSERT INTO daily_stats (date, source, count_calls, count_points)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(date, source) DO UPDATE SET
            count_calls = count_calls + excluded.count_calls,
            count_points = count_points + excluded.count_points;
    """
    cursor.executemany(daily_insert_sql, daily_stats_data)
    conn.commit()

    # Insert usage cells in batch
    print(f"Inserting {len(monthly_usage_cells):,} aggregated 1 km cells into 'usage' table...")
    usage_data = [
        (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], count)
        for k, count in monthly_usage_cells.items()
    ]
    usage_insert_sql = """
        INSERT INTO usage (year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source)
        DO UPDATE SET count = count + excluded.count;
    """
    # Execute in chunks of 5000 to prevent locking
    chunk_size = 5000
    for i in range(0, len(usage_data), chunk_size):
        cursor.executemany(usage_insert_sql, usage_data[i : i + chunk_size])
        conn.commit()

    print("Data insertion completed successfully.")


def verify_database(conn):
    """Verifies bounds and prints summary metrics."""
    cursor = conn.cursor()

    # Verify corrupt rows count
    cursor.execute("""
        SELECT COUNT(*) FROM usage
        WHERE lat_deg < 43 OR lat_deg > 48 OR lon_deg < 20 OR lon_deg > 29;
    """)
    corrupt = cursor.fetchone()[0]

    cursor.execute("SELECT MIN(date), MAX(date), COUNT(*) FROM daily_stats;")
    min_date, max_date, total_days = cursor.fetchone()

    cursor.execute("SELECT COUNT(*), SUM(count) FROM usage;")
    total_cells, total_hits = cursor.fetchone()

    cursor.execute("SELECT SUM(count_calls), SUM(count_points) FROM daily_stats;")
    total_calls, total_points = cursor.fetchone()

    print("\n" + "=" * 70)
    print("DATABASE VERIFICATION SUMMARY")
    print("=" * 70)
    print(f"Corrupt/Out-of-grid Rows in 'usage':   {corrupt} (Must be 0)")
    print(f"Daily Stats Date Range:                {min_date} to {max_date}")
    print(f"Total Daily Stat Records:              {total_days:,}")
    print(f"Total Calls Logged:                    {total_calls:,}")
    print(f"Total Points Logged:                   {total_points:,}")
    print(f"Total Unique Aggregated ~1km Cells:    {total_cells:,}")
    print(f"Total Spatial Hits in Heatmap:         {total_hits:,}")
    print("=" * 70)


def main():
    if not DB_PATH.exists():
        print(f"Database file not found: {DB_PATH}")
        return

    backup_database()

    conn = sqlite3.connect(str(DB_PATH), timeout=10.0)
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")

    try:
        clean_corrupt_records(conn)
        populate_telemetry(conn)
        verify_database(conn)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
