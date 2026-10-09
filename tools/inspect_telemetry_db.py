"""Inspect and display the contents of the PyTransDat usage telemetry SQLite database.
Shows sources lookup, daily call/points volume stats, and geographic usage cells.
"""

import sys
import sqlite3
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LOGS_DIR = REPO_ROOT / "logs"


def format_table(headers, rows):
    if not rows:
        return "  (No records found)"
    col_widths = [len(h) for h in headers]
    str_rows = []
    for row in rows:
        str_row = [str(val) if val is not None else "NULL" for val in row]
        for i, val in enumerate(str_row):
            col_widths[i] = max(col_widths[i], len(val))
        str_rows.append(str_row)

    header_line = " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    sep_line = "-+-".join("-" * col_widths[i] for i in range(len(headers)))
    content_lines = [" | ".join(val.ljust(col_widths[i]) for i, val in enumerate(r)) for r in str_rows]

    return "\n".join([header_line, sep_line] + content_lines)


def inspect_db(db_path):
    print("=" * 80)
    print(f"DATABASE FILE: {db_path}")
    print(f"FILE SIZE:     {db_path.stat().st_size:,} bytes ({db_path.stat().st_size / (1024*1024):.2f} MB)")
    print("=" * 80)

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # 1. Sources Dimension Table
    print(f"\n[1. Sources Lookup Table]")
    cursor.execute("SELECT id, name, description FROM sources ORDER BY id")
    rows = cursor.fetchall()
    headers = ["ID", "Name", "Description"]
    print(format_table(headers, rows))

    # 2. Daily Call & Point Volume Metrics (daily_stats)
    print(f"\n[2. Daily Call & Point Volume (daily_stats)]")
    cursor.execute("""
        SELECT 
            d.date,
            d.source,
            COALESCE(s.name, 'unknown') as source_name,
            d.count_calls,
            d.count_points
        FROM daily_stats d
        LEFT JOIN sources s ON d.source = s.id
        ORDER BY d.date DESC, d.source
    """)
    rows = cursor.fetchall()
    headers = ["Date", "Source ID", "Source Name", "Total Calls", "Total Points"]
    print(format_table(headers, rows))

    # 3. Spatial Usage Summary
    cursor.execute("SELECT COUNT(*), SUM(count) FROM usage")
    total_cells, total_hits = cursor.fetchone()
    print(f"\n[3. Spatial Usage Summary (usage)]")
    print(f"  * Total Aggregated ~1km Cells: {total_cells:,}")
    print(f"  * Total Spatial Hits:           {total_hits:,}")

    # 4. Breakdown by Direction & Source Name
    print(f"\n[4. Spatial Breakdown by Direction and Source]")
    cursor.execute("""
        SELECT 
            u.is_st70_to_etrs89,
            CASE u.is_st70_to_etrs89 
                WHEN 1 THEN 'Stereo70 -> ETRS89' 
                ELSE 'ETRS89 -> Stereo70' 
            END as direction_desc,
            u.source as source_id,
            COALESCE(s.name, 'unknown') as source_name,
            u.is_3d,
            COUNT(*) as unique_cells,
            SUM(u.count) as total_hits
        FROM usage u
        LEFT JOIN sources s ON u.source = s.id
        GROUP BY u.is_st70_to_etrs89, u.source, u.is_3d
        ORDER BY u.is_st70_to_etrs89 DESC, u.source
    """)
    rows = cursor.fetchall()
    headers = ["is_st70_to_etrs89", "Direction", "Source ID", "Source Name", "3D", "Unique Cells", "Total Hits"]
    print(format_table(headers, rows))

    # 5. Top 10 Geographic Cells by Hit Count
    print(f"\n[5. Top 10 Geographic Cells by Hit Count (~1 km)]")
    cursor.execute("""
        SELECT 
            u.lat_deg, u.lat_grad, u.lon_deg, u.lon_grad,
            printf('%02d.%02d', u.lat_deg, u.lat_grad) as lat_deg_str,
            printf('%02d.%02d', u.lon_deg, u.lon_grad) as lon_deg_str,
            u.is_st70_to_etrs89,
            COALESCE(s.name, 'unknown') as source_name,
            u.count
        FROM usage u
        LEFT JOIN sources s ON u.source = s.id
        ORDER BY u.count DESC
        LIMIT 10
    """)
    rows = cursor.fetchall()
    headers = ["Lat_Deg", "Lat_Grad", "Lon_Deg", "Lon_Grad", "Lat (deg)", "Lon (deg)", "Fwd", "Source", "Count"]
    print(format_table(headers, rows))

    # 6. Sample Rows from usage
    print(f"\n[6. Sample Rows from usage (First 10 records)]")
    cursor.execute("""
        SELECT year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source, count
        FROM usage
        LIMIT 10
    """)
    rows = cursor.fetchall()
    headers = ["Year", "Month", "Lat_Deg", "Lat_Grad", "Lon_Deg", "Lon_Grad", "Fwd", "3D", "Source", "Count"]
    print(format_table(headers, rows))

    conn.close()
    print("\n" + "=" * 80)


def main():
    db_path = LOGS_DIR / "usage.db"
    if not db_path.exists():
        print(f"No usage.db found in {LOGS_DIR}.")
        sys.exit(1)
    inspect_db(db_path)


if __name__ == "__main__":
    main()
