import sqlite3
import threading
import queue
import datetime
import calendar
import math
import atexit
from pathlib import Path


class SqliteUsageLogger:
    def __init__(self, log_dir="./logs", max_queue_size=10000):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self._queue = queue.Queue(maxsize=max_queue_size)

        self._stop_event = threading.Event()
        self._worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._worker_thread.start()

        atexit.register(self.close)

    def record(self, operation, lats, lons, is_3d, source=0):
        try:
            self._queue.put_nowait((operation, lats, lons, is_3d, source))
        except queue.Full:
            pass

    def _worker_loop(self):
        conn = self._get_db_connection()

        while not self._stop_event.is_set() or not self._queue.empty():
            try:
                # Wait for items, timeout allows checking the stop event
                item = self._queue.get(timeout=0.5)
            except queue.Empty:
                continue

            batch = [item]
            # Drain up to 1000 items from the queue
            while len(batch) < 1000:
                try:
                    batch.append(self._queue.get_nowait())
                except queue.Empty:
                    break

            now = datetime.datetime.now()
            today_date = now.date()

            self._process_batch(conn, batch, today_date)

        if conn:
            conn.close()

    def _get_db_connection(self):
        db_path = self.log_dir / "usage.db"
        conn = sqlite3.connect(str(db_path), timeout=5.0)
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")

        # 1. Sources dimension lookup table
        conn.execute("""
            CREATE TABLE IF NOT EXISTS sources (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            );
        """)
        conn.executemany("""
            INSERT OR IGNORE INTO sources (id, name, description)
            VALUES (?, ?, ?);
        """, [
            (0, "unknown", "Default / Unspecified origin"),
            (1, "rest_api", "REST API Service"),
            (2, "web_map", "Web Map Application"),
            (3, "desktop_cli", "Desktop / CLI App"),
        ])

        # 2. Operations dimension lookup table
        conn.execute("""
            CREATE TABLE IF NOT EXISTS operations (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            );
        """)
        conn.executemany("""
            INSERT OR IGNORE INTO operations (id, name, description)
            VALUES (?, ?, ?);
        """, [
            (1, "st70_to_etrs89", "Stereo70 to ETRS89"),
            (2, "etrs89_to_st70", "ETRS89 to Stereo70"),
        ])

        # 3. Telemetry window configuration registry table
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
        # Initialize default configs if table is empty
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM telemetry_config")
        if cursor.fetchone()[0] == 0:
            conn.executemany("""
                INSERT INTO telemetry_config (metric_type, window_days, anchor_date, valid_from, valid_to, description)
                VALUES (?, ?, ?, ?, ?, ?);
            """, [
                ("computation", 1, "2023-01-01", "2023-01-01", None, "Daily computation usage (calls & points)"),
                ("spatial", 30, "2023-01-01", "2023-01-01", None, "Monthly approximated spatial cell usage"),
            ])

        # 4. Computation volume metrics table
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
        conn.execute("CREATE INDEX IF NOT EXISTS idx_comp_usage_period ON computation_usage(period_start);")

        # 5. Spatial usage aggregation table (~1 km geographic cells)
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
                count INTEGER NOT NULL DEFAULT 1,
                PRIMARY KEY (period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
            );
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_spatial_usage_period ON spatial_usage(period_start);")

        # Automatic migration from legacy schema (daily_stats -> computation_usage, usage -> spatial_usage)
        self._migrate_legacy_tables_if_needed(conn)

        conn.commit()
        return conn

    def _migrate_legacy_tables_if_needed(self, conn):
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        existing_tables = {row[0] for row in cursor.fetchall()}

        # Migrate daily_stats -> computation_usage
        if "daily_stats" in existing_tables:
            cursor.execute("""
                INSERT OR IGNORE INTO computation_usage 
                    (period_start, period_end, operation, is_3d, source, count_calls, count_points)
                SELECT 
                    date, 
                    date, 
                    1, 
                    0, 
                    source, 
                    count_calls, 
                    count_points
                FROM daily_stats;
            """)
            cursor.execute("DROP TABLE daily_stats;")

        # Migrate usage -> spatial_usage
        if "usage" in existing_tables:
            # We migrate usage rows: mapping year, month to period_start (YYYY-MM-01) and period_end (YYYY-MM-last_day)
            # and mapping legacy is_st70_to_etrs89 (1 -> 1, 0 -> 2)
            cursor.execute("""
                SELECT year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source, count
                FROM usage
            """)
            legacy_rows = cursor.fetchall()
            if legacy_rows:
                migrated_rows = []
                for row in legacy_rows:
                    y, m, lat_d, lat_g, lon_d, lon_g, fwd, is_3d, src, cnt = row
                    op = 1 if fwd in (1, True) else 2
                    last_day = calendar.monthrange(y, m)[1]
                    p_start = f"{y:04d}-{m:02d}-01"
                    p_end = f"{y:04d}-{m:02d}-{last_day:02d}"
                    migrated_rows.append((p_start, p_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src, cnt))

                cursor.executemany("""
                    INSERT INTO spatial_usage 
                        (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
                    DO UPDATE SET count = count + excluded.count;
                """, migrated_rows)

            cursor.execute("DROP TABLE usage;")

        # Migrate source IDs if sources table has old assignment (1=web_map, 2=rest_api)
        cursor.execute("SELECT name FROM sources WHERE id = 1")
        row_1 = cursor.fetchone()
        if row_1 and row_1[0] == "web_map":
            # Temporary renumbering to avoid primary key conflict
            cursor.execute("UPDATE sources SET id = -1 WHERE id = 1;")
            cursor.execute("UPDATE sources SET id = -2 WHERE id = 2;")
            cursor.execute("UPDATE sources SET id = 2, name = 'web_map', description = 'Web Map Application' WHERE id = -1;")
            cursor.execute("UPDATE sources SET id = 1, name = 'rest_api', description = 'REST API Service' WHERE id = -2;")

            # Swap in computation_usage
            cursor.execute("UPDATE computation_usage SET source = -1 WHERE source = 1;")
            cursor.execute("UPDATE computation_usage SET source = -2 WHERE source = 2;")
            cursor.execute("UPDATE computation_usage SET source = 2 WHERE source = -1;")
            cursor.execute("UPDATE computation_usage SET source = 1 WHERE source = -2;")

            # Swap in spatial_usage
            cursor.execute("UPDATE spatial_usage SET source = -1 WHERE source = 1;")
            cursor.execute("UPDATE spatial_usage SET source = -2 WHERE source = 2;")
            cursor.execute("UPDATE spatial_usage SET source = 2 WHERE source = -1;")
            cursor.execute("UPDATE spatial_usage SET source = 1 WHERE source = -2;")

    def _get_active_config(self, conn, target_date_str):
        cursor = conn.cursor()
        cursor.execute("""
            SELECT metric_type, window_days, anchor_date
            FROM telemetry_config
            WHERE valid_from <= ? AND (valid_to IS NULL OR valid_to >= ?)
            ORDER BY id DESC
        """, (target_date_str, target_date_str))

        configs = {}
        for metric_type, window_days, anchor_date in cursor.fetchall():
            if metric_type not in configs:
                configs[metric_type] = (window_days, anchor_date)

        # Fallback defaults if not matched
        if "computation" not in configs:
            configs["computation"] = (1, "2023-01-01")
        if "spatial" not in configs:
            configs["spatial"] = (30, "2023-01-01")

        return configs

    def _calculate_period_window(self, date_obj, window_days, anchor_date_str):
        if window_days <= 1:
            date_str = date_obj.strftime("%Y-%m-%d")
            return date_str, date_str

        # If window_days is ~1 month (>= 28 days), align to calendar month boundaries
        if window_days >= 28:
            y, m = date_obj.year, date_obj.month
            last_day = calendar.monthrange(y, m)[1]
            return f"{y:04d}-{m:02d}-01", f"{y:04d}-{m:02d}-{last_day:02d}"

        # Otherwise align based on window_days starting from anchor_date
        anchor = datetime.date.fromisoformat(anchor_date_str)
        delta_days = (date_obj - anchor).days
        bucket_index = delta_days // window_days
        period_start = anchor + datetime.timedelta(days=bucket_index * window_days)
        period_end = period_start + datetime.timedelta(days=window_days - 1)
        return period_start.strftime("%Y-%m-%d"), period_end.strftime("%Y-%m-%d")

    def _process_batch(self, conn, batch, today_date):
        today_str = today_date.strftime("%Y-%m-%d")
        configs = self._get_active_config(conn, today_str)

        comp_window_days, comp_anchor = configs.get("computation", (1, "2023-01-01"))
        comp_start, comp_end = self._calculate_period_window(today_date, comp_window_days, comp_anchor)

        spatial_window_days, spatial_anchor = configs.get("spatial", (30, "2023-01-01"))
        spatial_start, spatial_end = self._calculate_period_window(today_date, spatial_window_days, spatial_anchor)

        computation_stats = {}  # (comp_start, comp_end, op, is_3d, src) -> [calls, points]
        spatial_counts = {}     # (spatial_start, spatial_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src) -> count

        for operation, lats, lons, is_3d, source in batch:
            try:
                source_int = int(source)
            except (ValueError, TypeError):
                source_int = 0

            # Normalize operation: 1 (st70_to_etrs89), 2 (etrs89_to_st70)
            if operation in (1, True, "st70_to_etrs89"):
                op_int = 1
            else:
                op_int = 2

            is_3d_int = 1 if is_3d in (1, True) else 0

            is_scalar = not isinstance(lats, (list, tuple))
            n_pts = 1 if is_scalar else len(lats)
            lat_iter = lats if not is_scalar else (lats,)
            lon_iter = lons if not is_scalar else (lons,)

            # Accumulate computation usage metrics
            comp_key = (comp_start, comp_end, op_int, is_3d_int, source_int)
            if comp_key not in computation_stats:
                computation_stats[comp_key] = [0, 0]
            computation_stats[comp_key][0] += 1
            computation_stats[comp_key][1] += n_pts

            # Accumulate spatial ~1km cells
            for lat, lon in zip(lat_iter, lon_iter):
                total_lat = int(round(math.degrees(lat) * 100))
                lat_deg = total_lat // 100
                lat_grad = total_lat % 100

                total_lon = int(round(math.degrees(lon) * 100))
                lon_deg = total_lon // 100
                lon_grad = total_lon % 100

                spatial_key = (
                    spatial_start, spatial_end,
                    lat_deg, lat_grad, lon_deg, lon_grad,
                    op_int, is_3d_int, source_int
                )
                spatial_counts[spatial_key] = spatial_counts.get(spatial_key, 0) + 1

        comp_params = [
            (k[0], k[1], k[2], k[3], k[4], v[0], v[1])
            for k, v in computation_stats.items()
        ]

        spatial_params = [
            (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], count)
            for k, count in spatial_counts.items()
        ]

        comp_query = """
        INSERT INTO computation_usage (period_start, period_end, operation, is_3d, source, count_calls, count_points)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(period_start, operation, is_3d, source)
        DO UPDATE SET 
            count_calls = count_calls + excluded.count_calls,
            count_points = count_points + excluded.count_points;
        """

        spatial_query = """
        INSERT INTO spatial_usage (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
        DO UPDATE SET count = count + excluded.count;
        """

        try:
            with conn:
                conn.executemany(comp_query, comp_params)
                conn.executemany(spatial_query, spatial_params)
        except Exception:
            pass  # Failsafe against write errors

    def close(self):
        self._stop_event.set()
        if self._worker_thread.is_alive():
            self._worker_thread.join(timeout=2.0)
