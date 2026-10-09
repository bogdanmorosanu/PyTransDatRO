import sqlite3
import threading
import queue
import datetime
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

    def record(self, is_st70_to_etrs89, lats, lons, is_3d, source=0):
        try:
            self._queue.put_nowait((is_st70_to_etrs89, lats, lons, is_3d, source))
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
            year = now.year
            month = now.month
            date_str = now.strftime("%Y-%m-%d")

            self._process_batch(conn, batch, year, month, date_str)

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
            (1, "web_map", "Web Map Application"),
            (2, "rest_api", "REST API Service"),
            (3, "desktop_cli", "Desktop / CLI App"),
        ])

        # 2. Daily call & point volume metrics table
        conn.execute("""
            CREATE TABLE IF NOT EXISTS daily_stats (
                date TEXT NOT NULL,
                source INTEGER NOT NULL,
                count_calls INTEGER NOT NULL DEFAULT 1,
                count_points INTEGER NOT NULL DEFAULT 1,
                PRIMARY KEY (date, source)
            );
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_daily_stats_date ON daily_stats(date);")

        # 3. Spatial aggregation table (1 km geographic cells)
        conn.execute("""
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
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_usage_year_month ON usage(year, month);")

        conn.commit()
        return conn

    def _process_batch(self, conn, batch, year, month, date_str):
        usage_counts = {}
        daily_stats = {}  # (date_str, source_int) -> [count_calls, count_points]

        for is_st70_to_etrs89, lats, lons, is_3d, source in batch:
            try:
                source_int = int(source)
            except (ValueError, TypeError):
                source_int = 0

            direction_int = 1 if is_st70_to_etrs89 in (1, True, "st70_to_etrs89") else 0

            is_scalar = not isinstance(lats, (list, tuple))
            n_pts = 1 if is_scalar else len(lats)
            lat_iter = lats if not is_scalar else (lats,)
            lon_iter = lons if not is_scalar else (lons,)

            # Accumulate daily call & point metrics
            daily_key = (date_str, source_int)
            if daily_key not in daily_stats:
                daily_stats[daily_key] = [0, 0]
            daily_stats[daily_key][0] += 1
            daily_stats[daily_key][1] += n_pts

            # Accumulate spatial ~1km cells
            for lat, lon in zip(lat_iter, lon_iter):
                total_lat = int(round(math.degrees(lat) * 100))
                lat_deg = total_lat // 100
                lat_grad = total_lat % 100

                total_lon = int(round(math.degrees(lon) * 100))
                lon_deg = total_lon // 100
                lon_grad = total_lon % 100

                key = (year, month, lat_deg, lat_grad, lon_deg, lon_grad, direction_int, is_3d, source_int)
                usage_counts[key] = usage_counts.get(key, 0) + 1

        usage_params = [
            (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], count)
            for k, count in usage_counts.items()
        ]

        daily_params = [
            (k[0], k[1], v[0], v[1])
            for k, v in daily_stats.items()
        ]

        usage_query = """
        INSERT INTO usage (year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source)
        DO UPDATE SET count = count + excluded.count;
        """

        daily_query = """
        INSERT INTO daily_stats (date, source, count_calls, count_points)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(date, source)
        DO UPDATE SET 
            count_calls = count_calls + excluded.count_calls,
            count_points = count_points + excluded.count_points;
        """

        try:
            with conn:
                conn.executemany(daily_query, daily_params)
                conn.executemany(usage_query, usage_params)
        except Exception:
            pass  # Failsafe against write errors

    def close(self):
        self._stop_event.set()
        if self._worker_thread.is_alive():
            self._worker_thread.join(timeout=2.0)
