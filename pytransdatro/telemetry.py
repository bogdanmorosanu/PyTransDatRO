import atexit
import calendar
import datetime
import logging
import math
import queue
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger("pytransdatro.telemetry")


class SqliteUsageLogger:
    def __init__(
        self,
        log_dir: str = "./logs",
        max_queue_size: int = 10000,
        fallback_path: Optional[str] = None,
    ):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.fallback_path = (
            Path(fallback_path) if fallback_path else self.log_dir / "telemetry_fallback.jsonl"
        )
        self._queue = queue.Queue(maxsize=max_queue_size)

        self._lock = threading.Lock()
        self._dropped_events_count = 0
        self._total_recorded_points = 0
        self._last_error: Optional[str] = None
        self._last_error_time: Optional[float] = None
        self._closed = False

        self._stop_event = threading.Event()
        self._worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._worker_thread.start()

        atexit.register(self.close)

    def record(self, operation, lats, lons, is_3d, source=0):
        if self._closed or self._stop_event.is_set():
            return

        is_scalar = not isinstance(lats, (list, tuple))
        n_pts = 1 if is_scalar else len(lats)

        try:
            self._queue.put_nowait((operation, lats, lons, is_3d, source))
        except queue.Full:
            with self._lock:
                self._dropped_events_count += 1
            logger.warning(
                "Telemetry queue saturated (size=%d/%d). Dropped %d points for operation %s (source=%s).",
                self._queue.qsize(),
                self._queue.maxsize,
                n_pts,
                operation,
                source,
            )

    def get_diagnostics(self) -> dict:
        """Returns runtime diagnostic metrics of the telemetry collector."""
        with self._lock:
            last_err = self._last_error
            last_err_time = self._last_error_time
            dropped = self._dropped_events_count
            total_points = self._total_recorded_points

        return {
            "is_alive": self._worker_thread.is_alive(),
            "queue_size": self._queue.qsize(),
            "max_queue_size": self._queue.maxsize,
            "dropped_events_count": dropped,
            "total_recorded_points": total_points,
            "last_error": last_err,
            "last_error_time": last_err_time,
            "fallback_path": str(self.fallback_path),
        }

    def _save_fallback_payload(
        self,
        *,
        reason: str,
        comp_params: Optional[list] = None,
        spatial_params: Optional[list] = None,
        raw_batch: Optional[list] = None,
    ) -> None:
        """Appends unwritten telemetry items to a structured JSONL fallback file."""
        import json

        record = {
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "reason": reason,
        }
        if comp_params is not None:
            record["comp_params"] = comp_params
        if spatial_params is not None:
            record["spatial_params"] = spatial_params
        if raw_batch is not None:
            record["raw_batch"] = raw_batch

        try:
            self.fallback_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.fallback_path, "a", encoding="utf-8") as f:
                f.write(json.dumps(record, ensure_ascii=False) + "\n")
        except Exception as file_exc:
            logger.error(
                "Failed to write to telemetry fallback file %s: %s",
                self.fallback_path,
                file_exc,
            )

    def _worker_loop(self):
        conn = None
        try:
            conn = self._get_db_connection()
        except Exception as exc:
            with self._lock:
                self._last_error = str(exc)
                self._last_error_time = time.time()
            logger.error("Failed to initialize SQLite telemetry database: %s", exc)

        while not self._stop_event.is_set() or not self._queue.empty():
            try:
                # Wait for items, timeout allows checking the stop event
                item = self._queue.get(timeout=0.2)
            except queue.Empty:
                continue

            batch = [item]
            # Drain up to 1000 items from the queue
            while len(batch) < 1000:
                try:
                    batch.append(self._queue.get_nowait())
                except queue.Empty:
                    break

            if conn is None:
                try:
                    conn = self._get_db_connection()
                except Exception as exc:
                    with self._lock:
                        self._last_error = str(exc)
                        self._last_error_time = time.time()
                        self._dropped_events_count += len(batch)
                    logger.error(
                        "SQLite connection unavailable: %s. Saved unwritten queue batch (%d items) to %s.",
                        exc,
                        len(batch),
                        self.fallback_path,
                    )
                    self._save_fallback_payload(
                        reason=f"SQLite connection unavailable: {exc}",
                        raw_batch=batch,
                    )
                    continue

            now = datetime.datetime.now()
            today_date = now.date()

            try:
                self._process_batch(conn, batch, today_date)
            except Exception as exc:
                with self._lock:
                    self._last_error = str(exc)
                    self._last_error_time = time.time()
                    self._dropped_events_count += len(batch)
                logger.error(
                    "Unexpected error in telemetry worker loop: %s. Saved unwritten queue batch (%d items) to %s.",
                    exc,
                    len(batch),
                    self.fallback_path,
                    exc_info=True,
                )
                self._save_fallback_payload(
                    reason=f"Unexpected error in worker loop: {exc}",
                    raw_batch=batch,
                )

        if conn:
            try:
                conn.close()
            except Exception:
                pass

    def _get_db_connection(self):
        self.log_dir.mkdir(parents=True, exist_ok=True)
        db_path = self.log_dir / "usage.db"
        conn = sqlite3.connect(str(db_path), timeout=5.0)
        try:
            conn.execute("PRAGMA journal_mode = WAL;")
        except sqlite3.OperationalError:
            conn.execute("PRAGMA journal_mode = DELETE;")
        conn.execute("PRAGMA synchronous = NORMAL;")

        # 1. Sources dimension lookup table
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

        # 2. Operations dimension lookup table
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
            conn.executemany(
                """
                INSERT INTO telemetry_config (metric_type, window_days, anchor_date, valid_from, valid_to, description)
                VALUES (?, ?, ?, ?, ?, ?);
            """,
                [
                    (
                        "computation",
                        1,
                        "2023-01-01",
                        "2023-01-01",
                        None,
                        "Daily computation usage (calls & points)",
                    ),
                    (
                        "spatial",
                        -1,
                        "2023-01-01",
                        "2023-01-01",
                        None,
                        "Monthly approximated spatial cell usage",
                    ),
                ],
            )

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
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_comp_usage_period ON computation_usage(period_start);"
        )

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
                count_calls INTEGER NOT NULL DEFAULT 1,
                count INTEGER NOT NULL DEFAULT 1,
                PRIMARY KEY (period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
            );
        """)
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_spatial_usage_period ON spatial_usage(period_start);"
        )

        # Auto-migration: ensure count_calls column exists on spatial_usage
        cursor.execute("PRAGMA table_info(spatial_usage)")
        spatial_cols = {col[1] for col in cursor.fetchall()}
        if "count_calls" not in spatial_cols:
            try:
                conn.execute(
                    "ALTER TABLE spatial_usage ADD COLUMN count_calls INTEGER NOT NULL DEFAULT 1;"
                )
            except sqlite3.OperationalError:
                pass

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
                    migrated_rows.append(
                        (p_start, p_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src, cnt)
                    )

                cursor.executemany(
                    """
                    INSERT INTO spatial_usage
                        (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count_calls, count)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 1, ?)
                    ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
                    DO UPDATE SET
                        count_calls = count_calls + excluded.count_calls,
                        count = count + excluded.count;
                """,
                    migrated_rows,
                )

            cursor.execute("DROP TABLE usage;")

        # Migrate source IDs if sources table has old assignment (1=web_map, 2=rest_api)
        cursor.execute("SELECT name FROM sources WHERE id = 1")
        row_1 = cursor.fetchone()
        if row_1 and row_1[0] == "web_map":
            # Temporary renumbering to avoid primary key conflict
            cursor.execute("UPDATE sources SET id = -1 WHERE id = 1;")
            cursor.execute("UPDATE sources SET id = -2 WHERE id = 2;")
            cursor.execute(
                "UPDATE sources SET id = 2, name = 'web_map', description = 'Web Map Application' WHERE id = -1;"
            )
            cursor.execute(
                "UPDATE sources SET id = 1, name = 'rest_api', description = 'REST API Service' WHERE id = -2;"
            )

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
        cursor.execute(
            """
            SELECT metric_type, window_days, anchor_date
            FROM telemetry_config
            WHERE valid_from <= ? AND (valid_to IS NULL OR valid_to >= ?)
            ORDER BY id DESC
        """,
            (target_date_str, target_date_str),
        )

        configs = {}
        for metric_type, window_days, anchor_date in cursor.fetchall():
            if metric_type not in configs:
                configs[metric_type] = (window_days, anchor_date)

        # Fallback defaults if not matched
        if "computation" not in configs:
            configs["computation"] = (1, "2023-01-01")
        if "spatial" not in configs:
            configs["spatial"] = (-1, "2023-01-01")

        return configs

    def _calculate_period_window(self, date_obj, window_days, anchor_date_str):
        if window_days == -2:
            # Calendar Year
            y = date_obj.year
            return f"{y:04d}-01-01", f"{y:04d}-12-31"

        if window_days == -1 or window_days == 30:
            # Calendar Month (-1 sentinel or legacy 30)
            y, m = date_obj.year, date_obj.month
            last_day = calendar.monthrange(y, m)[1]
            return f"{y:04d}-{m:02d}-01", f"{y:04d}-{m:02d}-{last_day:02d}"

        if window_days <= 1:
            date_str = date_obj.strftime("%Y-%m-%d")
            return date_str, date_str

        # Mathematical N-day steps from anchor_date
        anchor = datetime.date.fromisoformat(anchor_date_str)
        delta_days = (date_obj - anchor).days
        bucket_index = delta_days // window_days
        period_start = anchor + datetime.timedelta(days=bucket_index * window_days)
        period_end = period_start + datetime.timedelta(days=window_days - 1)
        return period_start.strftime("%Y-%m-%d"), period_end.strftime("%Y-%m-%d")

    @classmethod
    def update_window_config(
        cls,
        conn: sqlite3.Connection,
        metric_type: str,
        window_days: int,
        valid_from: Optional[str] = None,
        anchor_date: Optional[str] = None,
        description: Optional[str] = None,
    ) -> dict:
        """Transitions telemetry window configuration cleanly starting from valid_from.

        To prevent data loss and period overlap:
        - If valid_from is omitted, defaults to tomorrow.
        - Clamps any active open period's period_end in `{metric_type}_usage` to yesterday (valid_from - 1 day).
        - Closes the previous active telemetry_config row with valid_to = yesterday.
        - Inserts the new configuration starting on valid_from.
        """
        today = datetime.date.today()
        if valid_from is None:
            valid_from_date = today + datetime.timedelta(days=1)
            valid_from_str = valid_from_date.strftime("%Y-%m-%d")
        else:
            valid_from_date = datetime.date.fromisoformat(valid_from)
            valid_from_str = valid_from

        cutoff_date = valid_from_date - datetime.timedelta(days=1)
        cutoff_str = cutoff_date.strftime("%Y-%m-%d")

        if anchor_date is None:
            anchor_date = valid_from_str

        table_name = f"{metric_type}_usage"

        with conn:
            # 1. Clamp open periods ending after cutoff_str that started on or before cutoff_str
            conn.execute(
                f"""
                UPDATE {table_name}
                SET period_end = ?
                WHERE period_end > ? AND period_start <= ?;
                """,
                (cutoff_str, cutoff_str, cutoff_str),
            )

            # 2. Close previous active config(s)
            conn.execute(
                """
                UPDATE telemetry_config
                SET valid_to = ?
                WHERE metric_type = ? AND (valid_to IS NULL OR valid_to >= ?);
                """,
                (cutoff_str, metric_type, cutoff_str),
            )

            # 3. Insert new config
            desc = (
                description
                or f"{metric_type} interval {window_days} days starting {valid_from_str}"
            )
            cursor = conn.execute(
                """
                INSERT INTO telemetry_config (metric_type, window_days, anchor_date, valid_from, valid_to, description)
                VALUES (?, ?, ?, ?, NULL, ?);
                """,
                (metric_type, window_days, anchor_date, valid_from_str, desc),
            )
            new_id = cursor.lastrowid

        return {
            "id": new_id,
            "metric_type": metric_type,
            "window_days": window_days,
            "anchor_date": anchor_date,
            "valid_from": valid_from_str,
            "cutoff_date": cutoff_str,
        }

    def update_config(
        self,
        metric_type: str,
        window_days: int,
        valid_from: Optional[str] = None,
        anchor_date: Optional[str] = None,
        description: Optional[str] = None,
    ) -> dict:
        """Instance helper to update window config via the logger's database."""
        conn = self._get_db_connection()
        try:
            return self.update_window_config(
                conn, metric_type, window_days, valid_from, anchor_date, description
            )
        finally:
            conn.close()

    def _process_batch(self, conn, batch, today_date):
        today_str = today_date.strftime("%Y-%m-%d")
        configs = self._get_active_config(conn, today_str)

        comp_window_days, comp_anchor = configs.get("computation", (1, "2023-01-01"))
        comp_start, comp_end = self._calculate_period_window(
            today_date, comp_window_days, comp_anchor
        )

        spatial_window_days, spatial_anchor = configs.get("spatial", (-1, "2023-01-01"))
        spatial_start, spatial_end = self._calculate_period_window(
            today_date, spatial_window_days, spatial_anchor
        )

        computation_stats = {}  # (comp_start, comp_end, op, is_3d, src) -> [calls, points]
        spatial_stats = {}  # (spatial_start, spatial_end, lat_d, lat_g, lon_d, lon_g, op, is_3d, src) -> [calls, points]

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
            lat_iter = lats if not is_scalar else (lats,)
            lon_iter = lons if not is_scalar else (lons,)

            valid_spatial_cells = []
            for lat, lon in zip(lat_iter, lon_iter):
                try:
                    lat_f = float(lat)
                    lon_f = float(lon)
                    if not (math.isfinite(lat_f) and math.isfinite(lon_f)):
                        continue
                except (ValueError, TypeError):
                    continue

                total_lat = int(round(math.degrees(lat_f) * 100))
                lat_deg = total_lat // 100
                lat_grad = total_lat % 100

                total_lon = int(round(math.degrees(lon_f) * 100))
                lon_deg = total_lon // 100
                lon_grad = total_lon % 100

                spatial_key = (
                    spatial_start,
                    spatial_end,
                    lat_deg,
                    lat_grad,
                    lon_deg,
                    lon_grad,
                    op_int,
                    is_3d_int,
                    source_int,
                )
                valid_spatial_cells.append(spatial_key)

            if not valid_spatial_cells:
                continue

            # Accumulate computation usage metrics
            comp_key = (comp_start, comp_end, op_int, is_3d_int, source_int)
            if comp_key not in computation_stats:
                computation_stats[comp_key] = [0, 0]
            computation_stats[comp_key][0] += 1
            computation_stats[comp_key][1] += len(valid_spatial_cells)

            # Accumulate spatial usage: each distinct cell touched in this call gets 1 call, plus point counts
            distinct_cells = set(valid_spatial_cells)
            for sk in distinct_cells:
                if sk not in spatial_stats:
                    spatial_stats[sk] = [0, 0]  # [calls, points]
                spatial_stats[sk][0] += 1

            for sk in valid_spatial_cells:
                spatial_stats[sk][1] += 1

        comp_params = [
            (k[0], k[1], k[2], k[3], k[4], v[0], v[1]) for k, v in computation_stats.items()
        ]

        spatial_params = [
            (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], v[0], v[1])
            for k, v in spatial_stats.items()
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
        INSERT INTO spatial_usage (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count_calls, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(period_start, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source)
        DO UPDATE SET
            count_calls = count_calls + excluded.count_calls,
            count = count + excluded.count;
        """

        batch_total_points = sum(v[1] for v in computation_stats.values())
        if not comp_params and not spatial_params:
            return

        max_retries = 3
        for attempt in range(1, max_retries + 1):
            try:
                with conn:
                    conn.executemany(comp_query, comp_params)
                    conn.executemany(spatial_query, spatial_params)
                with self._lock:
                    self._last_error = None
                    self._total_recorded_points += batch_total_points
                return
            except sqlite3.OperationalError as exc:
                if attempt < max_retries:
                    time.sleep(0.05 * attempt)
                else:
                    with self._lock:
                        self._last_error = str(exc)
                        self._last_error_time = time.time()
                    logger.error(
                        "SQLite write failed after %d attempts: %s. Saved unwritten batch (%d comp, %d spatial items) to %s.",
                        max_retries,
                        exc,
                        len(comp_params),
                        len(spatial_params),
                        self.fallback_path,
                    )
                    self._save_fallback_payload(
                        reason=f"SQLite write failed after {max_retries} attempts: {exc}",
                        comp_params=comp_params,
                        spatial_params=spatial_params,
                    )
            except Exception as exc:
                with self._lock:
                    self._last_error = str(exc)
                    self._last_error_time = time.time()
                logger.error(
                    "Unexpected error during telemetry batch write: %s. Saved unwritten batch (%d comp, %d spatial items) to %s.",
                    exc,
                    len(comp_params),
                    len(spatial_params),
                    self.fallback_path,
                )
                self._save_fallback_payload(
                    reason=f"Unexpected error during batch write: {exc}",
                    comp_params=comp_params,
                    spatial_params=spatial_params,
                )
                return

    def close(self, timeout: float = 10.0):
        self._closed = True
        self._stop_event.set()
        if self._worker_thread.is_alive():
            self._worker_thread.join(timeout=timeout)


class BufferedTelemetryCollector:
    """Buffers individual coordinate points and flushes them to the underlying

    telemetry logger in chunks (default 500 points) to prevent queue saturation.
    """

    def __init__(self, target_collector: Any, chunk_size: int = 500) -> None:
        self.target = target_collector
        self.chunk_size = chunk_size
        self._buffers: dict[tuple[int, int, int], tuple[list[float], list[float]]] = {}
        self._lock = threading.Lock()

    def record(self, operation: Any, lats: Any, lons: Any, is_3d: Any = 0, source: int = 0) -> None:
        if self.target is None:
            return

        op_int = 1 if operation in (1, True, "st70_to_etrs89") else 2
        is_3d_int = 1 if is_3d in (1, True) else 0
        try:
            source_int = int(source)
        except (ValueError, TypeError):
            source_int = 0

        key = (op_int, is_3d_int, source_int)

        with self._lock:
            if key not in self._buffers:
                self._buffers[key] = ([], [])
            buf_lats, buf_lons = self._buffers[key]

            if isinstance(lats, (list, tuple)):
                buf_lats.extend(lats)
                buf_lons.extend(lons)
            else:
                buf_lats.append(lats)
                buf_lons.append(lons)

            while len(buf_lats) >= self.chunk_size:
                chunk_lats = buf_lats[: self.chunk_size]
                chunk_lons = buf_lons[: self.chunk_size]
                del buf_lats[: self.chunk_size]
                del buf_lons[: self.chunk_size]
                if self.target is not None:
                    self.target.record(key[0], chunk_lats, chunk_lons, is_3d=key[1], source=key[2])

    def _flush_key_locked(self, key: tuple[int, int, int]) -> None:
        if self.target is None:
            self._buffers.pop(key, None)
            return

        if key in self._buffers:
            buf_lats, buf_lons = self._buffers[key]
            while buf_lats:
                chunk_lats = buf_lats[: self.chunk_size]
                chunk_lons = buf_lons[: self.chunk_size]
                del buf_lats[: self.chunk_size]
                del buf_lons[: self.chunk_size]
                self.target.record(key[0], chunk_lats, chunk_lons, is_3d=key[1], source=key[2])
            del self._buffers[key]

    def flush_all(self) -> None:
        """Flushes all buffered points to the target collector."""
        with self._lock:
            if self.target is None:
                self._buffers.clear()
                return
            for key in list(self._buffers.keys()):
                self._flush_key_locked(key)

    def get_diagnostics(self) -> dict:
        if hasattr(self.target, "get_diagnostics"):
            return self.target.get_diagnostics()
        return {}
