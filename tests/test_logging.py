"""Tests for application rotating file logging infrastructure and telemetry log integration."""

import logging
from unittest.mock import MagicMock

from api.config import Settings, configure_logging
from pytransdatro.telemetry import SqliteUsageLogger


def test_rotating_file_handler_creation_and_formatting(tmp_path):
    """Verifies configure_logging creates the log file and formats logs with timestamp and level."""
    log_file = tmp_path / "logs" / "app.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
    )

    configure_logging(settings)

    api_logger = logging.getLogger("api.test")
    core_logger = logging.getLogger("pytransdatro.test")

    api_logger.info("API verification message 12345")
    core_logger.warning("Core warning message 67890")

    # Flush handlers
    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in api_logger.handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")

    assert "API verification message 12345" in content
    assert "[INFO] api.test" in content
    assert "Core warning message 67890" in content
    assert "[WARNING] pytransdatro.test" in content


def test_log_file_rotation(tmp_path):
    """Verifies log file rotates when exceeding max_bytes and retains backup_count files."""
    log_file = tmp_path / "rotate" / "app.log"
    max_bytes = 1024  # 1 KB threshold for testing rotation
    backup_count = 3

    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        app_log_max_bytes=max_bytes,
        app_log_backup_count=backup_count,
        log_level="INFO",
    )

    configure_logging(settings)
    test_logger = logging.getLogger("api.rotation")

    # Write enough lines to exceed 1 KB multiple times (~4 KB total)
    msg = "X" * 120
    for i in range(40):
        test_logger.info("Message %03d: %s", i, msg)

    for handler in test_logger.handlers:
        handler.flush()

    assert log_file.exists()
    # At least app.log.1 should exist
    backup_1 = tmp_path / "rotate" / "app.log.1"
    assert backup_1.exists()

    # Backups should not exceed backup_count (i.e. app.log.4 must NOT exist)
    backup_4 = tmp_path / "rotate" / f"app.log.{backup_count + 1}"
    assert not backup_4.exists()


def test_telemetry_queue_saturation_logged_to_file(tmp_path):
    """Verifies that telemetry queue saturation logs a warning into the application log file."""
    log_file = tmp_path / "app_telem.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="DEBUG",
    )
    configure_logging(settings)

    # Logger with tiny queue size
    telem_dir = tmp_path / "telem_db"
    telem = SqliteUsageLogger(log_dir=str(telem_dir), max_queue_size=2)

    try:
        # Pause worker thread processing briefly or flood immediately
        for i in range(20):
            telem.record(1, [45.0 + i * 0.01], [25.0 + i * 0.01], is_3d=0, source=1)

        # Flush logging handlers
        for handler in logging.getLogger().handlers:
            handler.flush()

        diag = telem.get_diagnostics()
        assert diag["dropped_events_count"] > 0

        assert log_file.exists()
        content = log_file.read_text(encoding="utf-8")
        assert "Telemetry queue saturated" in content
        assert "Dropped" in content
    finally:
        telem.close()


def test_telemetry_db_write_failure_fallback_logging(tmp_path):
    """Verifies that persistent SQLite write errors fallback log unwritten batch payloads to app.log."""
    log_file = tmp_path / "app_fallback.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="DEBUG",
    )
    configure_logging(settings)

    telem_dir = tmp_path / "telem_broken"
    telem = SqliteUsageLogger(log_dir=str(telem_dir), max_queue_size=100)

    try:
        # Mock connection inside _process_batch to raise sqlite3.OperationalError
        mock_conn = MagicMock()
        import sqlite3

        mock_conn.executemany.side_effect = sqlite3.OperationalError(
            "database disk image is malformed"
        )

        # Directly invoke _process_batch with sample items
        import datetime

        batch = [(1, [0.785398], [0.436332], 0, 1)]
        telem._process_batch(mock_conn, batch, datetime.date.today())

        for handler in logging.getLogger().handlers:
            handler.flush()

        assert log_file.exists()
        content = log_file.read_text(encoding="utf-8")
        assert "SQLite write failed after 3 attempts" in content
        assert "Logging unwritten batch payload" in content
        assert "Unwritten telemetry computation payload" in content
        assert "Unwritten telemetry spatial payload" in content

        diag = telem.get_diagnostics()
        assert "database disk image is malformed" in diag["last_error"]
    finally:
        telem.close()
