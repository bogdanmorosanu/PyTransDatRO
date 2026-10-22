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
        assert "Saved unwritten batch" in content
        assert "Unwritten telemetry computation payload" not in content
        assert "Unwritten telemetry spatial payload" not in content

        assert telem.fallback_path.exists()
        fallback_lines = telem.fallback_path.read_text(encoding="utf-8").strip().splitlines()
        assert len(fallback_lines) == 1
        import json

        record = json.loads(fallback_lines[0])
        assert "database disk image is malformed" in record["reason"]
        assert len(record["comp_params"]) == 1

        diag = telem.get_diagnostics()
        assert "database disk image is malformed" in diag["last_error"]
    finally:
        telem.close()


def test_feedback_smtp_failure_logging_and_jsonl_persistence(tmp_path, monkeypatch):
    """Verifies SMTP failures log a warning to app.log while preserving payload in feedback.jsonl."""
    import json
    import smtplib

    from api.services.mailer import dispatch_feedback

    log_file = tmp_path / "app_feedback.log"
    feedback_file = tmp_path / "feedback.jsonl"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        feedback_log_path=str(feedback_file),
        log_level="INFO",
        smtp_host="smtp.example.com",
        smtp_port=587,
        smtp_user="user@example.com",
        smtp_password="SuperSecretPassword123!",
        feedback_recipient_email="admin@example.com",
    )
    configure_logging(settings)

    def failing_smtp(*args, **kwargs):
        raise smtplib.SMTPConnectError(421, b"Service not available, closing transmission channel")

    monkeypatch.setattr(smtplib, "SMTP", failing_smtp)

    payload = {
        "report_type": "discrepancy",
        "description": "Coordinate mismatch observed in Constanta county.",
        "client_timestamp": "2026-10-22T21:00:00Z",
        "app_version": "1.0.0",
        "active_grid": "SPG",
    }

    delivered, msg = dispatch_feedback(payload, "192.168.1.100", settings=settings)
    assert delivered is False
    assert "SMTP notification failed" in msg

    for handler in logging.getLogger().handlers:
        handler.flush()

    # 1. Assert log file exists and warning is recorded without leaking password
    assert log_file.exists()
    log_content = log_file.read_text(encoding="utf-8")
    assert "SMTP feedback delivery failed (SMTPConnectError" in log_content
    assert "Feedback safely preserved in" in log_content
    assert "SuperSecretPassword123!" not in log_content

    # 2. Assert feedback.jsonl exists and contains the user submission
    assert feedback_file.exists()
    feedback_lines = feedback_file.read_text(encoding="utf-8").strip().splitlines()
    assert len(feedback_lines) == 1
    feedback_data = json.loads(feedback_lines[0])
    assert feedback_data["report_type"] == "discrepancy"
    assert "Coordinate mismatch" in feedback_data["description"]
    assert feedback_data["client_ip"] == "192.168.1.100"


def test_unhandled_500_exception_logged_with_traceback(tmp_path, monkeypatch):
    """Verifies unhandled 500 exceptions in API routes log complete tracebacks to app.log."""
    from fastapi.testclient import TestClient

    from api.main import create_app
    from pytransdatro.trans_ro import TransRO

    log_file = tmp_path / "app_500.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
    )
    configure_logging(settings)

    def crashing_transform(*args, **kwargs):
        raise RuntimeError("Simulated unexpected database or server crash")

    monkeypatch.setattr(TransRO, "st70_to_etrs89", crashing_transform)

    app = create_app(settings=settings)

    with TestClient(app, raise_server_exceptions=False) as client:
        resp = client.post(
            "/api/v1/transform/point",
            json={"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0]},
        )
        assert resp.status_code == 500
        data = resp.json()
        assert data["code"] == "INTERNAL_SERVER_ERROR"

    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in logging.getLogger("api.main").handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "Unhandled server exception on POST /api/v1/transform/point" in content
    assert "RuntimeError: Simulated unexpected database or server crash" in content
    assert "Traceback (most recent call last):" in content


def test_rate_limit_warning_logged_to_app_log(tmp_path):
    """Verifies that 429 rate limit triggers emit a structured WARNING to app.log."""
    from fastapi.testclient import TestClient

    from api.main import create_app

    log_file = tmp_path / "app_ratelimit.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
        rate_limit_point=2,  # Low limit for fast testing
        rate_limit_window_seconds=60,
    )
    configure_logging(settings)

    app = create_app(settings=settings)

    with TestClient(app) as client:
        payload = {"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0]}
        # Request 1 & 2: allowed
        r1 = client.post("/api/v1/transform/point", json=payload)
        assert r1.status_code == 200
        r2 = client.post("/api/v1/transform/point", json=payload)
        assert r2.status_code == 200

        # Request 3: throttled (429)
        r3 = client.post("/api/v1/transform/point", json=payload)
        assert r3.status_code == 429

    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in logging.getLogger("api.middleware.rate_limit").handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "[WARNING] api.middleware.rate_limit:" in content
    assert "Rate limit exceeded for client" in content
    assert "/api/v1/transform/point" in content


def test_normal_200_ok_request_does_not_log_to_app_log(tmp_path):
    """Verifies that normal 200 OK requests stay completely silent in app.log."""
    from fastapi.testclient import TestClient

    from api.main import create_app

    log_file = tmp_path / "app_silent.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
    )
    configure_logging(settings)

    app = create_app(settings=settings)

    # Empty log file before request
    if log_file.exists():
        log_file.write_text("", encoding="utf-8")

    with TestClient(app) as client:
        resp = client.post(
            "/api/v1/transform/point",
            json={"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0]},
        )
        assert resp.status_code == 200

    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in logging.getLogger("api.main").handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8").strip()
    # Content must only contain startup / shutdown or be empty, never per-request access lines
    assert "POST /api/v1/transform/point" not in content


def test_file_upload_invalid_delimiter_logged_to_app_log(tmp_path):
    """Verifies rejected file uploads due to invalid delimiter log a WARNING to app.log."""
    import io

    from fastapi.testclient import TestClient

    from api.main import create_app

    log_file = tmp_path / "app_file_delim.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
    )
    configure_logging(settings)

    app = create_app(settings=settings)

    with TestClient(app) as client:
        file_content = b"500000.0,300000.0\n"
        resp = client.post(
            "/api/v1/transform/file",
            files={"file": ("points.csv", io.BytesIO(file_content), "text/csv")},
            data={"op": "Stereo70ToETRS89", "delimiter": "pipe"},
        )
        assert resp.status_code == 422

    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in logging.getLogger("api.routes.file").handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "[WARNING] api.routes.file:" in content
    assert "File upload 'points.csv' rejected: invalid delimiter 'pipe'" in content


def test_file_upload_corrupt_line_warning_logged_to_app_log(tmp_path):
    """Verifies that invalid coordinate lines in uploaded files log a WARNING with line number."""
    import io

    from fastapi.testclient import TestClient

    from api.main import create_app

    log_file = tmp_path / "app_file_corrupt.log"
    settings = Settings(
        app_log_enabled=True,
        app_log_path=str(log_file),
        log_level="INFO",
    )
    configure_logging(settings)

    app = create_app(settings=settings)

    with TestClient(app) as client:
        file_content = b"P1,NOT_A_NUMBER,BAD_DATA\n500000.0,300000.0\n"
        resp = client.post(
            "/api/v1/transform/file",
            files={"file": ("corrupt_survey.csv", io.BytesIO(file_content), "text/csv")},
            data={"op": "Stereo70ToETRS89", "delimiter": "comma"},
        )
        assert resp.status_code == 200
        # Consume streaming body to trigger generator
        body = resp.content.decode("utf-8")
        assert "Invalid coordinate data" in body

    for handler in logging.getLogger().handlers:
        handler.flush()
    for handler in logging.getLogger("api.routes.file").handlers:
        handler.flush()

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "[WARNING] api.routes.file:" in content
    assert "File upload 'corrupt_survey.csv' issue at line 1: Invalid coordinate data" in content
