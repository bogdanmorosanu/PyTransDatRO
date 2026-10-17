"""Tests for anonymous, zero-PII feedback submission and mailer/logger services."""

import json
from pathlib import Path

from fastapi.testclient import TestClient

from api.main import app

client = TestClient(app)


def test_submit_feedback_honeypot_misdirection():
    """Honeypot filled by bot should return success silently without doing work."""
    response = client.post(
        "/api/v1/feedback",
        json={
            "report_type": "bug",
            "description": "Spam bot payload testing the form",
            "honeypot": "https://spam.example.com",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"
    assert data["delivered_email"] is False


def test_submit_feedback_validation_min_length():
    """Short descriptions should be rejected by Pydantic validation."""
    response = client.post(
        "/api/v1/feedback",
        json={
            "report_type": "bug",
            "description": "too short",
        },
    )
    assert response.status_code == 422


def test_submit_feedback_successful_local_logging(tmp_path, monkeypatch):
    """Legitimate feedback should be logged to logs/feedback.jsonl when SMTP unconfigured."""
    log_dir = tmp_path / "logs"
    monkeypatch.setattr("api.services.mailer.Path", lambda p: log_dir if p == "./logs" else Path(p))

    response = client.post(
        "/api/v1/feedback",
        json={
            "report_type": "discrepancy",
            "description": "Coordinates at point X differed by 0.005m from official mark.",
            "client_timestamp": "2026-10-14T19:00:00Z",
            "app_version": "v2.1",
            "active_grid": "RO_TransDatRO_ETRS89_Stereo70.spg",
            "honeypot": "",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "logged"
    assert data["delivered_email"] is False
    assert "Thank you" in data["message"]

    log_file = log_dir / "feedback.jsonl"
    assert log_file.exists()
    with open(log_file, "r", encoding="utf-8") as f:
        line = f.readline()
        record = json.loads(line)
        assert record["report_type"] == "discrepancy"
        assert record["active_grid"] == "RO_TransDatRO_ETRS89_Stereo70.spg"
        assert "differed by 0.005m" in record["description"]


def test_dispatch_feedback_with_custom_settings(tmp_path):
    """Verifies dispatch_feedback respects injected Settings."""
    from api.config import Settings
    from api.services.mailer import dispatch_feedback

    custom_log_dir = tmp_path / "custom_logs"
    settings = Settings(log_dir=str(custom_log_dir))

    payload = {
        "report_type": "bug",
        "description": "Custom settings test description.",
    }
    delivered, msg = dispatch_feedback(payload, "127.0.0.1", settings=settings)

    assert delivered is False
    assert "recorded in server logs" in msg
    assert (custom_log_dir / "feedback.jsonl").exists()


def test_submit_feedback_smtp_delivery_success_and_forwarded_ip(tmp_path, monkeypatch):
    """Verifies SMTP TLS/auth email dispatch, X-Forwarded-For IP extraction, and JSONL backup."""
    from api.config import Settings, get_settings

    sent_calls = []

    class DummySMTP:
        def __init__(self, host, port, timeout=10):
            self.host = host
            self.port = port
            self.tls_started = False
            self.logged_in = None
            self.quitted = False

        def starttls(self):
            self.tls_started = True

        def login(self, user, password):
            self.logged_in = (user, password)

        def sendmail(self, from_addr, to_addrs, msg_str):
            sent_calls.append(
                {
                    "host": self.host,
                    "port": self.port,
                    "tls": self.tls_started,
                    "login": self.logged_in,
                    "from": from_addr,
                    "to": to_addrs,
                    "msg": msg_str,
                }
            )

        def quit(self):
            self.quitted = True

    monkeypatch.setattr("api.services.mailer.smtplib.SMTP", DummySMTP)
    log_dir = tmp_path / "smtp_logs"
    smtp_settings = Settings(
        log_dir=str(log_dir),
        smtp_host="smtp.example.org",
        smtp_port=587,
        smtp_user="notifier@example.org",
        smtp_password="secret",
        smtp_use_tls=True,
        feedback_recipient_email="admin@example.org",
        feedback_sender_email="noreply@example.org",
    )
    app.dependency_overrides[get_settings] = lambda: smtp_settings
    try:
        res = client.post(
            "/api/v1/feedback",
            headers={"X-Forwarded-For": "203.0.113.42, 10.0.0.1"},
            json={
                "report_type": "feature",
                "description": "Please support additional metadata export formats in the UI.",
                "app_version": "v1.0.0",
            },
        )
        assert res.status_code == 200
        data = res.json()
        assert data["status"] == "success"
        assert data["delivered_email"] is True
        assert "sent successfully" in data["message"]

        assert len(sent_calls) == 1
        call = sent_calls[0]
        assert call["host"] == "smtp.example.org"
        assert call["port"] == 587
        assert call["tls"] is True
        assert call["login"] == ("notifier@example.org", "secret")
        assert call["from"] == "noreply@example.org"
        assert call["to"] == ["admin@example.org"]

        # Verify local JSONL backup recorded the first X-Forwarded-For IP
        record = json.loads(
            (log_dir / "feedback.jsonl").read_text(encoding="utf-8").splitlines()[0]
        )
        assert record["client_ip"] == "203.0.113.42"
        assert record["report_type"] == "feature"
    finally:
        app.dependency_overrides.pop(get_settings, None)


def test_submit_feedback_smtp_exception_fallback(tmp_path, monkeypatch):
    """Verifies that if SMTP raises an exception, feedback still succeeds with 'logged' status."""
    import smtplib

    from api.config import Settings, get_settings

    class FailingSMTP:
        def __init__(self, host, port, timeout=10):
            raise smtplib.SMTPConnectError(421, "Service not available")

    monkeypatch.setattr("api.services.mailer.smtplib.SMTP", FailingSMTP)
    log_dir = tmp_path / "fallback_logs"
    smtp_settings = Settings(
        log_dir=str(log_dir),
        smtp_host="smtp.broken.local",
        feedback_recipient_email="admin@example.org",
    )
    app.dependency_overrides[get_settings] = lambda: smtp_settings
    try:
        res = client.post(
            "/api/v1/feedback",
            json={
                "report_type": "api_down",
                "description": "Testing resilient local logging when SMTP server is unreachable.",
            },
        )
        assert res.status_code == 200
        data = res.json()
        assert data["status"] == "logged"
        assert data["delivered_email"] is False
        assert (log_dir / "feedback.jsonl").exists()
    finally:
        app.dependency_overrides.pop(get_settings, None)


def test_submit_feedback_validation_max_length_and_report_type():
    """Verifies invalid report_type and >5000 character descriptions return HTTP 422."""
    res_bad_type = client.post(
        "/api/v1/feedback",
        json={"report_type": "unknown_category", "description": "Valid length description string."},
    )
    assert res_bad_type.status_code == 422

    res_too_long = client.post(
        "/api/v1/feedback",
        json={"report_type": "grid_update", "description": "x" * 5001},
    )
    assert res_too_long.status_code == 422
