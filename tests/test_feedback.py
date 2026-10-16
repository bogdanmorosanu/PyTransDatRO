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
