"""Pure Python feedback dispatcher with SMTP email delivery and local JSONL logging fallback.

Follows project invariants: pure Python standard library only (smtplib, email, json, pathlib).
No third-party dependencies required.
"""
import os
import json
import smtplib
from datetime import datetime, timezone
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from pathlib import Path
from typing import Dict, Any, Tuple


def _get_env_config() -> Dict[str, Any]:
    """Extracts SMTP and recipient configuration from environment variables."""
    return {
        "smtp_host": os.environ.get("SMTP_HOST", "").strip(),
        "smtp_port": int(os.environ.get("SMTP_PORT", "587")),
        "smtp_user": os.environ.get("SMTP_USER", "").strip(),
        "smtp_password": os.environ.get("SMTP_PASSWORD", "").strip(),
        "smtp_use_tls": os.environ.get("SMTP_USE_TLS", "true").lower() in ("true", "1", "yes"),
        "recipient_email": os.environ.get("FEEDBACK_RECIPIENT_EMAIL", "").strip(),
        "sender_email": os.environ.get("FEEDBACK_SENDER_EMAIL", "").strip() or os.environ.get("SMTP_USER", "").strip(),
    }


def _log_feedback_locally(payload: Dict[str, Any], client_ip: str) -> None:
    """Appends feedback entry to a local JSONL log file under logs/feedback.jsonl."""
    log_dir = Path("./logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "feedback.jsonl"

    record = {
        "received_at": datetime.now(timezone.utc).isoformat(),
        "client_ip": client_ip,
        **payload,
    }

    with open(log_file, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def dispatch_feedback(payload: Dict[str, Any], client_ip: str) -> Tuple[bool, str]:
    """Dispatches user feedback.

    Attempts SMTP email dispatch if credentials are provided in the environment.
    Always safely logs to local `logs/feedback.jsonl`.

    Returns:
        (delivered_email: bool, confirmation_message: str)
    """
    # 1. Always record in local log file for auditing / backup
    try:
        _log_feedback_locally(payload, client_ip)
    except Exception as e:
        # Don't let disk logging error halt overall pipeline
        print(f"[Feedback] Warning: Failed to write to local feedback log: {e}")

    cfg = _get_env_config()

    # If SMTP is unconfigured, return success with logged status
    if not cfg["smtp_host"] or not cfg["recipient_email"]:
        return (
            False,
            "Feedback recorded in server logs. (SMTP unconfigured on server).",
        )

    # 2. Attempt SMTP transmission
    try:
        report_type = payload.get("report_type", "general").upper()
        subject = f"[PyTransDatRO Feedback] {report_type} Report"

        body_lines = [
            f"PyTransDatRO User Feedback",
            f"----------------------------------------",
            f"Report Type:      {report_type}",
            f"Timestamp (UTC):  {datetime.now(timezone.utc).isoformat()}",
            f"Client Timestamp: {payload.get('client_timestamp', 'N/A')}",
            f"App Version:      {payload.get('app_version', 'N/A')}",
            f"Active Grid:      {payload.get('active_grid', 'N/A')}",
            f"Client IP:        {client_ip}",
            f"----------------------------------------\n",
            "Description & Details:",
            payload.get("description", "").strip(),
            f"\n----------------------------------------",
        ]
        body = "\n".join(body_lines)

        msg = MIMEMultipart()
        msg["From"] = cfg["sender_email"] or "no-reply@pytransdatro.local"
        msg["To"] = cfg["recipient_email"]
        msg["Subject"] = subject
        msg.attach(MIMEText(body, "plain", "utf-8"))

        server = smtplib.SMTP(cfg["smtp_host"], cfg["smtp_port"], timeout=10)
        try:
            if cfg["smtp_use_tls"]:
                server.starttls()
            if cfg["smtp_user"] and cfg["smtp_password"]:
                server.login(cfg["smtp_user"], cfg["smtp_password"])
            server.sendmail(msg["From"], [cfg["recipient_email"]], msg.as_string())
        finally:
            server.quit()

        return (True, "Feedback delivered successfully via email.")

    except Exception as e:
        print(f"[Feedback] Warning: SMTP delivery failed ({e}); recorded in local log.")
        return (False, f"Feedback saved to server log. (SMTP notification failed: {e})")
