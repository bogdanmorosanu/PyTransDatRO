import json
import logging
import smtplib
from datetime import datetime, timezone
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from api.config import Settings, get_settings

logger = logging.getLogger("api.services.mailer")


def _get_smtp_config(settings: Optional[Settings] = None) -> Dict[str, Any]:
    """Extracts SMTP and recipient configuration from Settings."""
    cfg = settings or get_settings()
    return {
        "smtp_host": cfg.smtp_host.strip(),
        "smtp_port": cfg.smtp_port,
        "smtp_user": cfg.smtp_user.strip(),
        "smtp_password": cfg.smtp_password.strip(),
        "smtp_use_tls": cfg.smtp_use_tls,
        "recipient_email": cfg.feedback_recipient_email.strip(),
        "sender_email": cfg.feedback_sender_email.strip() or cfg.smtp_user.strip(),
    }


def _log_feedback_locally(
    payload: Dict[str, Any], client_ip: str, settings: Optional[Settings] = None
) -> Path:
    """Appends feedback entry to a local JSONL log file under logs/feedback.jsonl."""
    cfg = settings or get_settings()
    if cfg.feedback_log_path and cfg.feedback_log_path != "./logs/feedback.jsonl":
        log_file = Path(cfg.feedback_log_path)
    else:
        log_file = Path(cfg.log_dir) / "feedback.jsonl"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    record = {
        "received_at": datetime.now(timezone.utc).isoformat(),
        "client_ip": client_ip,
        **payload,
    }

    with open(log_file, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return log_file


def dispatch_feedback(
    payload: Dict[str, Any], client_ip: str, settings: Optional[Settings] = None
) -> Tuple[bool, str]:
    """Dispatches user feedback.

    Attempts SMTP email dispatch if credentials are provided in settings.
    Always safely logs to local `logs/feedback.jsonl`.

    Returns:
        (delivered_email: bool, confirmation_message: str)
    """
    cfg = settings or get_settings()

    # 1. Always record in local log file for auditing / backup
    saved_log_path = str(
        Path(cfg.feedback_log_path)
        if cfg.feedback_log_path != "./logs/feedback.jsonl"
        else Path(cfg.log_dir) / "feedback.jsonl"
    )
    try:
        log_file = _log_feedback_locally(payload, client_ip, settings=cfg)
        saved_log_path = str(log_file)
    except (OSError, TypeError, ValueError) as e:
        # Don't let disk logging error halt overall pipeline
        logger.error("Failed to write to local feedback log (%s): %s", saved_log_path, e)

    smtp_cfg = _get_smtp_config(cfg)

    # If SMTP is unconfigured, return success with logged status
    if not smtp_cfg["smtp_host"] or not smtp_cfg["recipient_email"]:
        return (
            False,
            "Feedback recorded in server logs. (SMTP unconfigured on server).",
        )

    # 2. Attempt SMTP transmission
    try:
        report_type = payload.get("report_type", "general").upper()
        subject = f"[PyTransDatRO Feedback] {report_type} Report"

        body_lines = [
            "PyTransDatRO User Feedback",
            "----------------------------------------",
            f"Report Type:      {report_type}",
            f"Timestamp (UTC):  {datetime.now(timezone.utc).isoformat()}",
            f"Client Timestamp: {payload.get('client_timestamp', 'N/A')}",
            f"App Version:      {payload.get('app_version', 'N/A')}",
            f"Active Grid:      {payload.get('active_grid', 'N/A')}",
            f"Client IP:        {client_ip}",
            "----------------------------------------\n",
            "Description & Details:",
            payload.get("description", "").strip(),
            "\n----------------------------------------",
        ]
        body = "\n".join(body_lines)

        msg = MIMEMultipart()
        msg["From"] = smtp_cfg["sender_email"] or "no-reply@pytransdatro.local"
        msg["To"] = smtp_cfg["recipient_email"]
        msg["Subject"] = subject
        msg.attach(MIMEText(body, "plain", "utf-8"))

        server = smtplib.SMTP(smtp_cfg["smtp_host"], smtp_cfg["smtp_port"], timeout=10)
        try:
            if smtp_cfg["smtp_use_tls"]:
                server.starttls()
            if smtp_cfg["smtp_user"] and smtp_cfg["smtp_password"]:
                server.login(smtp_cfg["smtp_user"], smtp_cfg["smtp_password"])
            server.sendmail(msg["From"], [smtp_cfg["recipient_email"]], msg.as_string())
        finally:
            server.quit()

        return (True, "Feedback delivered successfully via email.")

    except (smtplib.SMTPException, OSError, TimeoutError) as e:
        logger.warning(
            "SMTP feedback delivery failed (%s: %s). Feedback safely preserved in %s.",
            type(e).__name__,
            e,
            saved_log_path,
        )
        return (False, f"Feedback saved to server log. (SMTP notification failed: {e})")
