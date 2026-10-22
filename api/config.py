"""Centralized configuration module for PyTransDatRO Web Service.

Loads configuration parameters from environment variables (prefixed with PYTRANSDAT_)
or local .env file, providing sensible out-of-the-box defaults for local development.
"""

import json
from functools import lru_cache
from pathlib import Path
from typing import List, Union

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

# Repository root; relative paths to repo build artifacts resolve against it, not the CWD.
_REPO_ROOT = Path(__file__).resolve().parent.parent


class Settings(BaseSettings):
    """Application settings backed by pydantic-settings."""

    model_config = SettingsConfigDict(
        env_prefix="PYTRANSDAT_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # Server network & process settings
    host: str = Field(default="0.0.0.0", description="Host interface to bind uvicorn server")
    port: int = Field(default=8000, description="Network port for API service")
    log_level: str = Field(
        default="INFO", description="Standard logging level (DEBUG, INFO, WARNING, ERROR)"
    )

    # Security & CORS settings
    cors_origins: Union[List[str], str] = Field(
        default=["*"],
        description="Allowed CORS origins list or '*' for open access",
    )

    # In-memory sliding-window rate limiting
    rate_limit_window_seconds: float = Field(
        default=60.0,
        description="Time window in seconds for sliding-window rate calculations",
    )
    rate_limit_point: int = Field(
        default=60,
        description="Maximum single-point transformation requests per window per IP",
    )
    rate_limit_batch: int = Field(
        default=10,
        description="Maximum batch/file transformation requests per window per IP",
    )

    # Transformation throughput & safety limits
    max_batch_points: int = Field(
        default=20000,
        description="Maximum coordinate points allowed in a single in-memory batch request",
    )
    max_file_points: int = Field(
        default=100000,
        description="Maximum coordinate points allowed in streaming file transformation",
    )
    max_file_upload_bytes: int = Field(
        default=50 * 1024 * 1024,
        description="Maximum allowable file upload payload size in bytes (default 50 MB)",
    )

    # Telemetry & Logging
    log_dir: str = Field(
        default="./logs",
        description="Directory for SQLite WAL usage database and telemetry logs",
    )
    legacy_db_path: str = Field(
        default="./logs/legacy_stats.db",
        description="Path to SQLite legacy TransdatOnline statistics database",
    )
    feedback_log_path: str = Field(
        default="./logs/feedback.jsonl",
        description="Path to user feedback submissions log file",
    )
    telemetry_fallback_path: str = Field(
        default="./logs/telemetry_fallback.jsonl",
        description="Path to telemetry JSONL emergency fallback spillover file",
    )
    app_log_enabled: bool = Field(
        default=True,
        description="Enable rotating file logging for application logs",
    )
    app_log_path: str = Field(
        default="./logs/app.log",
        description="File path for rotating application log file",
    )
    app_log_max_bytes: int = Field(
        default=10 * 1024 * 1024,
        description="Maximum size in bytes before rotating application log file (default 10 MB)",
    )
    app_log_backup_count: int = Field(
        default=5,
        description="Number of rotated backup log files to retain",
    )

    # Static web application assets
    webapp_dist_dir: Path = Field(
        default=_REPO_ROOT / "webapp" / "dist",
        description=(
            "Directory of the built web application served at '/' (if it exists); "
            "relative paths resolve against the repository root"
        ),
    )

    # SMTP feedback delivery settings
    smtp_host: str = Field(
        default="",
        description="SMTP server hostname for feedback email delivery (empty disables SMTP)",
    )
    smtp_port: int = Field(default=587, description="SMTP server port")
    smtp_user: str = Field(default="", description="SMTP authentication username")
    smtp_password: str = Field(default="", description="SMTP authentication password")
    smtp_use_tls: bool = Field(default=True, description="Enable STARTTLS encryption for SMTP")
    feedback_recipient_email: str = Field(
        default="", description="Recipient email address for feedback notifications"
    )
    feedback_sender_email: str = Field(
        default="", description="Sender email address for feedback notifications"
    )

    @field_validator("webapp_dist_dir")
    @classmethod
    def _anchor_webapp_dist_dir(cls, value: Path) -> Path:
        return value if value.is_absolute() else _REPO_ROOT / value

    @property
    def parsed_cors_origins(self) -> List[str]:
        """Resolves cors_origins into a list of strings if provided as JSON or comma-separated string."""
        if isinstance(self.cors_origins, list):
            return self.cors_origins
        val = self.cors_origins.strip()
        if val.startswith("[") and val.endswith("]"):
            try:
                parsed = json.loads(val)
                if isinstance(parsed, list):
                    return [str(x) for x in parsed]
            except (json.JSONDecodeError, TypeError):
                pass
        return [item.strip() for item in val.split(",") if item.strip()]


def get_logging_config(settings: Settings | None = None) -> dict:
    """Generates standard dictConfig dictionary driven by application settings."""
    if settings is None:
        settings = get_settings()
    level = settings.log_level.upper()
    app_log_enabled = settings.app_log_enabled
    app_log_path = settings.app_log_path
    max_bytes = settings.app_log_max_bytes
    backup_count = settings.app_log_backup_count

    handlers_dict: dict = {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "standard",
            "stream": "ext://sys.stdout",
        },
    }

    active_handlers = ["console"]

    if app_log_enabled:
        log_file = Path(app_log_path)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers_dict["file"] = {
            "class": "logging.handlers.RotatingFileHandler",
            "formatter": "standard",
            "filename": str(log_file),
            "maxBytes": max_bytes,
            "backupCount": backup_count,
            "encoding": "utf-8",
        }
        active_handlers.append("file")

    return {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "standard": {
                "format": "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
        "handlers": handlers_dict,
        "loggers": {
            "api": {
                "handlers": list(active_handlers),
                "level": level,
                "propagate": False,
            },
            "pytransdatro": {
                "handlers": list(active_handlers),
                "level": level,
                "propagate": False,
            },
        },
        "root": {
            "handlers": list(active_handlers),
            "level": level,
        },
    }


def configure_logging(settings: Settings | None = None) -> None:
    """Initializes Python standard library logging via dictConfig."""
    import logging.config

    cfg = get_logging_config(settings)
    logging.config.dictConfig(cfg)


@lru_cache()
def get_settings() -> Settings:
    """Returns a singleton cached instance of the application Settings."""
    return Settings()
