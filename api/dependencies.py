from pathlib import Path
from typing import Optional

from fastapi import Depends, Request

from api.config import Settings, get_settings
from pytransdatro.telemetry import SqliteUsageLogger
from pytransdatro.trans_ro import TransRO

__all__ = [
    "Settings",
    "get_client_source",
    "get_db_path",
    "get_legacy_db_path",
    "get_settings",
    "get_telemetry_logger",
    "get_trans_engine",
]


def get_trans_engine(request: Request) -> TransRO:
    """Retrieve the singleton TransRO coordinate transformation engine from app state."""
    return request.app.state.trans


def get_telemetry_logger(request: Request) -> Optional[SqliteUsageLogger]:
    """Retrieve the optional SqliteUsageLogger from app state."""
    return getattr(request.app.state, "telemetry", None)


def get_client_source(request: Request) -> int:
    """Resolves client origin identifier from request headers.

    Returns:
        1: rest_api (default for API endpoints)
        2: web_map (when X-Client-Source indicates web_map or web_app)
        3: desktop_cli (when X-Client-Source indicates desktop_cli)
    """
    client_src = request.headers.get("X-Client-Source", "").strip().lower()
    if client_src in ("web_map", "web_app", "web"):
        return 2
    if client_src in ("desktop_cli", "cli", "desktop"):
        return 3
    return 1


def get_db_path(settings: Settings = Depends(get_settings)) -> Path:
    """Resolves the SQLite telemetry database path from application settings."""
    return Path(settings.log_dir) / "usage.db"


def get_legacy_db_path(settings: Settings = Depends(get_settings)) -> Path:
    """Resolves the SQLite legacy TransdatOnline database path from application settings."""
    return Path(settings.legacy_db_path)
