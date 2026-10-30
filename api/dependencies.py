from pathlib import Path
from typing import Optional

from fastapi import Depends, Request

from api.client_sources import CLIENT_SOURCE_LOOKUP
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

    Decisional Structure:
    1. If 'X-Client-Source' header is missing or empty:
       Defaults to 1 (rest_api) to maintain seamless compatibility with standard
       HTTP clients, QGIS plugins, CAD scripts, and legacy consumers.
    2. If 'X-Client-Source' matches a registered origin in CLIENT_SOURCE_LOOKUP:
       Returns the configured integer source ID (e.g., 1 for rest_api, 2 for web_map,
       3 for desktop_cli, or 4+ for future partner sources).
    3. If 'X-Client-Source' is provided but unrecognized:
       Safely falls back to 0 (unknown).
    """
    raw_header = request.headers.get("X-Client-Source")
    if raw_header is None or not raw_header.strip():
        return 1

    normalized = raw_header.strip().lower()
    return CLIENT_SOURCE_LOOKUP.get(normalized, 0)


def get_db_path(settings: Settings = Depends(get_settings)) -> Path:
    """Resolves the SQLite telemetry database path from application settings."""
    return Path(settings.log_dir) / "usage.db"


def get_legacy_db_path(settings: Settings = Depends(get_settings)) -> Path:
    """Resolves the SQLite legacy TransdatOnline database path from application settings."""
    return settings.resolved_legacy_db_path
