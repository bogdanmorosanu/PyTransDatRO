"""FastAPI dependency providers for accessing application state.
"""
from fastapi import Request
from pytransdatro.trans_ro import TransRO
from pytransdatro.telemetry import SqliteUsageLogger
from typing import Optional


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
