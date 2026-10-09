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
