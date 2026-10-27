"""Client origin sources registry and lookup configuration.

Defines the mapping between incoming HTTP request identifiers (via the
'X-Client-Source' header) and the integer source IDs persisted in
the telemetry SQLite database (usage.db).

Baseline Source IDs (0 to 3) are system constants initialized in usage.db:
    0: unknown     - Unspecified origin or unrecognized client header
    1: rest_api    - REST API Service (default for callers without header)
    2: web_map     - Web Map Application (the official web frontend)
    3: desktop_cli - Desktop applications or CLI tools

To register a new community / partner source in the future:
    1. Add a new tuple (ID, name, description) to REGISTERED_SOURCES (e.g. ID 4).
    2. Add the expected lowercase header alias(es) to CLIENT_SOURCE_LOOKUP.
"""

from typing import Dict, List, Tuple

# Baseline dimension entries registered in the SQLite 'sources' table:
# Format: (id, name, description)
REGISTERED_SOURCES: List[Tuple[int, str, str]] = [
    # --- System Baselines (Do NOT modify these IDs) ---
    (0, "unknown", "Default / Unspecified origin"),
    (1, "rest_api", "REST API Service"),
    (2, "web_map", "Web Map Application"),
    (3, "desktop_cli", "Desktop / CLI App"),
    # --- Future community / partner sources can be appended here ---
    # (4, "partner_portal", "Geo-Spatial Partner Portal"),
]

# Canonical lookup mapping from incoming header strings to numerical source IDs:
CLIENT_SOURCE_LOOKUP: Dict[str, int] = {
    # 0: Explicit unknown
    "unknown": 0,
    # 1: REST API clients
    "rest_api": 1,
    "api": 1,
    "rest": 1,
    # 2: Official Web Map application
    "web_map": 2,
    "web_app": 2,
    "web": 2,
    # 3: Desktop and CLI tools
    "desktop_cli": 3,
    "cli": 3,
    "desktop": 3,
    # --- Future community / partner source aliases ---
    # "partner_portal": 4,
    # "partner": 4,
}
