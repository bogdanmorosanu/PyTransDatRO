#!/usr/bin/env python3
"""Single-source-of-truth version management for PyTransDatRO.

Synchronizes the project version across:
  - pytransdatro/__init__.py (__version__)
  - pyproject.toml (project.version)
  - webapp/package.json (version)

Usage:
  python tools/bump_version.py --check
  python tools/bump_version.py --current
  python tools/bump_version.py 1.0.0-beta
  python tools/bump_version.py 1.0.0
"""

import json
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

SEMVER_REGEX = re.compile(
    r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)"
    r"(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    r"(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)

REPO_ROOT = Path(__file__).resolve().parent.parent
INIT_FILE = REPO_ROOT / "pytransdatro" / "__init__.py"
PYPROJECT_FILE = REPO_ROOT / "pyproject.toml"
PACKAGE_JSON_FILE = REPO_ROOT / "webapp" / "package.json"


def validate_semver(version: str) -> bool:
    """Validates that a version string strictly complies with SemVer 2.0.0."""
    return bool(SEMVER_REGEX.match(version.strip()))


def read_init_version() -> Optional[str]:
    """Reads __version__ from pytransdatro/__init__.py."""
    if not INIT_FILE.exists():
        return None
    content = INIT_FILE.read_text(encoding="utf-8")
    m = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', content, re.MULTILINE)
    return m.group(1).strip() if m else None


def read_pyproject_version() -> Optional[str]:
    """Reads version from pyproject.toml under [project]."""
    if not PYPROJECT_FILE.exists():
        return None
    content = PYPROJECT_FILE.read_text(encoding="utf-8")
    # Matches version = "..." under [project]
    m = re.search(r'^version\s*=\s*["\']([^"\']+)["\']', content, re.MULTILINE)
    return m.group(1).strip() if m else None


def read_package_json_version() -> Optional[str]:
    """Reads version from webapp/package.json."""
    if not PACKAGE_JSON_FILE.exists():
        return None
    try:
        data = json.loads(PACKAGE_JSON_FILE.read_text(encoding="utf-8"))
        return str(data.get("version", "")).strip() or None
    except Exception:
        return None


def get_all_versions() -> Dict[str, Optional[str]]:
    """Returns a dictionary of all recorded file versions."""
    return {
        "pytransdatro/__init__.py": read_init_version(),
        "pyproject.toml": read_pyproject_version(),
        "webapp/package.json": read_package_json_version(),
    }


def write_init_version(new_version: str) -> None:
    """Updates __version__ in pytransdatro/__init__.py."""
    content = INIT_FILE.read_text(encoding="utf-8")
    updated, count = re.subn(
        r'^(__version__\s*=\s*["\'])[^"\']+(["\'])',
        rf"\g<1>{new_version}\g<2>",
        content,
        count=1,
        flags=re.MULTILINE,
    )
    if count == 0:
        raise ValueError(f"Could not find __version__ assignment in {INIT_FILE}")
    INIT_FILE.write_text(updated, encoding="utf-8")


def write_pyproject_version(new_version: str) -> None:
    """Updates version in pyproject.toml."""
    content = PYPROJECT_FILE.read_text(encoding="utf-8")
    updated, count = re.subn(
        r'^(version\s*=\s*["\'])[^"\']+(["\'])',
        rf"\g<1>{new_version}\g<2>",
        content,
        count=1,
        flags=re.MULTILINE,
    )
    if count == 0:
        raise ValueError(f"Could not find version assignment under [project] in {PYPROJECT_FILE}")
    PYPROJECT_FILE.write_text(updated, encoding="utf-8")


def write_package_json_version(new_version: str) -> None:
    """Updates version in webapp/package.json while preserving indentation."""
    content = PACKAGE_JSON_FILE.read_text(encoding="utf-8")
    data = json.loads(content)
    data["version"] = new_version
    PACKAGE_JSON_FILE.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def check_versions() -> Tuple[bool, str]:
    """Validates that all versions exist, match, and adhere to SemVer."""
    versions = get_all_versions()
    missing = [f for f, v in versions.items() if not v]
    if missing:
        return False, f"Missing version in file(s): {', '.join(missing)}"

    distinct = set(versions.values())
    if len(distinct) > 1:
        details = "\n".join(f"  - {f}: {v}" for f, v in versions.items())
        return False, f"Version mismatch detected across files:\n{details}"

    current = next(iter(distinct))
    if not validate_semver(current):
        return (
            False,
            f"Version '{current}' does not comply with SemVer 2.0.0 (MAJOR.MINOR.PATCH[-PRERELEASE])",
        )

    return True, f"All versions synchronized and SemVer compliant: {current}"


def set_version(new_version: str) -> None:
    """Updates all 3 version files atomically to new_version."""
    clean_version = new_version.strip().lstrip("v")
    if not validate_semver(clean_version):
        sys.exit(
            f"Error: '{clean_version}' is not a valid SemVer 2.0.0 string (e.g. 1.0.0, 1.0.0-beta)."
        )

    print(f"Setting version to: {clean_version}")
    write_init_version(clean_version)
    print("  [OK] Updated pytransdatro/__init__.py")
    write_pyproject_version(clean_version)
    print("  [OK] Updated pyproject.toml")
    write_package_json_version(clean_version)
    print("  [OK] Updated webapp/package.json")

    ok, msg = check_versions()
    if not ok:
        sys.exit(f"Post-update verification failed: {msg}")
    print(f"Success: {msg}")


def main() -> None:
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0)

    if args[0] == "--check":
        ok, msg = check_versions()
        if ok:
            print(f"OK: {msg}")
            sys.exit(0)
        else:
            print(f"FAIL: {msg}", file=sys.stderr)
            sys.exit(1)

    if args[0] == "--current":
        versions = get_all_versions()
        print(versions.get("pytransdatro/__init__.py") or "unknown")
        sys.exit(0)

    target_version = args[0]
    set_version(target_version)


if __name__ == "__main__":
    main()
