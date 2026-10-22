"""Unit tests for centralized configuration and core constants."""

import math

from api.config import Settings
from pytransdatro import constants, utils
from pytransdatro.proj_stereo import StereoProj


def test_default_settings():
    """Verifies that Settings instantiate with expected default parameters."""
    settings = Settings()
    assert settings.host == "0.0.0.0"
    assert settings.port == 8000
    assert settings.log_level == "INFO"
    assert settings.cors_origins == ["*"]
    assert settings.parsed_cors_origins == ["*"]
    assert settings.rate_limit_window_seconds == 60.0
    assert settings.rate_limit_point == 60
    assert settings.rate_limit_batch == 10
    assert settings.max_batch_points == 20000
    assert settings.max_file_points == 100000
    assert settings.max_file_upload_bytes == 50 * 1024 * 1024
    assert settings.log_dir == "./logs"
    assert settings.feedback_log_path == "./logs/feedback.jsonl"
    assert settings.telemetry_fallback_path == "./logs/telemetry_fallback.jsonl"


def test_settings_cors_parsing():
    """Verifies that cors_origins parses both JSON arrays and comma-separated strings."""
    # JSON list string
    s_json = Settings(cors_origins='["https://example.com", "http://localhost:3000"]')
    assert s_json.parsed_cors_origins == ["https://example.com", "http://localhost:3000"]

    # Comma-separated string
    s_csv = Settings(cors_origins="https://example.com, http://localhost:5173")
    assert s_csv.parsed_cors_origins == ["https://example.com", "http://localhost:5173"]

    # Python list
    s_list = Settings(cors_origins=["https://custom.app"])
    assert s_list.parsed_cors_origins == ["https://custom.app"]


def test_settings_env_override(monkeypatch):
    """Verifies that PYTRANSDAT_ environment variables override defaults."""
    monkeypatch.setenv("PYTRANSDAT_PORT", "9050")
    monkeypatch.setenv("PYTRANSDAT_LOG_LEVEL", "DEBUG")
    monkeypatch.setenv("PYTRANSDAT_MAX_BATCH_POINTS", "5000")

    settings = Settings()
    assert settings.port == 9050
    assert settings.log_level == "DEBUG"
    assert settings.max_batch_points == 5000


def test_geodetic_constants_values():
    """Verifies exact physical, ellipsoid, and projection constants."""
    assert math.isclose(constants.DEG_TO_RAD_FACTOR, math.pi / 180.0, rel_tol=1e-15)
    assert math.isclose(constants.RAD_TO_DEG_FACTOR, 180.0 / math.pi, rel_tol=1e-15)

    # WGS84 Ellipsoid
    assert constants.WGS84_ELL_SMAJ_AXIS == 6378137.0
    assert constants.WGS84_ELL_INV_FLAT == 298.257223563

    # Stereo 70 (origin in radians, matching 46° 00' 00.0" N and 25° 00' 00.0" E)
    assert constants.STEREO70_ORIG_LAT == 46.0 * constants.DEG_TO_RAD_FACTOR
    assert constants.STEREO70_ORIG_LON == 25.0 * constants.DEG_TO_RAD_FACTOR
    assert constants.STEREO70_ORIG_LAT == utils.sexa_to_rad("46 0 0.0")
    assert constants.STEREO70_ORIG_LON == utils.sexa_to_rad("25 0 0.0")
    assert constants.STEREO70_FALSE_N == 500000.0
    assert constants.STEREO70_FALSE_E == 500000.0
    assert constants.STEREO70_SCALE_FACTOR == 0.99975
    assert constants.STEREO70_CONV_TOLERANCE == 0.0000000000484814
    assert constants.STEREO70_MAX_ITERATIONS == 50

    # Romania Envelope Bounding Box
    assert constants.ROMANIA_BBOX_LAT_MIN == 43.5
    assert constants.ROMANIA_BBOX_LAT_MAX == 48.3
    assert constants.ROMANIA_BBOX_LON_MIN == 20.2
    assert constants.ROMANIA_BBOX_LON_MAX == 29.8


def test_constants_backward_compatibility():
    """Verifies backward-compatible re-exports on StereoProj and utils."""
    assert StereoProj.WGS84_ELL_SMAJ_AXIS == constants.WGS84_ELL_SMAJ_AXIS
    assert StereoProj.WGS84_ELL_INV_FLAT == constants.WGS84_ELL_INV_FLAT
    assert StereoProj.ORIG_LAT == constants.STEREO70_ORIG_LAT
    assert StereoProj.ORIG_LON == constants.STEREO70_ORIG_LON
    assert StereoProj.FALSE_N == constants.STEREO70_FALSE_N
    assert StereoProj.FALSE_E == constants.STEREO70_FALSE_E
    assert StereoProj.PROJ_SCALE == constants.STEREO70_SCALE_FACTOR

    assert utils.DEG_TO_RAD_FACTOR == constants.DEG_TO_RAD_FACTOR


def test_smtp_settings_defaults():
    """Verifies SMTP settings defaults in Settings."""
    settings = Settings()
    assert settings.smtp_host == ""
    assert settings.smtp_port == 587
    assert settings.smtp_user == ""
    assert settings.smtp_password == ""
    assert settings.smtp_use_tls is True
    assert settings.feedback_recipient_email == ""
    assert settings.feedback_sender_email == ""


def test_logging_configuration():
    """Verifies that get_logging_config and configure_logging initialize correctly."""
    from api.config import configure_logging, get_logging_config

    settings = Settings(log_level="DEBUG")
    cfg = get_logging_config(settings)

    assert cfg["version"] == 1
    assert "formatters" in cfg
    assert "handlers" in cfg
    assert cfg["loggers"]["api"]["level"] == "DEBUG"
    assert cfg["loggers"]["pytransdatro"]["level"] == "DEBUG"

    # Should configure without errors
    configure_logging(settings)


def test_webapp_dist_dir_default_and_override(monkeypatch, tmp_path):
    """Verifies the static web assets directory defaults to webapp/dist and is overridable."""
    from pathlib import Path

    repo_root = Path(__file__).resolve().parent.parent
    monkeypatch.delenv("PYTRANSDAT_WEBAPP_DIST_DIR", raising=False)
    # _env_file=None: a developer's local .env must not leak into the default check.
    assert Settings(_env_file=None).webapp_dist_dir == repo_root / "webapp" / "dist"

    monkeypatch.setenv("PYTRANSDAT_WEBAPP_DIST_DIR", str(tmp_path))
    assert Settings().webapp_dist_dir == tmp_path

    # Relative overrides anchor to the repo root, independent of the working directory.
    monkeypatch.setenv("PYTRANSDAT_WEBAPP_DIST_DIR", "webapp/dist-e2e")
    monkeypatch.chdir(tmp_path)
    assert Settings().webapp_dist_dir == repo_root / "webapp" / "dist-e2e"
