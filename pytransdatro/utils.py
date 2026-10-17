"""Utility functions for angular, sexagesimal, and scale factor conversions.

Designed for internal geodetic processing within PyTransDatRO.
"""

import math
import re

from pytransdatro import constants, exceptions

# Re-export angular conversion factor for backward compatibility
DEG_TO_RAD_FACTOR: float = constants.DEG_TO_RAD_FACTOR
RAD_TO_DEG_FACTOR: float = constants.RAD_TO_DEG_FACTOR

_DMS_SPLIT_RE = re.compile(r"[°'\"`dms\s]+", re.IGNORECASE)


def sexa_to_rad(s: str) -> float:
    """Returns the value in radians of a space-separated sexagesimal representation of an angle.

    The expected format of the representation is: ``[DD [MM [SS.SSS]]]`` or ``[-DD [MM [SS.SSS]]]``.
    Enforces strict sexagesimal bounds (:math:`0 \\le m < 60`, :math:`0 \\le s < 60`) and finite values.

    :param s: Sexagesimal DMS angle string
    :type s: str
    :return: Angle value in radians
    :rtype: float
    :raises InvalidCoordinateErr: If the string is empty, non-numeric, non-finite, or minutes/seconds >= 60
    """
    return sexa_dms_chars_to_rad(s)


def sexa_dms_chars_to_rad(s: str) -> float:
    """Returns the value in radians of a sexagesimal representation of an angle with DMS symbols.

    Supports formats such as ``DD°MM'SS.SSSSS"N``, ``45d 59m 58.99s``, ``-45 30 00``,
    and cardinal indicators (``N``, ``E``, ``S``, ``W``, ``V``).
    Enforces strict sexagesimal bounds (:math:`0 \\le m < 60`, :math:`0 \\le s < 60`).

    :param s: DMS angle string with unit characters or whitespace
    :type s: str
    :return: Angle value in radians
    :rtype: float
    :raises InvalidCoordinateErr: If the string is empty, malformed, or minutes/seconds fall outside [0, 60)
    """
    if not isinstance(s, str):
        raise exceptions.InvalidCoordinateErr("DMS input must be a string.")
    raw = s.strip().replace("\ufffd", "°")
    if not raw:
        raise exceptions.InvalidCoordinateErr("Empty coordinate input.")

    sgn = 1.0
    clean = raw
    has_prefix_indicator = False
    has_suffix_indicator = False

    # Handle leading sign or leading cardinal indicator (N, E, S, W, V)
    if re.match(r"^-\s*", clean) or re.match(r"^[SWV]\s*\d", clean, re.IGNORECASE):
        sgn = -1.0
        has_prefix_indicator = True
        clean = re.sub(r"^-\s*|^[SWV]\s*", "", clean, count=1, flags=re.IGNORECASE)
    elif re.match(r"^\+\s*", clean) or re.match(r"^[NE]\s*\d", clean, re.IGNORECASE):
        has_prefix_indicator = True
        clean = re.sub(r"^\+\s*|^[NE]\s*", "", clean, count=1, flags=re.IGNORECASE)

    # Handle trailing cardinal indicator (S, W, V = negative; N, E = positive)
    if re.search(r"(?:\s+|[\"'`\d])[SWV]$", clean, re.IGNORECASE):
        has_suffix_indicator = True
        sgn = -1.0
        clean = re.sub(r"\s*[SWV]$", "", clean, count=1, flags=re.IGNORECASE)
    elif re.search(r"(?:\s+|[\"'`\d])[NE]$", clean, re.IGNORECASE):
        has_suffix_indicator = True
        clean = re.sub(r"\s*[NE]$", "", clean, count=1, flags=re.IGNORECASE)

    if has_prefix_indicator and has_suffix_indicator:
        raise exceptions.InvalidCoordinateErr(
            f"Conflicting prefix and suffix sign/cardinal indicators in '{s}'"
        )

    clean = clean.strip()
    if not clean or clean[0] in ("-", "+"):
        raise exceptions.InvalidCoordinateErr(f"Invalid coordinate input: '{s}'")

    tokens = [tok.strip() for tok in _DMS_SPLIT_RE.split(clean) if tok.strip()]
    if not tokens or len(tokens) > 3:
        raise exceptions.InvalidCoordinateErr(f"Invalid DMS token count in '{s}'")

    if len(tokens) > 1 and "." in tokens[0]:
        raise exceptions.InvalidCoordinateErr(
            f"Degrees must be integer when minutes/seconds are present in '{s}'"
        )
    if len(tokens) > 2 and "." in tokens[1]:
        raise exceptions.InvalidCoordinateErr(
            f"Minutes must be integer when seconds are present in '{s}'"
        )

    try:
        d_raw = float(tokens[0])
        m = float(tokens[1]) if len(tokens) > 1 else 0.0
        sec = float(tokens[2]) if len(tokens) > 2 else 0.0
    except ValueError as err:
        raise exceptions.InvalidCoordinateErr(f"Non-numeric DMS component in '{s}'") from err

    if not (math.isfinite(d_raw) and math.isfinite(m) and math.isfinite(sec)):
        raise exceptions.InvalidCoordinateErr(f"Non-finite DMS component in '{s}'")

    if len(tokens) > 1 and (tokens[1].startswith("-") or tokens[1].startswith("+")):
        raise exceptions.InvalidCoordinateErr(f"Minutes cannot be signed in '{s}'")
    if len(tokens) > 2 and (tokens[2].startswith("-") or tokens[2].startswith("+")):
        raise exceptions.InvalidCoordinateErr(f"Seconds cannot be signed in '{s}'")

    if m < 0.0 or m >= 60.0:
        raise exceptions.InvalidCoordinateErr(f"Minutes must be in [0, 60): {m}")
    if sec < 0.0 or sec >= 60.0:
        raise exceptions.InvalidCoordinateErr(f"Seconds must be in [0, 60): {sec}")

    d = abs(d_raw)
    return sgn * DEG_TO_RAD_FACTOR * (sec / 3600.0 + m / 60.0 + d)


def format_rad_to_dms(rad: float, precision: int = 5) -> str:
    """Formats an angle in radians into compacted sexagesimal DMS string ``DD°MM'SS.SSSSS"``.

    Includes automatic rollover when rounded seconds or minutes reach 60.

    :param rad: Angle in radians
    :type rad: float
    :param precision: Decimal places for seconds (default 5)
    :type precision: int
    :return: Formatted DMS string without spaces
    :rtype: str
    :raises InvalidCoordinateErr: If rad is NaN or infinite
    """
    if not math.isfinite(rad):
        raise exceptions.InvalidCoordinateErr(f"Cannot format non-finite radian value: {rad}")
    deg = math.degrees(rad)
    sign = "-" if deg < 0 else ""
    abs_deg = abs(deg)
    d = int(abs_deg)
    rem = (abs_deg - d) * 60.0
    m = int(rem)
    s = (rem - m) * 60.0
    factor = 10**precision
    s = round(s * factor) / factor
    if s >= 60.0:
        s = 0.0
        m += 1
        if m >= 60:
            m = 0
            d += 1
    s_str = f"{s:0{precision + 3}.{precision}f}" if precision > 0 else f"{round(s):02d}"
    return f"{sign}{d:d}°{m:02d}'{s_str}\""


def rad_to_sexa(a: float) -> str:
    """Converts an angle in radians to its sexagesimal space-separated string representation.

    Output format: ``[[-]DD MM SS.SSSSSS]``.

    :param a: Angle in radians
    :type a: float
    :return: Sexagesimal space-separated string
    :rtype: str
    """
    sgn = "-" if a < 0 else ""
    tmp_d = abs(a / DEG_TO_RAD_FACTOR)
    d = int(tmp_d)
    tmp_reminder = tmp_d - d
    m = int(tmp_reminder * 60)
    tmp_reminder = (tmp_reminder * 60) - m
    sec = tmp_reminder * 60
    return f"{sgn}{d} {m} {sec}"


def rad_to_deg(a: float) -> float:
    """Converts an angle in radians to decimal degrees.

    :param a: Angle in radians
    :type a: float
    :return: Angle in decimal degrees
    :rtype: float
    """
    return math.degrees(a)


def deg_to_rad(a: float) -> float:
    """Converts an angle in decimal degrees to radians.

    :param a: Angle in decimal degrees
    :type a: float
    :return: Angle in radians
    :rtype: float
    """
    return math.radians(a)


def ppm_to_unity(ppm: float) -> float:
    """Converts a parts-per-million (ppm) scale value to its dimensionless unity scale equivalent.

    :param ppm: Scale distortion in parts-per-million (ppm)
    :type ppm: float
    :return: Dimensionless scale factor ($1 + \text{ppm} \\cdot 10^{-6}$)
    :rtype: float
    """
    return 1.0 + ppm * 1e-6
