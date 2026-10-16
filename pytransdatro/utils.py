"""Utility functions for angular, sexagesimal, and scale factor conversions.

Designed for internal geodetic processing within PyTransDatRO.
"""

import math

from pytransdatro import constants

# Re-export angular conversion factor for backward compatibility
DEG_TO_RAD_FACTOR: float = constants.DEG_TO_RAD_FACTOR
RAD_TO_DEG_FACTOR: float = constants.RAD_TO_DEG_FACTOR


def sexa_to_rad(s: str) -> float:
    """Returns the value in radians of a sexagesimal representation of an angle.

    The expected format of the representation is: ``[DD MM SS.SSS]`` or ``[-DD MM SS.SSS]``.

    :param s: Sexagesimal DMS angle string
    :type s: str
    :return: Angle value in radians
    :rtype: float

    .. warning::
        This function does not perform strict validation checks on input structure.
        It is designed for internal coordinate processing where valid formatting is guaranteed.
    """
    sgn = -1.0 if s[0] == "-" else 1.0
    values = s.split(" ")
    d = abs(float(values[0]))
    m = float(values[1])
    sec = float(values[2])
    return sgn * DEG_TO_RAD_FACTOR * (sec / 3600.0 + m / 60.0 + d)


def sexa_dms_chars_to_rad(s: str) -> float:
    """Returns the value in radians of a sexagesimal representation of an angle with DMS symbols.

    The format of the representation is: ``[DD°MM'SS.SSSSS"N]`` or with cardinal suffix characters.
    Standardizes character separators to spaces before conversion.

    :param s: DMS angle string with unit characters
    :type s: str
    :return: Angle value in radians
    :rtype: float
    """
    for c in "°'\"NESV":
        s = s.replace(c, " ", 1)
    return sexa_to_rad(s)


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
