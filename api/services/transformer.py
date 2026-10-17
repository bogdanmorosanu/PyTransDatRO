"""Unified coordinate validation, unit conversion, and exception-to-warning adapter.

Centralizes API-level validation and transformation dispatch across Legacy (/cooOpService)
and Modern (/api/v1/transform/*) endpoints while delegating geodetic math and domain
checks to pytransdatro.TransRO.
"""

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pytransdatro.exceptions as exc
from api.schemas.point import AngleUnitEnum, OperationEnum, PointResponse
from pytransdatro.trans_ro import TransRO

LEGACY_SUPPORTED_OPS = {
    "Stereo70ToETRS89",
    "ETRS89ToStereo70",
    "Stereo30ToETRS89",
    "ETRS89ToStereo30",
}

STEREO30_WARNING = "Stereo30 transformation is obsolete and unsupported"
WARNING_WRONG_DIMENSION = "Wrong dimension"
WARNING_INVALID_COORDINATE = "Invalid coordinate data"
WARNING_OUT_OF_GRID = "Out of grid"
WARNING_NO_DATA = "No data on grid"


def validate_legacy_operation(coo_op: str) -> None:
    """Validates that the legacy operation name is recognized.

    :param coo_op: Operation identifier string
    :type coo_op: str
    :raises ValueError: If coo_op is not a defined legacy transformation
    """
    if coo_op not in LEGACY_SUPPORTED_OPS:
        raise ValueError(f"There is no defined transformation with name {coo_op}")


def parse_legacy_semicolon_coos(coos_str: str) -> Tuple[List[float], Optional[str]]:
    """Parses a semicolon-delimited coordinate string replicating Java String.split(';') semantics.

    Trailing empty segments (e.g. ``500000;500000;``) are stripped, whereas leading or
    intermediate empty segments (e.g. ``;500000`` or ``500000;;500000``) are rejected as
    invalid coordinate data.

    :param coos_str: Raw semicolon-separated coordinate query string
    :type coos_str: str
    :return: Tuple of ``(parsed_floats, warning_or_none)``
    :rtype: Tuple[List[float], Optional[str]]
    """
    raw_tokens = coos_str.split(";")
    while raw_tokens and raw_tokens[-1].strip() == "":
        raw_tokens.pop()

    if not raw_tokens:
        return [], WARNING_INVALID_COORDINATE

    parsed_coords: List[float] = []
    for token in raw_tokens:
        stripped = token.strip()
        if not stripped:
            return [], WARNING_INVALID_COORDINATE
        try:
            val = float(stripped)
        except ValueError:
            return [], WARNING_INVALID_COORDINATE
        if not math.isfinite(val):
            return [], WARNING_INVALID_COORDINATE
        parsed_coords.append(val)

    return parsed_coords, None


def parse_numeric_coord_list(raw_coords: Any) -> Tuple[List[float], Optional[str]]:
    """Validates and converts a sequence of raw coordinate values into finite floats.

    :param raw_coords: Input sequence of coordinate values
    :type raw_coords: Any
    :return: Tuple of ``(numeric_coords, warning_or_none)``
    :rtype: Tuple[List[float], Optional[str]]
    """
    if not isinstance(raw_coords, (list, tuple)):
        return [], WARNING_INVALID_COORDINATE

    numeric_coords: List[float] = []
    for c in raw_coords:
        if isinstance(c, bool) or c is None:
            return [], WARNING_INVALID_COORDINATE
        try:
            val = float(c)
        except (ValueError, TypeError):
            return [], WARNING_INVALID_COORDINATE
        if not math.isfinite(val):
            return [], WARNING_INVALID_COORDINATE
        numeric_coords.append(val)

    return numeric_coords, None


def execute_point_transformation(
    trans: TransRO,
    op: str,
    coords: Sequence[float],
    unit: AngleUnitEnum = AngleUnitEnum.RADIANS,
    source: int = 1,
) -> Tuple[List[float], Optional[str]]:
    """Validates coordinate dimensionality and executes single-point transformation.

    Enforces:
        -Finite numeric coordinate values (rejecting NaN / Infinity with empty coos).
        -Strict 2D or 3D dimensionality (:math:`2 \\le \\text{dim} \\le 3`), echoing input
         coordinates with ``"Wrong dimension"`` when :math:`\\text{dim} < 2` or :math:`\\text{dim} > 3`.
        -Unit conversion between radians and decimal degrees.
        -Domain exception mapping (``OutOfGridErr``, ``NoDataGridErr``, ``InvalidCoordinateErr``,
         ``WrongDimensionErr``) to standard warning strings.

    :param trans: Pre-warmed TransRO instance
    :type trans: TransRO
    :param op: Operation name (e.g. Stereo70ToETRS89, ETRS89ToStereo70, Stereo30ToETRS89)
    :type op: str
    :param coords: Sequence of numeric coordinate values
    :type coords: Sequence[float]
    :param unit: Angular unit for ETRS89 coordinates (radians or degrees)
    :type unit: AngleUnitEnum
    :param source: Telemetry client source ID
    :type source: int
    :return: Tuple of ``(output_coords, warning_or_none)``
    :rtype: Tuple[List[float], Optional[str]]
    """
    coord_list, parse_warn = parse_numeric_coord_list(coords)
    if parse_warn is not None:
        return [], parse_warn

    dim = len(coord_list)
    if dim < 2 or dim > 3:
        return coord_list, WARNING_WRONG_DIMENSION

    if op in ("Stereo30ToETRS89", "ETRS89ToStereo30"):
        return coord_list, STEREO30_WARNING

    is_degrees = unit == AngleUnitEnum.DEGREES

    try:
        if op == OperationEnum.STEREO70_TO_ETRS89.value:
            n, e = coord_list[0], coord_list[1]
            z = coord_list[2] if dim == 3 else None
            res = trans.st70_to_etrs89(n, e, z=z, source=source)

            if is_degrees:
                out_coords = [math.degrees(res[0]), math.degrees(res[1])]
                if dim == 3:
                    out_coords.append(res[2])  # type: ignore[index]
                return out_coords, None

            return list(res), None

        elif op == OperationEnum.ETRS89_TO_STEREO70.value:
            lat, lon = coord_list[0], coord_list[1]
            if is_degrees:
                lat = math.radians(lat)
                lon = math.radians(lon)

            h = coord_list[2] if dim == 3 else None
            res = trans.etrs89_to_st70(lat, lon, h=h, source=source)
            return list(res), None

        else:
            raise ValueError(f"There is no defined transformation with name {op}")

    except exc.OutOfGridErr:
        return coord_list, WARNING_OUT_OF_GRID
    except exc.NoDataGridErr:
        return coord_list, WARNING_NO_DATA
    except exc.WrongDimensionErr:
        return coord_list, WARNING_WRONG_DIMENSION
    except exc.InvalidCoordinateErr:
        return [], WARNING_INVALID_COORDINATE


def _map_domain_exception_to_warning(domain_err: exc.PyTransDatError) -> str:
    """Maps a pytransdatro domain exception to its canonical API warning string."""
    if isinstance(domain_err, exc.OutOfGridErr):
        return WARNING_OUT_OF_GRID
    if isinstance(domain_err, exc.NoDataGridErr):
        return WARNING_NO_DATA
    if isinstance(domain_err, exc.WrongDimensionErr):
        return WARNING_WRONG_DIMENSION
    return WARNING_INVALID_COORDINATE


def execute_batch_transformation(
    trans: TransRO,
    op: str,
    points: Sequence[Any],
    unit: AngleUnitEnum = AngleUnitEnum.RADIANS,
    source: int = 1,
) -> List[Tuple[List[float], Optional[str]]]:
    """Executes fault-isolated batch coordinate transformation via TransRO.transform_batch.

    :param trans: Pre-warmed TransRO instance
    :type trans: TransRO
    :param op: Operation name (e.g. Stereo70ToETRS89, ETRS89ToStereo70, Stereo30ToETRS89)
    :type op: str
    :param points: Sequence of coordinate rows
    :type points: Sequence[Any]
    :param unit: Angular unit for ETRS89 coordinates (radians or degrees)
    :type unit: AngleUnitEnum
    :param source: Telemetry client source ID
    :type source: int
    :return: Ordered list of ``(output_coords, warning_or_none)`` tuples
    :rtype: List[Tuple[List[float], Optional[str]]]
    """
    if op in ("Stereo30ToETRS89", "ETRS89ToStereo30"):
        st30_results: List[Tuple[List[float], Optional[str]]] = []
        for raw_pt in points:
            coord_list, parse_warn = parse_numeric_coord_list(raw_pt)
            if parse_warn is not None:
                st30_results.append(([], parse_warn))
            elif len(coord_list) < 2 or len(coord_list) > 3:
                st30_results.append((coord_list, WARNING_WRONG_DIMENSION))
            else:
                st30_results.append((coord_list, STEREO30_WARNING))
        return st30_results

    if op not in (
        OperationEnum.STEREO70_TO_ETRS89.value,
        OperationEnum.ETRS89_TO_STEREO70.value,
    ):
        raise ValueError(f"There is no defined transformation with name {op}")

    from_st70 = op == OperationEnum.STEREO70_TO_ETRS89.value
    is_degrees = unit == AngleUnitEnum.DEGREES

    raw_results = trans.transform_batch(
        points, from_st70=from_st70, is_degrees=is_degrees, source=source
    )

    out: List[Tuple[List[float], Optional[str]]] = []
    for coords_out, domain_err in raw_results:
        if domain_err is None:
            out.append((coords_out, None))
        else:
            out.append((coords_out, _map_domain_exception_to_warning(domain_err)))
    return out


def format_legacy_dict(coos: List[float], warning: Optional[str]) -> Dict[str, Any]:
    """Formats a transformation result into the legacy TransDatOnline dictionary contract.

    Omits the ``warning`` key when ``warning`` is None.
    """
    if warning is None:
        return {"coos": coos}
    return {"coos": coos, "warning": warning}


def format_modern_response(coos: List[float], warning: Optional[str]) -> PointResponse:
    """Formats a transformation result into a modern PointResponse schema instance."""
    return PointResponse(coos=coos, warning=warning)
