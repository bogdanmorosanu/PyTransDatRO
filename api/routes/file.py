"""Streaming delimited file transformation endpoint.
Supports up to 100,000 coordinates in CSV/TXT formats with robust Point ID and DMS handling.
"""

import math
from typing import AsyncGenerator, List, Optional

from fastapi import APIRouter, Depends, File, Form, UploadFile
from fastapi.responses import StreamingResponse

import pytransdatro.exceptions as exc
import pytransdatro.utils as geoutils
from api.config import Settings, get_settings
from api.dependencies import get_client_source, get_trans_engine
from api.schemas.point import AngleUnitEnum, OperationEnum
from pytransdatro.trans_ro import TransRO

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])

DELIMITER_MAP = {
    "comma": ",",
    "semicolon": ";",
    "tab": "\t",
    "space": " ",
}


def _detect_delimiter(line: str) -> str:
    """Auto-detects the delimiter in a line from common candidates."""
    for delim in [";", "\t", ","]:
        if delim in line:
            return delim
    return " "


def rad_to_dms_str(rad: float, precision: int = 5) -> str:
    """Format radians into DD°MM'SS.SSSSS" with no spaces."""
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
    s_str = f"{s:0{precision + 3}.{precision}f}"
    return f"{sign}{d:d}°{m:02d}'{s_str}\""


def _parse_coordinate_token(token: str, is_geo: bool, is_degrees: bool) -> float:
    t = token.strip().strip('"').strip("'")
    if is_geo:
        if any(c in t for c in "°'\"\ufffd") or len(t.split()) == 3:
            clean_t = t.replace("\ufffd", "°")
            return geoutils.sexa_dms_chars_to_rad(clean_t)
        val = float(t)
        return math.radians(val) if is_degrees else val
    return float(t)


def _transform_line(
    line: str,
    delim_char: str,
    op: OperationEnum,
    unit: AngleUnitEnum,
    trans: TransRO,
    source: int = 1,
) -> str:
    # Split tokens
    tokens = [t.strip() for t in (line.split() if delim_char == " " else line.split(delim_char))]
    tokens = [t for t in tokens if t]

    if len(tokens) < 2:
        return line + " # Wrong dimension\n"

    is_source_geo = op == OperationEnum.ETRS89_TO_STEREO70
    is_degrees = unit == AngleUnitEnum.DEGREES

    point_id: Optional[str] = None
    coord_tokens: List[str] = tokens

    if len(tokens) >= 4:
        point_id = tokens[0]
        coord_tokens = tokens[1:4]
    elif len(tokens) == 3:
        t0, t1, t2 = tokens
        if is_source_geo:
            has_angle0 = any(c in t0 for c in "°'\"")
            has_angle1 = any(c in t1 for c in "°'\"")
            try:
                num0 = float(t0.strip("\"'"))
            except ValueError:
                num0 = None
            if not has_angle0 and (
                has_angle1 or (num0 is not None and abs(num0) > 90) or num0 is None
            ):
                point_id = t0
                coord_tokens = [t1, t2]
            else:
                coord_tokens = tokens
        else:
            try:
                n0 = float(t0.strip("\"'"))
                n1 = float(t1.strip("\"'"))
                n2 = float(t2.strip("\"'"))
                if n1 > 50000 and n2 > 50000 and n0 < 50000:
                    point_id = t0
                    coord_tokens = [t1, t2]
                else:
                    coord_tokens = tokens
            except ValueError:
                point_id = t0
                coord_tokens = [t1, t2]

    try:
        c1 = _parse_coordinate_token(coord_tokens[0], is_source_geo, is_degrees)
        c2 = _parse_coordinate_token(coord_tokens[1], is_source_geo, is_degrees)
        c3 = float(coord_tokens[2].strip("\"'")) if len(coord_tokens) >= 3 else None
    except (ValueError, TypeError, IndexError):
        return line + " # Invalid coordinate data\n"

    dim = len(coord_tokens)

    try:
        if op == OperationEnum.STEREO70_TO_ETRS89:
            res = trans.st70_to_etrs89(c1, c2, z=c3, source=source)
            if unit == AngleUnitEnum.DMS:
                lat_str = rad_to_dms_str(res[0])
                lon_str = rad_to_dms_str(res[1])
                out_coords = f"{lat_str}{delim_char}{lon_str}"
            elif unit == AngleUnitEnum.DEGREES:
                out_coords = f"{math.degrees(res[0]):.10f}{delim_char}{math.degrees(res[1]):.10f}"
            else:
                out_coords = f"{res[0]:.16f}{delim_char}{res[1]:.16f}"

            if dim >= 3:
                out_coords += f"{delim_char}{res[2]:.4f}"

        else:
            res = trans.etrs89_to_st70(c1, c2, h=c3, source=source)
            out_coords = f"{res[0]:.4f}{delim_char}{res[1]:.4f}"
            if dim >= 3:
                out_coords += f"{delim_char}{res[2]:.4f}"

        if point_id is not None:
            return f"{point_id}{delim_char}{out_coords}\n"
        return f"{out_coords}\n"

    except exc.OutOfGridErr:
        return line + " # Out of grid\n"
    except exc.NoDataGridErr:
        return line + " # No data on grid\n"


@router.post("/file")
async def transform_file(
    file: UploadFile = File(..., description="CSV or TXT file containing coordinates"),
    op: OperationEnum = Form(..., description="Transformation direction"),
    unit: AngleUnitEnum = Form(
        default=AngleUnitEnum.RADIANS, description="Angular unit for ETRS89 coordinates"
    ),
    delimiter: str = Form(
        default="auto", description="Delimiter: auto, comma, semicolon, tab, or space"
    ),
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
    settings: Settings = Depends(get_settings),
):
    """Transforms a delimited file (CSV, TXT) streaming results back using chunked HTTP transfer.

    Processes points with negligible memory overhead up to the configured limit.
    """
    chosen_delim = DELIMITER_MAP.get(delimiter.lower())

    async def line_generator() -> AsyncGenerator[bytes, None]:
        nonlocal chosen_delim
        row_count = 0
        max_rows = settings.max_file_points

        # Read line by line from uploaded stream
        buffer = ""
        while True:
            chunk = await file.read(64 * 1024)  # 64 KB chunks
            if not chunk:
                break
            text = chunk.decode("utf-8", errors="replace")
            buffer += text
            lines = buffer.splitlines(keepends=True)
            if lines:
                buffer = ""
                # If last line doesn't end with a newline, buffer it
                if not lines[-1].endswith(("\n", "\r")):
                    buffer = lines.pop()

                for raw_line in lines:
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith("#"):
                        yield (line + "\n").encode("utf-8")
                        continue

                    # Auto-detect delimiter on first valid data line if needed
                    if chosen_delim is None:
                        chosen_delim = _detect_delimiter(line)

                    delim_char = chosen_delim if chosen_delim != " " else " "
                    row_count += 1
                    if row_count > max_rows:
                        yield f"# Limit exceeded: maximum {max_rows:,} points allowed\n".encode(
                            "utf-8"
                        )
                        return

                    out = _transform_line(line, delim_char, op, unit, trans, source=source)
                    yield out.encode("utf-8")

        # Process any remaining line in buffer
        if buffer.strip():
            line = buffer.strip()
            if line.startswith("#"):
                yield (line + "\n").encode("utf-8")
            else:
                delim_char = chosen_delim or ","
                out = _transform_line(line, delim_char, op, unit, trans, source=source)
                yield out.encode("utf-8")

    return StreamingResponse(
        line_generator(),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="transformed_coordinates.csv"'},
    )
