"""Streaming delimited file transformation endpoint.
Supports up to 100,000 coordinates in CSV/TXT formats with robust Point ID and DMS handling.
"""

import logging
import math
from typing import AsyncGenerator, List, Optional

from fastapi import APIRouter, Depends, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import StreamingResponse

import pytransdatro.exceptions as exc
import pytransdatro.utils as geoutils
from api.config import Settings, get_settings
from api.dependencies import get_client_source, get_trans_engine
from api.schemas.point import FileAngleUnitEnum, OperationEnum
from pytransdatro.telemetry import BufferedTelemetryCollector
from pytransdatro.trans_ro import TransRO

logger = logging.getLogger("api.routes.file")

router = APIRouter(prefix="/api/v1/transform", tags=["Coordinate Transformation"])

DELIMITER_MAP = {
    "comma": ",",
    "semicolon": ";",
    "tab": "\t",
    "space": " ",
}

ALLOWED_DELIMITERS = {"auto", *DELIMITER_MAP.keys()}


def _detect_delimiter(line: str) -> str:
    """Auto-detects the delimiter in a line from common candidates."""
    for delim in [";", "\t", ","]:
        if delim in line:
            return delim
    return " "


def rad_to_dms_str(rad: float, precision: int = 5) -> str:
    """Format radians into DD°MM'SS.SSSSS" with no spaces."""
    return geoutils.format_rad_to_dms(rad, precision=precision)


def _strip_balanced_quotes(token: str) -> str:
    """Strips outer CSV quotes only when symmetrically wrapped, preserving trailing DMS seconds (\")."""
    t = token.strip()
    if len(t) >= 2 and ((t[0] == '"' and t[-1] == '"') or (t[0] == "'" and t[-1] == "'")):
        return t[1:-1].strip()
    return t


def _parse_coordinate_token(token: str, is_geo: bool, is_degrees: bool) -> float:
    t = _strip_balanced_quotes(token)
    if is_geo:
        if any(c in t for c in "°'\"\ufffd") or len(t.split()) == 3:
            clean_t = t.replace("\ufffd", "°")
            return geoutils.sexa_dms_chars_to_rad(clean_t)
        val = float(t)
        if not math.isfinite(val):
            raise exc.InvalidCoordinateErr(f"Non-finite coordinate: {t}")
        return math.radians(val) if is_degrees else val
    val = float(t)
    if not math.isfinite(val):
        raise exc.InvalidCoordinateErr(f"Non-finite coordinate: {t}")
    return val


def _transform_line(
    line: str,
    delim_char: str,
    op: OperationEnum,
    unit: FileAngleUnitEnum,
    trans: TransRO,
    source: int = 1,
) -> str:
    # Split tokens
    if delim_char == " ":
        tokens = [t.strip() for t in line.split() if t.strip()]
    else:
        tokens = [t.strip() for t in line.split(delim_char)]
        while tokens and tokens[-1] == "":
            tokens.pop()
        if any(t == "" for t in tokens):
            return line + " # Invalid coordinate data\n"

    if len(tokens) < 2 or len(tokens) > 4:
        return line + " # Wrong dimension\n"

    is_source_geo = op == OperationEnum.ETRS89_TO_STEREO70
    is_degrees = unit in (FileAngleUnitEnum.DEGREES, FileAngleUnitEnum.DMS)

    point_id: Optional[str] = None
    coord_tokens: List[str] = tokens

    if len(tokens) == 4:
        point_id = tokens[0]
        coord_tokens = tokens[1:4]
    elif len(tokens) == 3:
        t0, t1, t2 = tokens
        if is_source_geo:
            has_angle0 = any(c in t0 for c in "°'\"")
            has_angle1 = any(c in t1 for c in "°'\"")
            try:
                num0 = float(t0.strip("\"'"))
                if not math.isfinite(num0):
                    return line + " # Invalid coordinate data\n"
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
                if not (math.isfinite(n0) and math.isfinite(n1) and math.isfinite(n2)):
                    return line + " # Invalid coordinate data\n"
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
        c3: Optional[float] = None
        if len(coord_tokens) == 3:
            c3 = float(coord_tokens[2].strip("\"'"))
            if not math.isfinite(c3):
                raise exc.InvalidCoordinateErr("Non-finite elevation coordinate")
    except (ValueError, TypeError, IndexError, exc.InvalidCoordinateErr):
        return line + " # Invalid coordinate data\n"

    dim = len(coord_tokens)

    try:
        if op == OperationEnum.STEREO70_TO_ETRS89:
            res = trans.st70_to_etrs89(c1, c2, z=c3, source=source)
            if unit == FileAngleUnitEnum.DMS:
                lat_str = rad_to_dms_str(res[0])
                lon_str = rad_to_dms_str(res[1])
                out_coords = f"{lat_str}{delim_char}{lon_str}"
            elif unit == FileAngleUnitEnum.DEGREES:
                out_coords = f"{math.degrees(res[0]):.10f}{delim_char}{math.degrees(res[1]):.10f}"
            else:
                out_coords = f"{res[0]:.16f}{delim_char}{res[1]:.16f}"

            if dim == 3:
                out_coords += f"{delim_char}{res[2]:.4f}"

        else:
            res = trans.etrs89_to_st70(c1, c2, h=c3, source=source)
            out_coords = f"{res[0]:.4f}{delim_char}{res[1]:.4f}"
            if dim == 3:
                out_coords += f"{delim_char}{res[2]:.4f}"

        if point_id is not None:
            return f"{point_id}{delim_char}{out_coords}\n"
        return f"{out_coords}\n"

    except exc.OutOfGridErr:
        return line + " # Out of grid\n"
    except exc.NoDataGridErr:
        return line + " # No data on grid\n"
    except exc.WrongDimensionErr:
        return line + " # Wrong dimension\n"
    except exc.InvalidCoordinateErr:
        return line + " # Invalid coordinate data\n"


@router.post("/file")
async def transform_file(
    request: Request,
    file: UploadFile = File(..., description="CSV or TXT file containing coordinates"),
    op: OperationEnum = Form(..., description="Transformation direction"),
    unit: FileAngleUnitEnum = Form(
        default=FileAngleUnitEnum.RADIANS, description="Angular unit for ETRS89 coordinates"
    ),
    delimiter: str = Form(
        default="auto", description="Delimiter: auto, comma, semicolon, tab, or space"
    ),
    skip_header: bool = Form(
        default=False, description="Whether to ignore the first header line without transformation"
    ),
    trans: TransRO = Depends(get_trans_engine),
    source: int = Depends(get_client_source),
    settings: Settings = Depends(get_settings),
):
    """Transforms a delimited file (CSV, TXT) streaming results back using chunked HTTP transfer.

    Processes points with negligible memory overhead up to the configured limit.
    """
    upload_name = file.filename or "unknown"
    delim_key = delimiter.strip().lower()
    if delim_key not in ALLOWED_DELIMITERS:
        logger.warning(
            "File upload '%s' rejected: invalid delimiter '%s'",
            upload_name,
            delimiter,
        )
        raise HTTPException(
            status_code=422,
            detail=f"Invalid delimiter '{delimiter}'. Allowed values: auto, comma, semicolon, tab, space.",
        )

    max_bytes = settings.max_file_upload_bytes
    content_length = request.headers.get("content-length")
    if content_length and content_length.isdigit() and int(content_length) > max_bytes:
        logger.warning(
            "File upload '%s' rejected: content length %s exceeds max %d bytes",
            upload_name,
            content_length,
            max_bytes,
        )
        raise HTTPException(
            status_code=413,
            detail=f"File upload exceeds maximum allowed size of {max_bytes:,} bytes.",
        )

    if file.size is not None and file.size > max_bytes:
        logger.warning(
            "File upload '%s' rejected: file size %s exceeds max %d bytes",
            upload_name,
            file.size,
            max_bytes,
        )
        raise HTTPException(
            status_code=413,
            detail=f"File upload exceeds maximum allowed size of {max_bytes:,} bytes.",
        )

    chosen_delim = DELIMITER_MAP.get(delim_key)
    buffered_collector = BufferedTelemetryCollector(trans._collector, chunk_size=500)
    file_trans = trans.with_collector(buffered_collector)

    async def line_generator() -> AsyncGenerator[bytes, None]:
        nonlocal chosen_delim
        row_count = 0
        total_bytes_read = 0
        max_rows = settings.max_file_points
        logged_errors = 0
        max_logged_errors = 5
        header_skipped = False

        # Read line by line from uploaded stream
        buffer = ""
        try:
            while True:
                chunk = await file.read(64 * 1024)  # 64 KB chunks
                if not chunk:
                    break
                total_bytes_read += len(chunk)
                if total_bytes_read > max_bytes:
                    yield f"# Limit exceeded: maximum upload size of {max_bytes:,} bytes exceeded\n".encode(
                        "utf-8"
                    )
                    return

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

                        if skip_header and not header_skipped:
                            header_skipped = True
                            yield (line + "\n").encode("utf-8")
                            continue

                        delim_char = chosen_delim if chosen_delim != " " else " "
                        row_count += 1
                        if row_count > max_rows:
                            yield f"# Limit exceeded: maximum {max_rows:,} points allowed\n".encode(
                                "utf-8"
                            )
                            return

                        out = _transform_line(line, delim_char, op, unit, file_trans, source=source)
                        if " #" in out:
                            logged_errors += 1
                            if logged_errors <= max_logged_errors:
                                reason = out.split(" #")[-1].strip()
                                logger.warning(
                                    "File upload '%s' issue at line %d: %s",
                                    upload_name,
                                    row_count,
                                    reason,
                                )
                            elif logged_errors == max_logged_errors + 1:
                                logger.warning(
                                    "File upload '%s': additional line error warnings suppressed (exceeded %d errors).",
                                    upload_name,
                                    max_logged_errors,
                                )
                        yield out.encode("utf-8")

            # Process any remaining line in buffer
            if buffer.strip():
                line = buffer.strip()
                if line.startswith("#"):
                    yield (line + "\n").encode("utf-8")
                else:
                    if chosen_delim is None:
                        chosen_delim = _detect_delimiter(line)
                    if skip_header and not header_skipped:
                        header_skipped = True
                        yield (line + "\n").encode("utf-8")
                    else:
                        delim_char = chosen_delim if chosen_delim != " " else " "
                        row_count += 1
                        if row_count > max_rows:
                            yield f"# Limit exceeded: maximum {max_rows:,} points allowed\n".encode(
                                "utf-8"
                            )
                            return
                        out = _transform_line(line, delim_char, op, unit, file_trans, source=source)
                        yield out.encode("utf-8")
        finally:
            buffered_collector.flush_all()

    return StreamingResponse(
        line_generator(),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="transformed_coordinates.csv"'},
    )
