"""Streaming delimited file transformation endpoint.
Supports up to 100,000 coordinates in CSV/TXT formats.
"""
import math
from typing import AsyncGenerator
from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import StreamingResponse
from api.dependencies import get_trans_engine
from api.schemas.point import OperationEnum, AngleUnitEnum
from pytransdatro.trans_ro import TransRO
import pytransdatro.exceptions as exc

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


@router.post("/file")
async def transform_file(
    file: UploadFile = File(..., description="CSV or TXT file containing coordinates"),
    op: OperationEnum = Form(..., description="Transformation direction"),
    unit: AngleUnitEnum = Form(default=AngleUnitEnum.RADIANS, description="Angular unit for ETRS89 coordinates"),
    delimiter: str = Form(default="auto", description="Delimiter: auto, comma, semicolon, tab, or space"),
    trans: TransRO = Depends(get_trans_engine)
):
    """Transforms a delimited file (CSV, TXT) streaming results back using chunked HTTP transfer.
    
    Processes up to 100,000 points with negligible memory overhead.
    """
    is_degrees = (unit == AngleUnitEnum.DEGREES)
    chosen_delim = DELIMITER_MAP.get(delimiter.lower())

    async def line_generator() -> AsyncGenerator[bytes, None]:
        nonlocal chosen_delim
        row_count = 0
        max_rows = 100000

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

                    # Split tokens
                    tokens = [t.strip() for t in (line.split() if chosen_delim == " " else line.split(chosen_delim))]
                    tokens = [t for t in tokens if t]

                    if len(tokens) < 2:
                        yield (line + " # Wrong dimension\n").encode("utf-8")
                        continue

                    try:
                        coords = [float(t) for t in tokens[:3]]
                    except ValueError:
                        yield (line + " # Invalid coordinate data\n").encode("utf-8")
                        continue

                    dim = len(coords)
                    row_count += 1
                    if row_count > max_rows:
                        yield b"# Limit exceeded: maximum 100,000 points allowed\n"
                        return

                    delim_char = chosen_delim if chosen_delim != " " else " "

                    try:
                        if op == OperationEnum.STEREO70_TO_ETRS89:
                            n, e = coords[0], coords[1]
                            z = coords[2] if dim >= 3 else None
                            res = trans.st70_to_etrs89(n, e, z=z, source=2)

                            if is_degrees:
                                lat_deg = math.degrees(res[0])
                                lon_deg = math.degrees(res[1])
                                out_line = f"{lat_deg:.10f}{delim_char}{lon_deg:.10f}"
                                if dim >= 3:
                                    out_line += f"{delim_char}{res[2]:.4f}"
                            else:
                                out_line = f"{res[0]:.16f}{delim_char}{res[1]:.16f}"
                                if dim >= 3:
                                    out_line += f"{delim_char}{res[2]:.4f}"

                            yield (out_line + "\n").encode("utf-8")

                        elif op == OperationEnum.ETRS89_TO_STEREO70:
                            lat, lon = coords[0], coords[1]
                            if is_degrees:
                                lat = math.radians(lat)
                                lon = math.radians(lon)

                            h = coords[2] if dim >= 3 else None
                            res = trans.etrs89_to_st70(lat, lon, h=h, source=2)

                            out_line = f"{res[0]:.4f}{delim_char}{res[1]:.4f}"
                            if dim >= 3:
                                out_line += f"{delim_char}{res[2]:.4f}"

                            yield (out_line + "\n").encode("utf-8")

                    except exc.OutOfGridErr:
                        yield (line + " # Out of grid\n").encode("utf-8")
                    except exc.NoDataGridErr:
                        yield (line + " # No data on grid\n").encode("utf-8")

        # Process any remaining line in buffer
        if buffer.strip():
            line = buffer.strip()
            if not line.startswith("#"):
                tokens = [t.strip() for t in (line.split() if chosen_delim == " " else line.split(chosen_delim))]
                tokens = [t for t in tokens if t]
                if len(tokens) >= 2:
                    try:
                        coords = [float(t) for t in tokens[:3]]
                        dim = len(coords)
                        delim_char = chosen_delim or ","
                        if op == OperationEnum.STEREO70_TO_ETRS89:
                            n, e = coords[0], coords[1]
                            z = coords[2] if dim >= 3 else None
                            res = trans.st70_to_etrs89(n, e, z=z, source=2)
                            if is_degrees:
                                out_line = f"{math.degrees(res[0]):.10f}{delim_char}{math.degrees(res[1]):.10f}"
                            else:
                                out_line = f"{res[0]:.16f}{delim_char}{res[1]:.16f}"
                            if dim >= 3:
                                out_line += f"{delim_char}{res[2]:.4f}"
                            yield (out_line + "\n").encode("utf-8")
                        else:
                            lat, lon = coords[0], coords[1]
                            if is_degrees:
                                lat, lon = math.radians(lat), math.radians(lon)
                            h = coords[2] if dim >= 3 else None
                            res = trans.etrs89_to_st70(lat, lon, h=h, source=2)
                            out_line = f"{res[0]:.4f}{delim_char}{res[1]:.4f}"
                            if dim >= 3:
                                out_line += f"{delim_char}{res[2]:.4f}"
                            yield (out_line + "\n").encode("utf-8")
                    except (exc.OutOfGridErr, exc.NoDataGridErr, ValueError):
                        yield (line + "\n").encode("utf-8")

    return StreamingResponse(
        line_generator(),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="transformed_coordinates.csv"'}
    )
