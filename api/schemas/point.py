"""Pydantic schemas for single point transformations.
"""
from enum import Enum
from typing import List, Optional
from pydantic import BaseModel, Field


class OperationEnum(str, Enum):
    STEREO70_TO_ETRS89 = "Stereo70ToETRS89"
    ETRS89_TO_STEREO70 = "ETRS89ToStereo70"


class AngleUnitEnum(str, Enum):
    RADIANS = "radians"
    DEGREES = "degrees"


class PointRequest(BaseModel):
    op: OperationEnum = Field(..., description="Transformation direction: Stereo70ToETRS89 or ETRS89ToStereo70")
    coos: List[float] = Field(
        ...,
        min_length=2,
        max_length=3,
        description="Coordinates: [Northing, Easting[, H]] for Stereo70 or [Latitude, Longitude[, h]] for ETRS89"
    )
    unit: AngleUnitEnum = Field(
        default=AngleUnitEnum.RADIANS,
        description="Angular unit for ETRS89 coordinates: 'radians' (default) or 'degrees'"
    )

    model_config = {
        "json_schema_extra": {
            "example": {
                "op": "Stereo70ToETRS89",
                "coos": [500000.0, 500000.0, 100.0],
                "unit": "degrees"
            }
        }
    }


class PointResponse(BaseModel):
    coos: List[float] = Field(..., description="Transformed coordinates [N, E[, H]] or [Lat, Lon[, h]]")
    warning: Optional[str] = Field(default=None, description="Warning message if an issue occurred, null on success")

    model_config = {
        "json_schema_extra": {
            "example": {
                "coos": [46.0, 25.0, 139.6084],
                "warning": None
            }
        }
    }
