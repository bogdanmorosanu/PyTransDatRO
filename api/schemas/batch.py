"""Pydantic schemas for in-memory batch coordinate transformations.
"""
from typing import List
from pydantic import BaseModel, Field
from api.schemas.point import OperationEnum, AngleUnitEnum, PointResponse


class BatchRequest(BaseModel):
    op: OperationEnum = Field(..., description="Transformation direction: Stereo70ToETRS89 or ETRS89ToStereo70")
    points: List[List[float]] = Field(
        ...,
        min_length=1,
        max_length=20000,
        description="List of 2D or 3D coordinate arrays"
    )
    unit: AngleUnitEnum = Field(
        default=AngleUnitEnum.RADIANS,
        description="Angular unit for ETRS89 coordinates: 'radians' (default) or 'degrees'"
    )

    model_config = {
        "json_schema_extra": {
            "example": {
                "op": "Stereo70ToETRS89",
                "points": [
                    [500000.0, 500000.0, 100.0],
                    [505000.0, 505000.0, 120.0]
                ],
                "unit": "degrees"
            }
        }
    }


class BatchResponse(BaseModel):
    results: List[PointResponse] = Field(..., description="Ordered list of transformation results matching input index")
    count: int = Field(..., description="Total count of points processed")
