"""Pydantic schema for uniform API error responses."""

from typing import Optional

from pydantic import BaseModel, Field


class ErrorResponse(BaseModel):
    """Uniform error response structure for REST API endpoints."""

    detail: str = Field(..., description="Human-readable error description")
    code: Optional[str] = Field(None, description="Machine-readable error classification code")

    model_config = {
        "json_schema_extra": {
            "example": {
                "detail": "Requested coordinate point is outside the Romanian national grid boundary.",
                "code": "OUT_OF_GRID",
            }
        }
    }
