"""Pydantic schemas for anonymous, zero-PII feedback and discrepancy reporting.
"""
from typing import Optional, Literal
from pydantic import BaseModel, Field


FeedbackReportType = Literal[
    "discrepancy",
    "grid_update",
    "bug",
    "feature",
    "api_down",
]


class FeedbackSubmitRequest(BaseModel):
    """Payload for submitting anonymous feedback."""
    report_type: FeedbackReportType = Field(
        ...,
        description="Category of the report (discrepancy, grid_update, bug, feature, api_down)",
    )
    description: str = Field(
        ...,
        min_length=10,
        max_length=5000,
        description="Detailed description or filled template for the feedback",
    )
    client_timestamp: Optional[str] = Field(
        None,
        description="Client-side ISO timestamp when the report was composed",
    )
    app_version: Optional[str] = Field(
        None,
        description="Client application version identifier",
    )
    active_grid: Optional[str] = Field(
        None,
        description="Active grid identifier or name reported by client",
    )
    honeypot: Optional[str] = Field(
        None,
        description="Hidden honeypot field. Must be empty for legitimate human submissions",
    )


class FeedbackSubmitResponse(BaseModel):
    """Response returned upon feedback submission."""
    status: Literal["success", "logged"] = Field(
        ...,
        description="Result status ('success' if emailed, 'logged' if recorded locally)",
    )
    delivered_email: bool = Field(
        ...,
        description="Whether the report was dispatched via configured SMTP email",
    )
    message: str = Field(
        ...,
        description="Human-readable confirmation message",
    )
