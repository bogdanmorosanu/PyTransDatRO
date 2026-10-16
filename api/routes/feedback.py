"""Feedback submission API endpoint."""

from fastapi import APIRouter, Depends, Request

from api.config import Settings, get_settings
from api.schemas.feedback import FeedbackSubmitRequest, FeedbackSubmitResponse
from api.services.mailer import dispatch_feedback

router = APIRouter(prefix="/api/v1", tags=["Feedback"])


@router.post(
    "/feedback",
    response_model=FeedbackSubmitResponse,
    summary="Submit anonymous user feedback",
    description="Submits an anonymous feedback, discrepancy, or bug report with zero PII.",
)
async def submit_feedback(
    payload: FeedbackSubmitRequest,
    request: Request,
    settings: Settings = Depends(get_settings),
):
    # 1. Zero-friction bot check (Honeypot)
    if payload.honeypot and payload.honeypot.strip():
        # Silently accept without action to misdirect automated spambots
        return FeedbackSubmitResponse(
            status="success",
            delivered_email=False,
            message="Feedback received. Thank you!",
        )

    # 2. Extract client IP safely (respecting potential reverse proxies)
    forwarded_for = request.headers.get("x-forwarded-for")
    if forwarded_for:
        client_ip = forwarded_for.split(",")[0].strip()
    elif request.client and request.client.host:
        client_ip = request.client.host
    else:
        client_ip = "127.0.0.1"

    # 3. Dispatch to mailer / logger service
    delivered_email, msg = dispatch_feedback(payload.model_dump(), client_ip, settings=settings)

    return FeedbackSubmitResponse(
        status="success" if delivered_email else "logged",
        delivered_email=delivered_email,
        message=(
            "Thank you! Your feedback has been recorded."
            if not delivered_email
            else "Thank you! Your feedback has been sent successfully."
        ),
    )
