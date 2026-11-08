use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::Json;

#[derive(Debug)]
pub enum AppError {
    Internal(anyhow::Error),
    BadRequest(String),
    /// The request was understood and refused by configuration - the caller
    /// cannot fix it by changing the request, so it is not a 400.
    Forbidden(String),
    /// The request is well formed and the server could produce an answer, but
    /// not one it can stand behind: ACMG classification was asked of a genome
    /// that lacks sources the classifier reads. The caller can fix it only by
    /// choosing another genome or by not asking for ACMG, which is why the body
    /// carries the verdict (`acmg`) beside the message.
    AcmgUnavailable {
        message: String,
        acmg: serde_json::Value,
    },
}

impl IntoResponse for AppError {
    fn into_response(self) -> Response {
        let (status, body) = match self {
            AppError::Internal(e) => {
                // Log the full error chain (which can include file paths and
                // parse internals) server-side only; clients get a generic
                // message so internal details aren't leaked in API responses.
                tracing::error!("Internal error: {:?}", e);
                (
                    StatusCode::INTERNAL_SERVER_ERROR,
                    serde_json::json!({"error": "Internal server error"}),
                )
            }
            AppError::BadRequest(msg) => {
                (StatusCode::BAD_REQUEST, serde_json::json!({"error": msg}))
            }
            AppError::Forbidden(msg) => (StatusCode::FORBIDDEN, serde_json::json!({"error": msg})),
            AppError::AcmgUnavailable { message, acmg } => (
                StatusCode::UNPROCESSABLE_ENTITY,
                serde_json::json!({"error": message, "acmg": acmg}),
            ),
        };
        (status, Json(body)).into_response()
    }
}

impl From<anyhow::Error> for AppError {
    fn from(e: anyhow::Error) -> Self {
        AppError::Internal(e)
    }
}

impl From<tokio::task::JoinError> for AppError {
    fn from(e: tokio::task::JoinError) -> Self {
        AppError::Internal(anyhow::anyhow!("Task join error: {}", e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    async fn body_json(resp: Response) -> serde_json::Value {
        let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
            .await
            .unwrap();
        serde_json::from_slice(&bytes).unwrap()
    }

    #[tokio::test]
    async fn internal_error_does_not_leak_details_to_client() {
        // Regression test: AppError::Internal used to forward the raw
        // anyhow error chain (which can contain file paths, parse internals,
        // lock-poison messages) directly into the JSON response body.
        let err = AppError::Internal(anyhow::anyhow!(
            "failed to open /secret/internal/path.gff3: permission denied"
        ));
        let resp = err.into_response();
        assert_eq!(resp.status(), StatusCode::INTERNAL_SERVER_ERROR);
        let body = body_json(resp).await;
        let message = body["error"].as_str().unwrap();
        assert_eq!(message, "Internal server error");
        assert!(!message.contains("/secret/internal/path.gff3"));
    }

    #[tokio::test]
    async fn an_unavailable_classification_is_a_422_that_carries_the_verdict() {
        let err = AppError::AcmgUnavailable {
            message: "not available".into(),
            acmg: serde_json::json!({"ready": false}),
        };
        let resp = err.into_response();
        assert_eq!(resp.status(), StatusCode::UNPROCESSABLE_ENTITY);
        let body = body_json(resp).await;
        assert_eq!(body["error"], "not available");
        assert_eq!(body["acmg"]["ready"], false);
    }

    #[tokio::test]
    async fn bad_request_message_is_passed_through() {
        // BadRequest is for caller-facing validation messages, so unlike
        // Internal, its text should reach the client unchanged.
        let err = AppError::BadRequest("No VCF data provided".to_string());
        let resp = err.into_response();
        assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
        let body = body_json(resp).await;
        assert_eq!(body["error"], "No VCF data provided");
    }
}
