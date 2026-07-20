use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::Json;

#[derive(Debug)]
pub enum AppError {
    Internal(anyhow::Error),
    BadRequest(String),
}

impl IntoResponse for AppError {
    fn into_response(self) -> Response {
        let (status, message) = match self {
            AppError::Internal(e) => {
                // Log the full error chain (which can include file paths and
                // parse internals) server-side only; clients get a generic
                // message so internal details aren't leaked in API responses.
                tracing::error!("Internal error: {:?}", e);
                (
                    StatusCode::INTERNAL_SERVER_ERROR,
                    "Internal server error".to_string(),
                )
            }
            AppError::BadRequest(msg) => (StatusCode::BAD_REQUEST, msg),
        };
        (status, Json(serde_json::json!({"error": message}))).into_response()
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
