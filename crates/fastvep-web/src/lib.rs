//! Library surface for the fastVEP web server.
//!
//! The binary in `main.rs` is a thin wrapper over this: it parses flags, loads
//! the annotation context, and hands both to [`build_router`]. The routing
//! table lives here so tests can drive it through `tower::ServiceExt::oneshot`
//! without binding a socket, which is what lets the JSON contract documented in
//! `docs/API.md` be asserted rather than described.

pub mod context;
pub mod errors;
pub mod handlers;

use axum::extract::DefaultBodyLimit;
use axum::http::{header, Method};
use axum::routing::{get, post};
use axum::Router;
use tower::limit::ConcurrencyLimitLayer;
use tower::ServiceBuilder;
use tower_http::cors::{Any, CorsLayer};
use tower_http::trace::TraceLayer;

use crate::handlers::AppState;

/// Build the full application router, including the middleware stack.
pub fn build_router(state: AppState, max_body_size: usize, max_concurrent: usize) -> Router {
    Router::new()
        .route("/", get(handlers::index_html))
        .route("/index.html", get(handlers::index_html))
        .route("/assets/logo.png", get(handlers::logo_png))
        .route("/favicon.ico", get(handlers::logo_png))
        .route("/apple-touch-icon.png", get(handlers::logo_png))
        .route("/api/status", get(handlers::status))
        .route("/api/genomes", get(handlers::list_genomes))
        .route("/api/load-genome", post(handlers::load_genome))
        .route("/api/annotate", post(handlers::annotate))
        .route("/api/upload-gff3", post(handlers::upload_gff3))
        .with_state(state)
        .layer(DefaultBodyLimit::max(max_body_size))
        .layer(
            ServiceBuilder::new()
                .layer(TraceLayer::new_for_http())
                .layer(
                    // Origin stays open (this is a public annotation API by
                    // design, per DEPLOYMENT.md), but methods/headers are
                    // scoped to what the routes above actually use instead
                    // of blanket `Any` on every axis.
                    CorsLayer::new()
                        .allow_origin(Any)
                        .allow_methods([Method::GET, Method::POST])
                        .allow_headers([header::CONTENT_TYPE]),
                )
                .layer(ConcurrencyLimitLayer::new(max_concurrent)),
        )
}
