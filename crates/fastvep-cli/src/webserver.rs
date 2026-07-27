use anyhow::Result;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;
use std::time::Duration;

use fastvep_annotate::AnnotationContext;

const INDEX_HTML: &str = include_str!("../../../web/index.html");

/// Upper bound on an untrusted client's `Content-Length` header. Requests
/// claiming a larger body are rejected before any allocation is attempted.
/// Generous enough for real annotation payloads, bounded enough that a
/// malicious/corrupted header can't force a huge `Vec::with_capacity`.
const MAX_CONTENT_LENGTH: usize = 100 * 1024 * 1024; // 100 MiB

/// Read timeout applied to every accepted connection before any reads
/// happen, so a connected-but-silent client can't hang the (single-
/// threaded, blocking) server indefinitely.
const READ_TIMEOUT: Duration = Duration::from_secs(30);

pub fn run_server(port: u16, gff3: Option<String>, fasta: Option<String>) -> Result<()> {
    let mut ctx = AnnotationContext::new(
        gff3.as_deref(),
        fasta.as_deref(),
        None, // no SA directory for embedded web server
        5000,
    )?;

    let addr = format!("0.0.0.0:{}", port);
    let listener = TcpListener::bind(&addr)?;

    eprintln!("fastVEP web interface running at http://localhost:{}", port);
    eprintln!("Press Ctrl+C to stop.");

    for stream in listener.incoming() {
        match stream {
            Ok(mut stream) => {
                // A connected-but-silent (or slow-drip) client must not be
                // able to hang this single-threaded server forever.
                if let Err(e) = stream.set_read_timeout(Some(READ_TIMEOUT)) {
                    eprintln!("Failed to set read timeout: {}", e);
                }
                if let Err(e) = handle_request(&mut stream, &mut ctx) {
                    eprintln!("Request error: {}", e);
                }
            }
            Err(e) => {
                eprintln!("Connection error: {}", e);
            }
        }
    }

    Ok(())
}

fn send_json(stream: &mut std::net::TcpStream, status: u16, body: &str) -> Result<()> {
    let status_text = match status {
        200 => "OK",
        400 => "Bad Request",
        413 => "Payload Too Large",
        500 => "Internal Server Error",
        _ => "OK",
    };
    let response = format!(
        "HTTP/1.1 {} {}\r\n\
         Content-Type: application/json\r\n\
         Content-Length: {}\r\n\
         Access-Control-Allow-Origin: *\r\n\
         Connection: close\r\n\
         \r\n\
         {}",
        status, status_text, body.len(), body
    );
    stream.write_all(response.as_bytes())?;
    Ok(())
}

fn handle_request(stream: &mut std::net::TcpStream, ctx: &mut AnnotationContext) -> Result<()> {
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut request_line = String::new();
    reader.read_line(&mut request_line)?;

    let parts: Vec<&str> = request_line.split_whitespace().collect();
    let method = parts.first().unwrap_or(&"GET");
    let path = parts.get(1).unwrap_or(&"/");

    // Read headers, extract Content-Length
    let mut content_length: usize = 0;
    loop {
        let mut line = String::new();
        reader.read_line(&mut line)?;
        if line.trim().is_empty() {
            break;
        }
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("content-length:") {
            content_length = lower.trim_start_matches("content-length:").trim().parse().unwrap_or(0);
        }
    }

    // Reject an oversized Content-Length up front, before any request
    // handler attempts `vec![0u8; content_length]` — an untrusted client
    // could otherwise claim an enormous body and force a huge allocation.
    if content_length > MAX_CONTENT_LENGTH {
        send_json(
            stream,
            413,
            r#"{"error":"Request body exceeds maximum allowed size"}"#,
        )?;
        return Ok(());
    }

    match (*method, *path) {
        ("GET", "/" | "/index.html") => {
            let response = format!(
                "HTTP/1.1 200 OK\r\n\
                 Content-Type: text/html; charset=utf-8\r\n\
                 Content-Length: {}\r\n\
                 Connection: close\r\n\
                 \r\n\
                 {}",
                INDEX_HTML.len(),
                INDEX_HTML
            );
            stream.write_all(response.as_bytes())?;
        }
        ("GET", "/api/status") => {
            let tr_count = ctx.transcript_provider.transcript_count();
            let status_json = serde_json::json!({
                "status": "ok",
                "backend": true,
                "transcripts": tr_count,
                "gff3_source": ctx.gff3_source,
                "has_fasta": ctx.seq_provider.is_some(),
            });
            send_json(stream, 200, &serde_json::to_string(&status_json)?)?;
        }
        ("POST", "/api/upload-gff3") => {
            let mut body = vec![0u8; content_length];
            reader.read_exact(&mut body)?;
            let gff3_text = String::from_utf8_lossy(&body);

            let start = std::time::Instant::now();
            match ctx.update_gff3_text(&gff3_text) {
                Ok((genes, transcripts)) => {
                    let elapsed = start.elapsed().as_millis();
                    let resp = serde_json::json!({
                        "genes": genes,
                        "transcripts": transcripts,
                        "time_ms": elapsed,
                    });
                    send_json(stream, 200, &serde_json::to_string(&resp)?)?;
                }
                Err(e) => {
                    let resp = serde_json::json!({"error": format!("{}", e)});
                    send_json(stream, 500, &serde_json::to_string(&resp)?)?;
                }
            }
        }
        ("POST", "/api/annotate") => {
            let mut body = vec![0u8; content_length];
            reader.read_exact(&mut body)?;
            let body_str = String::from_utf8_lossy(&body);

            let request: serde_json::Value = serde_json::from_str(&body_str)
                .unwrap_or(serde_json::json!({}));
            let vcf_text = request["vcf"].as_str().unwrap_or("");
            let pick = request["pick"].as_bool().unwrap_or(false);

            if vcf_text.is_empty() {
                send_json(stream, 400, r#"{"error":"No VCF data provided"}"#)?;
            } else {
                let start = std::time::Instant::now();
                match ctx.annotate_vcf_text(vcf_text, pick) {
                    Ok(results) => {
                        let elapsed = start.elapsed().as_millis();
                        let resp = serde_json::json!({
                            "results": results,
                            "count": results.len(),
                            "time_ms": elapsed,
                        });
                        send_json(stream, 200, &serde_json::to_string(&resp)?)?;
                    }
                    Err(e) => {
                        let resp = serde_json::json!({"error": format!("{}", e)});
                        send_json(stream, 500, &serde_json::to_string(&resp)?)?;
                    }
                }
            }
        }
        ("OPTIONS", _) => {
            let response = "HTTP/1.1 204 No Content\r\n\
                 Access-Control-Allow-Origin: *\r\n\
                 Access-Control-Allow-Methods: POST, GET, OPTIONS\r\n\
                 Access-Control-Allow-Headers: Content-Type\r\n\
                 Connection: close\r\n\
                 \r\n";
            stream.write_all(response.as_bytes())?;
        }
        _ => {
            let response = "HTTP/1.1 404 Not Found\r\n\
                 Content-Type: text/plain\r\n\
                 Connection: close\r\n\
                 \r\n\
                 404 Not Found";
            stream.write_all(response.as_bytes())?;
        }
    }
    stream.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::TcpStream;

    /// Spawn a one-shot server thread that accepts a single connection and
    /// runs `handle_request` on it, returning the thread's join handle and
    /// the address to connect to.
    fn spawn_one_shot_server() -> (std::thread::JoinHandle<()>, std::net::SocketAddr) {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind ephemeral port");
        let addr = listener.local_addr().unwrap();
        let handle = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("accept connection");
            let mut ctx = AnnotationContext::new(None, None, None, 5000)
                .expect("build empty AnnotationContext");
            handle_request(&mut stream, &mut ctx).expect("handle_request should not error");
        });
        (handle, addr)
    }

    #[test]
    fn content_length_over_cap_is_rejected_without_allocating() {
        // A Content-Length far beyond MAX_CONTENT_LENGTH must be rejected
        // immediately with a clean error response — and critically, the
        // server must never attempt to read (or allocate for) that many
        // body bytes. We only send the headers and no body at all; if the
        // server tried `vec![0u8; content_length]` + `read_exact`, this test
        // would hang (and eventually fail on the client's read timeout)
        // instead of getting an immediate response.
        let (server, addr) = spawn_one_shot_server();
        let mut client = TcpStream::connect(addr).expect("connect to test server");
        client
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();

        let huge = MAX_CONTENT_LENGTH + 1;
        let request = format!(
            "POST /api/annotate HTTP/1.1\r\nContent-Length: {}\r\n\r\n",
            huge
        );
        client.write_all(request.as_bytes()).unwrap();

        let mut response = String::new();
        client
            .read_to_string(&mut response)
            .expect("server should respond promptly instead of hanging");

        server.join().expect("server thread should not panic");

        assert!(
            response.starts_with("HTTP/1.1 413"),
            "expected a 413 response, got: {}",
            response
        );
        assert!(response.contains("Request body exceeds maximum allowed size"));
    }

    #[test]
    fn normal_request_parsing_still_works() {
        // Baseline regression: a well-formed GET request with no body must
        // still be handled exactly as before these hardening changes.
        let (server, addr) = spawn_one_shot_server();
        let mut client = TcpStream::connect(addr).expect("connect to test server");
        client
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();

        client
            .write_all(b"GET /api/status HTTP/1.1\r\n\r\n")
            .unwrap();

        let mut response = String::new();
        client.read_to_string(&mut response).unwrap();

        server.join().expect("server thread should not panic");

        assert!(
            response.starts_with("HTTP/1.1 200"),
            "expected a 200 response, got: {}",
            response
        );
        assert!(response.contains("\"status\":\"ok\""));
    }
}
