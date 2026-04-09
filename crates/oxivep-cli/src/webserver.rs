use anyhow::Result;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;

use crate::pipeline::AnnotationContext;

const INDEX_HTML: &str = include_str!("../../../web/index.html");

pub fn run_server(port: u16, gff3: Option<String>, fasta: Option<String>) -> Result<()> {
    let mut ctx = AnnotationContext::new(
        gff3.as_deref(),
        fasta.as_deref(),
        5000,
    )?;

    let addr = format!("0.0.0.0:{}", port);
    let listener = TcpListener::bind(&addr)?;

    eprintln!("OxiVEP web interface running at http://localhost:{}", port);
    eprintln!("Press Ctrl+C to stop.");

    for stream in listener.incoming() {
        match stream {
            Ok(mut stream) => {
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
