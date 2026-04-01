use anyhow::Result;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;

const INDEX_HTML: &str = include_str!("../../../web/index.html");

pub fn run_server(port: u16, _gff3: Option<String>, _fasta: Option<String>) -> Result<()> {
    let addr = format!("0.0.0.0:{}", port);
    let listener = TcpListener::bind(&addr)?;

    eprintln!("OxiVEP web interface running at http://localhost:{}", port);
    eprintln!("Press Ctrl+C to stop.");

    for stream in listener.incoming() {
        match stream {
            Ok(mut stream) => {
                let mut reader = BufReader::new(stream.try_clone()?);
                let mut request_line = String::new();
                reader.read_line(&mut request_line)?;

                // Read remaining headers
                loop {
                    let mut line = String::new();
                    reader.read_line(&mut line)?;
                    if line.trim().is_empty() {
                        break;
                    }
                }

                let path = request_line
                    .split_whitespace()
                    .nth(1)
                    .unwrap_or("/");

                match path {
                    "/" | "/index.html" => {
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
                    "/api/annotate" => {
                        // Read POST body for future API endpoint
                        let response = "HTTP/1.1 200 OK\r\n\
                            Content-Type: application/json\r\n\
                            Connection: close\r\n\
                            \r\n\
                            {\"status\":\"ok\",\"message\":\"Use the client-side annotation engine in the web UI\"}";
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
            }
            Err(e) => {
                eprintln!("Connection error: {}", e);
            }
        }
    }

    Ok(())
}
