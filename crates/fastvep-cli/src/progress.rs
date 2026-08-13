use std::io::{self, Read};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};
use std::time::Instant;

const REPORT_INTERVAL_SECS: f64 = 10.0;

/// Wraps any `Read` and counts bytes consumed, accessible via a shared `Arc<AtomicU64>`.
pub struct CountingReader<R: Read> {
    pub inner: R,
    pub bytes: Arc<AtomicU64>,
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

pub struct ProgressMeter {
    enabled: bool,
    start: Instant,
    last_report: Instant,
    pub records: u64,
    header_printed: bool,
    total_bytes: u64,
    bytes_read: Option<Arc<AtomicU64>>,
}

impl ProgressMeter {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            start: Instant::now(),
            last_report: Instant::now(),
            records: 0,
            header_printed: false,
            total_bytes: 0,
            bytes_read: None,
        }
    }

    /// Like `new`, but also tracks byte-based progress for ETA output.
    /// `total_bytes` is the compressed file size; `bytes_read` is a counter
    /// updated by a `CountingReader` wrapping the underlying file.
    pub fn with_progress(enabled: bool, total_bytes: u64, bytes_read: Arc<AtomicU64>) -> Self {
        Self {
            enabled,
            start: Instant::now(),
            last_report: Instant::now(),
            records: 0,
            header_printed: false,
            total_bytes,
            bytes_read: Some(bytes_read),
        }
    }

    pub fn update(&mut self) {
        self.update_n(1);
    }

    pub fn update_n(&mut self, n: u64) {
        if !self.enabled {
            return;
        }
        self.records += n;
        if self.last_report.elapsed().as_secs_f64() >= REPORT_INTERVAL_SECS {
            self.print_row();
        }
    }

    fn has_eta(&self) -> bool {
        self.total_bytes > 0 && self.bytes_read.is_some()
    }

    fn print_row(&mut self) {
        if !self.header_printed {
            if self.has_eta() {
                eprintln!(
                    "INFO  ProgressMeter -  {:>13}   {:>20}   {:>14}   {:>7}   {:>9}",
                    "Elapsed (min)", "Records Processed", "Records/Min", "% Done", "ETA (min)"
                );
            } else {
                eprintln!(
                    "INFO  ProgressMeter -  {:>13}   {:>20}   {:>14}",
                    "Elapsed (min)", "Records Processed", "Records/Min"
                );
            }
            self.header_printed = true;
        }

        let elapsed_min = self.start.elapsed().as_secs_f64() / 60.0;
        let rpm = if elapsed_min > 0.0 {
            self.records as f64 / elapsed_min
        } else {
            0.0
        };

        if let Some(ref counter) = self.bytes_read {
            let done_frac = counter.load(Ordering::Relaxed) as f64 / self.total_bytes as f64;
            let pct = done_frac * 100.0;
            let eta_min = if done_frac > 0.0 {
                elapsed_min * (1.0 / done_frac - 1.0)
            } else {
                0.0
            };
            eprintln!(
                "INFO  ProgressMeter -  {:>13.1}   {:>20}   {:>14.1}   {:>6.1}%   {:>9.1}",
                elapsed_min,
                fmt_count(self.records),
                rpm,
                pct,
                eta_min,
            );
        } else {
            eprintln!(
                "INFO  ProgressMeter -  {:>13.1}   {:>20}   {:>14.1}",
                elapsed_min,
                fmt_count(self.records),
                rpm
            );
        }
        self.last_report = Instant::now();
    }

    pub fn finish(mut self) {
        if !self.enabled {
            return;
        }
        // The gz decoder can stop a few trailing compressed bytes short of EOF,
        // so the final byte-based row would otherwise read e.g. 99.9% / ETA>0.
        // Pin the counter to the full size so the closing line reads 100% / 0.
        if let Some(ref counter) = self.bytes_read {
            counter.store(self.total_bytes, Ordering::Relaxed);
        }
        self.print_row();
        let elapsed_min = self.start.elapsed().as_secs_f64() / 60.0;
        eprintln!(
            "INFO  ProgressMeter - Done. Processed {} records in {:.1} min.",
            fmt_count(self.records),
            elapsed_min
        );
    }
}

/// Human-readable byte count with a binary-prefix unit (`1.5 GiB`). For
/// reporting file sizes in build/convert summaries, where an exact byte count is
/// noise.
pub fn fmt_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{} {}", bytes, UNITS[0])
    } else {
        format!("{:.1} {}", value, UNITS[unit])
    }
}

pub fn fmt_count(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out.chars().rev().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmt_bytes_scales_and_stops_at_tib() {
        assert_eq!(fmt_bytes(0), "0 B");
        assert_eq!(fmt_bytes(512), "512 B");
        assert_eq!(fmt_bytes(1024), "1.0 KiB");
        assert_eq!(fmt_bytes(1536), "1.5 KiB");
        assert_eq!(fmt_bytes(5 * 1024 * 1024), "5.0 MiB");
        assert_eq!(fmt_bytes(3 * 1024 * 1024 * 1024), "3.0 GiB");
        assert_eq!(fmt_bytes(2 * 1024_u64.pow(4)), "2.0 TiB");
        // TiB is the largest unit, so beyond it the number keeps growing rather
        // than a unit being invented.
        assert_eq!(fmt_bytes(2048 * 1024_u64.pow(4)), "2048.0 TiB");
        assert_eq!(fmt_bytes(u64::MAX), "16777216.0 TiB");
    }

    #[test]
    fn fmt_count_commas() {
        assert_eq!(fmt_count(0), "0");
        assert_eq!(fmt_count(999), "999");
        assert_eq!(fmt_count(1_000), "1,000");
        assert_eq!(fmt_count(1_234_567), "1,234,567");
        assert_eq!(fmt_count(1_000_000_000), "1,000,000,000");
    }
}
