//! Optional gzip transparent read.

use std::io::Read;

use flate2::read::GzDecoder;
use std::path::Path;

/// Read full file; if gzip magic is present, decompress to bytes (Phase-1 small-input path).
pub fn read_maybe_gzip(path: &Path) -> std::io::Result<Vec<u8>> {
    let raw = std::fs::read(path)?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut decoder = GzDecoder::new(&raw[..]);
        let mut out = Vec::new();
        decoder.read_to_end(&mut out)?;
        Ok(out)
    } else {
        Ok(raw)
    }
}
