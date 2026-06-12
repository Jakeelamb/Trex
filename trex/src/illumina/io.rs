//! Optional gzip transparent read.

use std::io::Read;

use flate2::read::MultiGzDecoder;
use std::path::Path;

/// Read full file; if gzip magic is present, decompress to bytes (Phase-1 small-input path).
pub fn read_maybe_gzip(path: &Path) -> std::io::Result<Vec<u8>> {
    let raw = std::fs::read(path)?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut decoder = MultiGzDecoder::new(&raw[..]);
        let mut out = Vec::new();
        decoder.read_to_end(&mut out)?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::read_maybe_gzip;

    fn gzip_member(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(bytes).expect("write gzip member");
        encoder.finish().expect("finish gzip member")
    }

    #[test]
    fn reads_concatenated_gzip_members() {
        let mut payload = gzip_member(b"@r1\nAC");
        payload.extend(gzip_member(b"GT\n+\nIIII\n"));
        let path = std::env::temp_dir().join(format!(
            "trex-multigzip-{}-{}.fq.gz",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        std::fs::write(&path, payload).expect("write gzip fixture");

        let decoded = read_maybe_gzip(&path).expect("read gzip");

        assert_eq!(decoded, b"@r1\nACGT\n+\nIIII\n");
        let _ = std::fs::remove_file(path);
    }
}
