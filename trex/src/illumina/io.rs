//! Optional gzip-transparent input helpers.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;
use std::path::Path;

use crate::error::IngestError;
use crate::illumina::fastq::{fastq_record_id, RawFastqRecord};

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

pub fn read_fastq_records_maybe_gzip(path: &Path) -> Result<Vec<RawFastqRecord>, IngestError> {
    let mut reader = open_maybe_gzip(path).map_err(IngestError::Io)?;
    let mut records = Vec::new();
    while let Some(header) = read_trimmed_line(&mut reader).map_err(IngestError::Io)? {
        if header.is_empty() {
            continue;
        }
        let seq_line = read_trimmed_line(&mut reader)
            .map_err(IngestError::Io)?
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after header".into()))?;
        let plus = read_trimmed_line(&mut reader)
            .map_err(IngestError::Io)?
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after sequence".into()))?;
        if !plus.starts_with(b"+") {
            return Err(IngestError::FastqFormat(format!(
                "expected `+` line, got `{}`",
                String::from_utf8_lossy(&plus)
            )));
        }
        let qual = read_trimmed_line(&mut reader)
            .map_err(IngestError::Io)?
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after `+`".into()))?;
        if qual.len() != seq_line.len() {
            return Err(IngestError::FastqFormat(format!(
                "quality length {} != sequence length {}",
                qual.len(),
                seq_line.len()
            )));
        }
        records.push(RawFastqRecord {
            id: fastq_record_id(&header)?,
            sequence: seq_line,
        });
    }
    Ok(records)
}

pub fn read_fasta_records_maybe_gzip(path: &Path) -> Result<Vec<RawFastqRecord>, IngestError> {
    let mut reader = open_maybe_gzip(path).map_err(IngestError::Io)?;
    let mut records = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    while let Some(line) = read_trimmed_line(&mut reader).map_err(IngestError::Io)? {
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(id) = cur_id.take() {
                push_fasta_record(&mut records, id, std::mem::take(&mut cur_seq))?;
            }
            let mut pseudo = Vec::with_capacity(header.len().saturating_add(1));
            pseudo.push(b'@');
            pseudo.extend_from_slice(header);
            cur_id = Some(fastq_record_id(&pseudo)?);
        } else if cur_id.is_some() {
            cur_seq.extend_from_slice(&line);
        } else {
            return Err(IngestError::FastqFormat(
                "FASTA sequence line before first header".into(),
            ));
        }
    }
    if let Some(id) = cur_id {
        push_fasta_record(&mut records, id, cur_seq)?;
    }
    Ok(records)
}

fn push_fasta_record(
    records: &mut Vec<RawFastqRecord>,
    id: String,
    sequence: Vec<u8>,
) -> Result<(), IngestError> {
    if sequence.is_empty() {
        return Err(IngestError::FastqFormat(format!(
            "empty FASTA sequence for `{id}`"
        )));
    }
    records.push(RawFastqRecord { id, sequence });
    Ok(())
}

fn open_maybe_gzip(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    if file_has_gzip_magic(path)? {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

fn file_has_gzip_magic(path: &Path) -> std::io::Result<bool> {
    let mut raw = File::open(path)?;
    let mut magic = [0u8; 2];
    match raw.read_exact(&mut magic) {
        Ok(()) => Ok(magic == [0x1f, 0x8b]),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    }
}

fn read_trimmed_line(reader: &mut dyn BufRead) -> std::io::Result<Option<Vec<u8>>> {
    let mut line = Vec::new();
    let n = reader.read_until(b'\n', &mut line)?;
    if n == 0 {
        return Ok(None);
    }
    if line.ends_with(b"\n") {
        line.pop();
    }
    if line.ends_with(b"\r") {
        line.pop();
    }
    Ok(Some(line))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::{read_fasta_records_maybe_gzip, read_fastq_records_maybe_gzip, read_maybe_gzip};

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

    #[test]
    fn streams_fastq_records_from_gzip() {
        let path = std::env::temp_dir().join(format!(
            "trex-stream-fastq-{}-{}.fq.gz",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        std::fs::write(&path, gzip_member(b"@r1 x\nACGT\n+\nIIII\n")).expect("write gzip");

        let records = read_fastq_records_maybe_gzip(&path).expect("read fastq");

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "r1");
        assert_eq!(records[0].sequence, b"ACGT");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn streams_fasta_records_from_plain_file() {
        let path = std::env::temp_dir().join(format!(
            "trex-stream-fasta-{}-{}.fa",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        std::fs::write(&path, b">r1 extra\nAC\nGT\n").expect("write fasta");

        let records = read_fasta_records_maybe_gzip(&path).expect("read fasta");

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "r1");
        assert_eq!(records[0].sequence, b"ACGT");
        let _ = std::fs::remove_file(path);
    }
}
