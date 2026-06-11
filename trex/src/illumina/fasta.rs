//! **FASTA** ingest (no qualities). Headers reuse the same **id token** rules as FASTQ (**Phase-1
//! paired read layout**): first whitespace-separated token, with a synthetic `@` prefix for the
//! shared parser.

use crate::error::IngestError;
use crate::illumina::fastq::{fastq_record_id, RawFastqRecord};

fn record_id_from_header(header: &[u8]) -> Result<String, IngestError> {
    let mut pseudo = Vec::with_capacity(header.len().saturating_add(1));
    pseudo.push(b'@');
    pseudo.extend_from_slice(header);
    fastq_record_id(&pseudo)
}

/// Parse FASTA from bytes (`>` headers, optional line breaks in sequence). Supports LF/CRLF.
pub fn parse_fasta(content: &[u8]) -> Result<Vec<RawFastqRecord>, IngestError> {
    let mut records = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line in content.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b">") {
            if let Some(id) = cur_id.take() {
                if cur_seq.is_empty() {
                    return Err(IngestError::FastqFormat(format!(
                        "empty FASTA sequence for `{id}`"
                    )));
                }
                records.push(RawFastqRecord {
                    id,
                    sequence: std::mem::take(&mut cur_seq),
                });
            }
            let header = &line[1..];
            cur_id = Some(record_id_from_header(header)?);
        } else if cur_id.is_some() {
            cur_seq.extend_from_slice(line);
        } else {
            return Err(IngestError::FastqFormat(
                "FASTA sequence line before first header".into(),
            ));
        }
    }
    if let Some(id) = cur_id.take() {
        if cur_seq.is_empty() {
            return Err(IngestError::FastqFormat(format!(
                "empty FASTA sequence for `{id}`"
            )));
        }
        records.push(RawFastqRecord {
            id,
            sequence: cur_seq,
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_two_records() {
        let b = b">r1 extra\nACGT\n>r2\nTT\n";
        let r = parse_fasta(b).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].id, "r1");
        assert_eq!(r[0].sequence, b"ACGT");
        assert_eq!(r[1].sequence, b"TT");
    }
}
