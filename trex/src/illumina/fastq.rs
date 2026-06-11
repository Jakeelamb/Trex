//! FASTQ parsing (four-line records). Quality lines are read and discarded for Phase-1 graph input.

use crate::error::IngestError;

/// Raw sequence line (still needs **Phase-1** preprocess: case + IUPAC).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawFastqRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

/// Extract FASTQ record id: first whitespace-separated token on the header line (without leading `@`).
pub fn fastq_record_id(header_line: &[u8]) -> Result<String, IngestError> {
    let s = header_line
        .strip_prefix(b"@")
        .ok_or_else(|| IngestError::FastqFormat("header must start with `@`".into()))?;
    let id_end = s
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(s.len());
    let id = s[..id_end].iter().map(|&b| b as char).collect::<String>();
    if id.is_empty() {
        return Err(IngestError::FastqFormat("empty FASTQ id".into()));
    }
    Ok(id)
}

/// Parse FASTQ from bytes (LF or CRLF). Sequence line is raw ASCII; must be non-empty per record.
pub fn parse_fastq(content: &[u8]) -> Result<Vec<RawFastqRecord>, IngestError> {
    let mut records = Vec::new();
    let mut lines = content.split(|&b| b == b'\n');
    loop {
        let header = match lines.next() {
            Some(h) => h,
            None => break,
        };
        let header = trim_cr(header);
        if header.is_empty() {
            continue;
        }
        let seq_line = lines
            .next()
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after header".into()))?;
        let seq_line = trim_cr(seq_line);
        let plus = lines
            .next()
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after sequence".into()))?;
        let plus = trim_cr(plus);
        if !plus.starts_with(b"+") {
            return Err(IngestError::FastqFormat(format!(
                "expected `+` line, got `{}`",
                String::from_utf8_lossy(plus)
            )));
        }
        let qual = lines
            .next()
            .ok_or_else(|| IngestError::FastqFormat("truncated FASTQ after `+`".into()))?;
        let qual = trim_cr(qual);
        if qual.len() != seq_line.len() {
            return Err(IngestError::FastqFormat(format!(
                "quality length {} != sequence length {}",
                qual.len(),
                seq_line.len()
            )));
        }
        let id = fastq_record_id(header)?;
        records.push(RawFastqRecord {
            id,
            sequence: seq_line.to_vec(),
        });
    }
    Ok(records)
}

fn trim_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_one_record() {
        let fq = b"@r1/1\nACGT\n+\nIIII\n";
        let r = parse_fastq(fq).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].id, "r1/1");
        assert_eq!(r[0].sequence, b"ACGT");
    }
}
