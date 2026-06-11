//! Uppercase, **IUPAC → N** (documented table in `preprocess_sequence` doc), **Phase-1 empty read** checks.

use crate::error::IngestError;
use crate::illumina::fastq::RawFastqRecord;
use crate::illumina::read::Read;

/// Documented **Phase-1** IUPAC ambiguity handling: any symbol outside **ACGT** (ASCII) becomes **`N`**,
/// including **`N`** and IUPAC letters. Lowercase letters are uppercased first.
pub fn preprocess_sequence(id: &str, raw: &[u8]) -> Result<Vec<u8>, IngestError> {
    let mut out = Vec::with_capacity(raw.len());
    for (pos, &b) in raw.iter().enumerate() {
        let u = b.to_ascii_uppercase();
        let mapped = match u {
            b'A' | b'C' | b'G' | b'T' => u,
            b'N' => b'N',
            // IUPAC and anything else → N (Phase-1 glossary: map non-ACGT IUPAC to N).
            _ if u.is_ascii_alphabetic() => b'N',
            _ => {
                return Err(IngestError::InvalidNucleotide {
                    id: id.to_string(),
                    pos,
                    byte: b,
                });
            }
        };
        out.push(mapped);
    }
    if out.is_empty() {
        return Err(IngestError::EmptyRead(id.to_string()));
    }
    Ok(out)
}

pub fn preprocess_records(raw: Vec<RawFastqRecord>) -> Result<Vec<Read>, IngestError> {
    let mut reads = Vec::with_capacity(raw.len());
    for r in raw {
        let seq = preprocess_sequence(&r.id, &r.sequence)?;
        reads.push(Read {
            id: r.id,
            sequence: seq,
        });
    }
    Ok(reads)
}

/// Split a preprocessed read into **N**-free **ACGT** segments (**Phase-1 N policy**).
pub fn n_free_acgt_segments(seq: &[u8]) -> Vec<&[u8]> {
    let mut segments: Vec<&[u8]> = Vec::new();
    let mut start = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        if b == b'N' {
            if start < i {
                segments.push(&seq[start..i]);
            }
            start = i + 1;
        }
    }
    if start < seq.len() {
        segments.push(&seq[start..]);
    }
    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_r_maps_to_n() {
        let s = preprocess_sequence("x", b"ACGr").unwrap();
        assert_eq!(s, b"ACGN");
    }

    #[test]
    fn empty_after_strip_is_error() {
        let e = preprocess_sequence("empty", b"").unwrap_err();
        assert!(matches!(e, IngestError::EmptyRead(_)));
    }

    #[test]
    fn segments_split_on_n() {
        let s = b"ACNGT";
        let segs: Vec<_> = n_free_acgt_segments(s);
        assert_eq!(segs, vec![&b"AC"[..], &b"GT"[..]]);
    }
}
