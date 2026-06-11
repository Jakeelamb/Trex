//! Canonical *k*-mers: **A < C < G < T** lexicographic order, strand collapse via
//! `min(forward, reverse_complement(forward)))` per **Phase-1 k-mer identity**.

use std::cmp::Ordering;

use crate::error::KmerError;

/// Lexicographic rank for **A/C/G/T** only (`N` must not appear in counted *k*-mers).
#[inline]
pub fn dna_byte_rank(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Compare two DNA strings using **Phase-1 canonical alphabet** order.
pub fn cmp_dna(a: &[u8], b: &[u8]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()) {
        let ox = dna_byte_rank(x).expect("cmp_dna: non-ACGT");
        let oy = dna_byte_rank(y).expect("cmp_dna: non-ACGT");
        match ox.cmp(&oy) {
            Ordering::Equal => {}
            o => return o,
        }
    }
    a.len().cmp(&b.len())
}

/// Reverse complement in **IUPAC-free** counted alphabet (input must be ACGT upper).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("reverse_complement: non-ACGT byte {b}"),
        });
    }
    out
}

/// Canonical *k*-mer key for a window of length `k` (**ACGT** only).
pub fn canonical_kmer(window: &[u8]) -> Vec<u8> {
    debug_assert!(!window.is_empty());
    let rc = reverse_complement(window);
    if cmp_dna(window, &rc) == Ordering::Greater {
        rc
    } else {
        window.to_vec()
    }
}

/// Enumerate canonical *k*-mers from one **N**-free segment (ASCII **ACGT**).
pub fn kmers_from_segment(segment: &[u8], k: usize) -> Result<Vec<Vec<u8>>, KmerError> {
    if k == 0 {
        return Err(KmerError::KZero(k));
    }
    if segment.len() < k {
        return Err(KmerError::KLongerThanSegment {
            k,
            segment_len: segment.len(),
        });
    }
    let mut out = Vec::with_capacity(segment.len().saturating_sub(k - 1));
    for w in segment.windows(k) {
        out.push(canonical_kmer(w));
    }
    Ok(out)
}

/// Sort (non-stable) then merge equal adjacent keys into counts.
pub fn sort_and_merge_counts(mut kmers: Vec<Vec<u8>>) -> Vec<(Vec<u8>, u64)> {
    #[cfg(feature = "parallel-kmer-sort")]
    {
        if kmers.len() > 4096 {
            use rayon::slice::ParallelSliceMut;
            kmers.par_sort_unstable_by(|a, b| cmp_dna(a, b));
        } else {
            kmers.sort_unstable_by(|a, b| cmp_dna(a, b));
        }
    }
    #[cfg(not(feature = "parallel-kmer-sort"))]
    {
        kmers.sort_unstable_by(|a, b| cmp_dna(a, b));
    }
    let mut merged: Vec<(Vec<u8>, u64)> = Vec::new();
    for kmer in kmers {
        match merged.last_mut() {
            Some((last, c)) if cmp_dna(last, &kmer) == Ordering::Equal => {
                *c += 1;
            }
            _ => merged.push((kmer, 1)),
        }
    }
    merged
}

/// Keep *k*-mers whose total multiplicity is **>= trusted_threshold** *T*.
pub fn apply_trusted_threshold(
    rows: Vec<(Vec<u8>, u64)>,
    trusted_threshold: u64,
) -> Vec<(Vec<u8>, u64)> {
    rows.into_iter()
        .filter(|(_, c)| *c >= trusted_threshold)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_prefers_min_under_acgt_order() {
        // RC("AG") = "CT"; lex min is "AG".
        let w = b"AG";
        assert_eq!(canonical_kmer(w), b"AG");
        // "TA" is its own reverse complement (2bp), so canonical equals forward.
        assert_eq!(canonical_kmer(b"TA"), b"TA");
    }

    #[test]
    fn merge_counts_after_unsorted_input() {
        let v = vec![b"AAA".to_vec(), b"TTT".to_vec(), b"AAA".to_vec()];
        let m = sort_and_merge_counts(v);
        assert_eq!(m, vec![(b"AAA".to_vec(), 2), (b"TTT".to_vec(), 1)]);
    }

    #[test]
    fn trusted_threshold_filters() {
        let rows = vec![(b"AAA".to_vec(), 1), (b"CCC".to_vec(), 3)];
        let f = apply_trusted_threshold(rows, 2);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].0, b"CCC");
    }
}
