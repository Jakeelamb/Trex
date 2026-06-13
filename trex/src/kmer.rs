//! Canonical *k*-mers: **A < C < G < T** lexicographic order, strand collapse via
//! `min(forward, reverse_complement(forward)))` per **Phase-1 k-mer identity**.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

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
        let ox = total_dna_byte_rank(x);
        let oy = total_dna_byte_rank(y);
        match ox.cmp(&oy) {
            Ordering::Equal => {}
            o => return o,
        }
    }
    a.len().cmp(&b.len())
}

#[inline]
fn total_dna_byte_rank(b: u8) -> u16 {
    dna_byte_rank(b).map(u16::from).unwrap_or(4 + u16::from(b))
}

/// Fallible reverse complement in **IUPAC-free** counted alphabet.
pub fn try_reverse_complement(seq: &[u8]) -> Result<Vec<u8>, KmerError> {
    let mut out: Vec<u8> = Vec::with_capacity(seq.len());
    for (pos, &b) in seq.iter().enumerate().rev() {
        out.push(match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => return Err(KmerError::NonAcgtByte { pos, byte: b }),
        });
    }
    Ok(out)
}

/// Reverse complement in **IUPAC-free** counted alphabet.
///
/// Non-ACGT bytes are emitted as `N`; use [`try_reverse_complement`] when invalid bases should be
/// rejected.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
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

#[derive(Debug, Clone, Default)]
pub struct CanonicalKmerSet {
    inner: HashSet<Vec<u8>>,
}

impl CanonicalKmerSet {
    pub fn new() -> Self {
        Self {
            inner: HashSet::new(),
        }
    }

    pub fn from_sequences(seqs: &[Vec<u8>], k: usize) -> Self {
        let mut set = Self::new();
        if k == 0 {
            return set;
        }
        for seq in seqs {
            for window in acgt_windows(seq, k) {
                set.insert_window(window);
            }
        }
        set
    }

    pub fn insert_window(&mut self, window: &[u8]) {
        self.inner.insert(canonical_kmer(window));
    }

    pub fn contains_window(&self, window: &[u8]) -> bool {
        self.inner.contains(&canonical_kmer(window))
    }

    pub fn contains_key(&self, key: &[u8]) -> bool {
        self.inner.contains(key)
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

#[derive(Debug, Clone, Default)]
pub struct CanonicalKmerCounts {
    inner: HashMap<Vec<u8>, u64>,
}

impl CanonicalKmerCounts {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    pub fn from_kmers(kmers: Vec<Vec<u8>>) -> Self {
        let mut counts = Self::new();
        for kmer in kmers {
            *counts.inner.entry(kmer).or_insert(0) += 1;
        }
        counts
    }

    pub fn from_sequences(seqs: &[Vec<u8>], k: usize) -> Self {
        let mut counts = Self::new();
        if k == 0 {
            return counts;
        }
        for seq in seqs {
            for window in acgt_windows(seq, k) {
                *counts.inner.entry(canonical_kmer(window)).or_insert(0) += 1;
            }
        }
        counts
    }

    pub fn get_key(&self, key: &[u8]) -> Option<u64> {
        self.inner.get(key).copied()
    }

    pub fn contains_key(&self, key: &[u8]) -> bool {
        self.inner.contains_key(key)
    }

    pub fn values(&self) -> impl Iterator<Item = u64> + '_ {
        self.inner.values().copied()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&[u8], u64)> + '_ {
        self.inner
            .iter()
            .map(|(kmer, count)| (kmer.as_slice(), *count))
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn total(&self) -> u64 {
        self.inner.values().sum()
    }

    pub fn into_sorted_rows(self) -> Vec<(Vec<u8>, u64)> {
        let mut rows: Vec<(Vec<u8>, u64)> = self.inner.into_iter().collect();
        rows.sort_unstable_by(|a, b| cmp_dna(&a.0, &b.0));
        rows
    }
}

pub fn acgt_windows(seq: &[u8], k: usize) -> impl Iterator<Item = &[u8]> {
    seq.windows(k).filter(|window| {
        window
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    })
}

pub fn packed_kmer_supported(k: usize) -> bool {
    (1..=32).contains(&k)
}

pub fn canonical_packed_kmer(window: &[u8], k: usize) -> Option<u64> {
    if window.len() != k || !packed_kmer_supported(k) {
        return None;
    }

    let mut forward = 0u64;
    for base in window {
        forward = (forward << 2) | dna_byte_rank(*base)? as u64;
    }

    let mut reverse = 0u64;
    for base in window.iter().rev() {
        reverse = (reverse << 2) | (dna_byte_rank(*base)? as u64 ^ 0b11);
    }

    Some(forward.min(reverse))
}

#[derive(Debug, Clone, Default)]
pub struct PackedCanonicalKmerSet {
    inner: HashSet<u64>,
}

impl PackedCanonicalKmerSet {
    pub fn from_sequences(seqs: &[Vec<u8>], k: usize) -> Option<Self> {
        if !packed_kmer_supported(k) {
            return None;
        }
        let mut set = Self {
            inner: HashSet::new(),
        };
        for seq in seqs {
            for window in acgt_windows(seq, k) {
                if let Some(kmer) = canonical_packed_kmer(window, k) {
                    set.inner.insert(kmer);
                }
            }
        }
        Some(set)
    }

    pub fn contains_window(&self, window: &[u8], k: usize) -> bool {
        canonical_packed_kmer(window, k)
            .map(|kmer| self.inner.contains(&kmer))
            .unwrap_or(false)
    }
}

#[derive(Debug, Clone, Default)]
pub struct PackedCanonicalKmerCounts {
    inner: HashMap<u64, u64>,
}

impl PackedCanonicalKmerCounts {
    pub fn from_sequences(seqs: &[Vec<u8>], k: usize) -> Option<Self> {
        if !packed_kmer_supported(k) {
            return None;
        }
        let mut counts = Self {
            inner: HashMap::new(),
        };
        for seq in seqs {
            for window in acgt_windows(seq, k) {
                if let Some(kmer) = canonical_packed_kmer(window, k) {
                    *counts.inner.entry(kmer).or_insert(0) += 1;
                }
            }
        }
        Some(counts)
    }

    pub fn values(&self) -> impl Iterator<Item = u64> + '_ {
        self.inner.values().copied()
    }

    pub fn iter(&self) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.inner.iter().map(|(kmer, count)| (*kmer, *count))
    }

    pub fn contains_key(&self, key: u64) -> bool {
        self.inner.contains_key(&key)
    }

    pub fn get_key(&self, key: u64) -> Option<u64> {
        self.inner.get(&key).copied()
    }

    pub fn keys(&self) -> impl Iterator<Item = u64> + '_ {
        self.inner.keys().copied()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn total(&self) -> u64 {
        self.inner.values().sum()
    }
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

    #[test]
    fn canonical_counts_emit_dna_sorted_rows() {
        let kmers = vec![b"TTT".to_vec(), b"AAA".to_vec(), b"AAA".to_vec()];
        let rows = CanonicalKmerCounts::from_kmers(kmers).into_sorted_rows();
        assert_eq!(rows, vec![(b"AAA".to_vec(), 2), (b"TTT".to_vec(), 1)]);
    }

    #[test]
    fn canonical_set_contains_reverse_complement_windows() {
        let set = CanonicalKmerSet::from_sequences(&[b"ACCGTTAA".to_vec()], 4);
        assert!(set.contains_window(b"AACG"));
        assert_eq!(set.len(), 5);
    }

    #[test]
    fn public_reverse_complement_does_not_panic_on_non_acgt() {
        assert_eq!(reverse_complement(b"ACNX"), b"NNGT");
        let err = try_reverse_complement(b"ACNX").unwrap_err();
        assert!(matches!(err, KmerError::NonAcgtByte { pos: 3, byte: b'X' }));
    }

    #[test]
    fn cmp_dna_has_total_order_for_non_acgt_bytes() {
        assert_eq!(cmp_dna(b"ACG", b"ACN"), Ordering::Less);
        assert_eq!(cmp_dna(b"ACN", b"ACX"), Ordering::Less);
    }
}
