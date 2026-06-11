//! *k*-mer enumeration from preprocessed reads and **Phase-1 k-mer count representation**.

use crate::error::KmerError;
use crate::illumina::preprocess::n_free_acgt_segments;
use crate::illumina::read::Read;
use crate::kmer::{self, kmers_from_segment};

/// Enumerate multiset of canonical *k*-mers, then sort-merge into counts.
pub fn enumerate_sorted_counts(reads: &[Read], k: usize) -> Result<Vec<(Vec<u8>, u64)>, KmerError> {
    let mut all: Vec<Vec<u8>> = Vec::new();
    for read in reads {
        for seg in n_free_acgt_segments(&read.sequence) {
            if seg.len() >= k {
                all.extend(kmers_from_segment(seg, k)?);
            }
        }
    }
    Ok(kmer::sort_and_merge_counts(all))
}
