//! **Phase-2 Illumina primary FASTA** — deterministic **A/C/G/T** collapse using trusted *k*-mer
//! multiplicity over overlapping windows (no **IUPAC** / **`N`** in primary stream).

use std::collections::HashMap;

use crate::kmer::{canonical_kmer, dna_byte_rank};

/// Per-base vote: sum trusted multiplicity over every length-*k* window that contains the base,
/// trying **A/C/G/T** at that base. Tie-break at equal score: **A < C < G < T** (same as Phase-1
/// canonical lex order). Positions with **zero** total score for all trials keep the stitched base.
pub fn collapse_primary_contig_by_trusted_kmers(
    sequence: &mut [u8],
    k: usize,
    trusted: &HashMap<Vec<u8>, u64>,
) {
    let len = sequence.len();
    if len == 0 || k == 0 || len < k {
        return;
    }
    let max_s = len - k;
    let mut buf = vec![0u8; k];
    for i in 0..len {
        let s_lo = i.saturating_sub(k - 1);
        let s_hi = max_s.min(i);
        if s_lo > s_hi {
            continue;
        }
        let mut best_score = 0u64;
        let mut best_rank = 4u8;
        let mut best_b = sequence[i];
        for &trial in b"ACGT" {
            let mut total = 0u64;
            for s in s_lo..=s_hi {
                buf.copy_from_slice(&sequence[s..s + k]);
                let oi = i - s;
                buf[oi] = trial;
                let key = canonical_kmer(&buf);
                total += *trusted.get(&key).unwrap_or(&0);
            }
            let trial_rank = dna_byte_rank(trial).expect("ACGT");
            if total > best_score || (total == best_score && trial_rank < best_rank) {
                best_score = total;
                best_rank = trial_rank;
                best_b = trial;
            }
        }
        if best_score > 0 {
            sequence[i] = best_b;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn collapse_prefers_trusted_base() {
        let mut seq = b"AAAA".to_vec();
        let mut m = HashMap::new();
        m.insert(crate::kmer::canonical_kmer(b"AAAT"), 10u64);
        m.insert(crate::kmer::canonical_kmer(b"AAAC"), 2u64);
        collapse_primary_contig_by_trusted_kmers(&mut seq, 4, &m);
        assert_eq!(seq[3], b'T');
    }
}
