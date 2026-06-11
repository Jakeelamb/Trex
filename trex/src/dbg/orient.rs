//! Pick one **forward** *k*-mer representative per canonical node, consistent with observed reads.

use std::collections::HashMap;

use crate::error::GraphError;
use crate::illumina::preprocess::n_free_acgt_segments;
use crate::illumina::read::Read;
use crate::kmer::{canonical_kmer, reverse_complement};

/// First-seen forward window wins; later windows must match it or its reverse complement.
pub fn forward_representatives(
    reads: &[Read],
    k: usize,
) -> Result<HashMap<Vec<u8>, Vec<u8>>, GraphError> {
    let mut m: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for read in reads {
        for seg in n_free_acgt_segments(&read.sequence) {
            if seg.len() < k {
                continue;
            }
            for i in 0..=(seg.len() - k) {
                let f = seg[i..i + k].to_vec();
                let c = canonical_kmer(&f);
                match m.get(&c) {
                    None => {
                        m.insert(c, f);
                    }
                    Some(existing) => {
                        let rc_f = reverse_complement(&f);
                        if existing != &f && existing != &rc_f {
                            return Err(GraphError::OrientationConflict);
                        }
                    }
                }
            }
        }
    }
    Ok(m)
}
