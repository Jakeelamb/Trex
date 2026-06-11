//! **Phase-2 Illumina mate usage**: strengthen existing **DBG** edges using paired-end layout
//! (does not alter **Phase-1 k-mer count representation** or trusted vertex weights).

use crate::dbg::graph::DbgGraph;
use crate::illumina::preprocess::n_free_acgt_segments;
use crate::illumina::read::Read;
use crate::kmer::canonical_kmer;

fn first_canonical_kmer_forward(seq: &[u8], k: usize) -> Option<Vec<u8>> {
    for seg in n_free_acgt_segments(seq) {
        if seg.len() >= k {
            return Some(canonical_kmer(&seg[..k]));
        }
    }
    None
}

fn last_canonical_kmer_forward(seq: &[u8], k: usize) -> Option<Vec<u8>> {
    let mut last: Option<Vec<u8>> = None;
    for seg in n_free_acgt_segments(seq) {
        if seg.len() >= k {
            let w = &seg[seg.len() - k..];
            last = Some(canonical_kmer(w));
        }
    }
    last
}

/// For each mate pair, if **R1**'s last forward *k*-mer and **R2**'s first forward *k*-mer are both
/// trusted vertices and already adjacent in the **DBG**, increment that undirected edge weight by **1**.
///
/// This is a conservative **Phase-2 Illumina** bridge: no new edges, no changes to `node_mul`.
/// It runs only when the operator supplied an **insert-size prior** (see `assemble_illumina` caller).
pub fn boost_mate_pairs_on_existing_dbg_edges(
    graph: &mut DbgGraph,
    reads: &[Read],
    r1_count: usize,
    k: usize,
) -> usize {
    if r1_count == 0 || reads.len() < r1_count * 2 {
        return 0;
    }
    let mut boosted = 0usize;
    for i in 0..r1_count {
        let r1 = &reads[i];
        let r2 = &reads[i + r1_count];
        let Some(a) = last_canonical_kmer_forward(&r1.sequence, k) else {
            continue;
        };
        let Some(b) = first_canonical_kmer_forward(&r2.sequence, k) else {
            continue;
        };
        if a == b {
            continue;
        }
        if !graph.node_mul.contains_key(&a) || !graph.node_mul.contains_key(&b) {
            continue;
        }
        let w0 = graph
            .adj
            .get(&a)
            .and_then(|m| m.get(&b))
            .copied()
            .unwrap_or(0);
        if w0 == 0 {
            continue;
        }
        if graph.add_undirected_edge(&a, &b, 1).is_ok() {
            boosted += 1;
        }
    }
    boosted
}
