//! Maximal **non-branching** chains in the simplified DBG, as ordered lists of canonical *k*-mers.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::dbg::graph::DbgGraph;
use crate::error::GraphError;
use crate::kmer::{cmp_dna, reverse_complement};

fn norm_edge(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    if cmp_dna(a, b) == Ordering::Greater {
        (b.to_vec(), a.to_vec())
    } else {
        (a.to_vec(), b.to_vec())
    }
}

fn mark_used(a: &[u8], b: &[u8], used: &mut BTreeSet<(Vec<u8>, Vec<u8>)>) {
    used.insert(norm_edge(a, b));
}

fn min_neighbor(graph: &DbgGraph, u: &[u8]) -> Vec<u8> {
    graph.adj[u]
        .keys()
        .min_by(|a, b| cmp_dna(a, b))
        .expect("neighbor")
        .clone()
}

/// Extract vertex chains (canonical *k*-mers along each unitig).
pub fn extract_unitigs(graph: &DbgGraph) -> Vec<Vec<Vec<u8>>> {
    let mut used: BTreeSet<(Vec<u8>, Vec<u8>)> = BTreeSet::new();
    let mut out: Vec<Vec<Vec<u8>>> = Vec::new();

    let breakpoints: Vec<Vec<u8>> = graph
        .adj
        .keys()
        .filter(|u| graph.degree(u) != 2)
        .cloned()
        .collect();

    if breakpoints.is_empty() && !graph.adj.is_empty() {
        extract_one_cycle(graph, &mut used, &mut out);
        return out;
    }

    for u in &breakpoints {
        for v in graph.adj.get(u).into_iter().flat_map(|m| m.keys()) {
            let e = norm_edge(u, v);
            if used.contains(&e) {
                continue;
            }
            let path = walk_from_break(graph, u, v, &mut used);
            if path.len() >= 2 {
                out.push(path);
            }
        }
    }
    out
}

fn walk_from_break(
    graph: &DbgGraph,
    u: &[u8],
    v: &[u8],
    used: &mut BTreeSet<(Vec<u8>, Vec<u8>)>,
) -> Vec<Vec<u8>> {
    let mut path = vec![u.to_vec(), v.to_vec()];
    mark_used(u, v, used);
    let mut prev = u.to_vec();
    let mut cur = v.to_vec();
    loop {
        if graph.degree(&cur) != 2 {
            break;
        }
        let nbrs: Vec<Vec<u8>> = graph.adj[&cur]
            .keys()
            .filter(|n| *n != &prev)
            .cloned()
            .collect();
        if nbrs.len() != 1 {
            break;
        }
        let nxt = nbrs.into_iter().next().unwrap();
        let e = norm_edge(&cur, &nxt);
        if used.contains(&e) {
            break;
        }
        mark_used(&cur, &nxt, used);
        path.push(nxt.clone());
        prev = cur;
        cur = nxt;
    }
    path
}

/// Single **2-regular** component (one or more disjoint cycles): emit one cyclic unitig from lex-min start.
fn extract_one_cycle(
    graph: &DbgGraph,
    used: &mut BTreeSet<(Vec<u8>, Vec<u8>)>,
    out: &mut Vec<Vec<Vec<u8>>>,
) {
    let start = graph
        .adj
        .keys()
        .min_by(|a, b| cmp_dna(a, b))
        .unwrap()
        .clone();
    let v0 = min_neighbor(graph, &start);
    let e0 = norm_edge(&start, &v0);
    if used.contains(&e0) {
        return;
    }
    let cyc = walk_cycle(graph, &start, &v0, used);
    if cyc.len() >= 3 {
        out.push(cyc);
    }
}

fn walk_cycle(
    graph: &DbgGraph,
    start: &[u8],
    first: &[u8],
    used: &mut BTreeSet<(Vec<u8>, Vec<u8>)>,
) -> Vec<Vec<u8>> {
    let mut path = vec![start.to_vec(), first.to_vec()];
    mark_used(start, first, used);
    let mut prev = start.to_vec();
    let mut cur = first.to_vec();
    let limit = graph.adj.len().saturating_mul(4).max(8);
    for _ in 0..limit {
        let nbrs: Vec<Vec<u8>> = graph.adj[&cur]
            .keys()
            .filter(|n| *n != &prev)
            .cloned()
            .collect();
        if nbrs.len() != 1 {
            break;
        }
        let nxt = nbrs.into_iter().next().unwrap();
        let e = norm_edge(&cur, &nxt);
        if used.contains(&e) {
            break;
        }
        mark_used(&cur, &nxt, used);
        path.push(nxt.clone());
        if nxt == *start {
            break;
        }
        prev = cur;
        cur = nxt;
    }
    path
}

/// Stitch canonical *k*-mers into one sequence using per-node forward representatives.
///
/// When a path stitches **read-consistent** edges that may join different reads, a single global
/// forward strand choice may not exist; this implementation tries both orientations of the first
/// *k*-mer and, at each step, both orientations of the next, picking the **lexicographically smallest**
/// complete sequence among valid completions (**Phase-1 contig tie-break** extended for multi-read paths).
pub fn stitch_sequence(
    path: &[Vec<u8>],
    forward: &std::collections::HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
) -> Result<Vec<u8>, GraphError> {
    if path.is_empty() {
        return Ok(Vec::new());
    }
    if path.len() == 1 {
        return pick_forward(path[0].as_slice(), forward);
    }
    let f0 = pick_forward(path[0].as_slice(), forward)?;
    let rc0 = reverse_complement(&f0);
    let mut starters: Vec<Vec<u8>> = if f0 == rc0 { vec![f0] } else { vec![f0, rc0] };
    starters.sort_by(|a, b| cmp_dna(a, b));
    starters.dedup();

    let mut best: Option<Vec<u8>> = None;
    for s0 in starters {
        if let Ok(seq) = stitch_from(path, forward, k, s0) {
            let replace = match &best {
                None => true,
                Some(b) => cmp_dna(&seq, b) == Ordering::Less,
            };
            if replace {
                best = Some(seq);
            }
        }
    }
    best.ok_or(GraphError::OrientationConflict)
}

fn stitch_from(
    path: &[Vec<u8>],
    forward: &std::collections::HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    mut cur: Vec<u8>,
) -> Result<Vec<u8>, GraphError> {
    let mut seq = cur.clone();
    for cn in &path[1..] {
        let cand = pick_forward(cn.as_slice(), forward)?;
        let rc_c = reverse_complement(&cand);
        let mut opts: Vec<Vec<u8>> = Vec::new();
        if cur[1..] == cand[..k - 1] {
            opts.push(cand.clone());
        }
        if cand != rc_c && cur[1..] == rc_c[..k - 1] {
            opts.push(rc_c);
        }
        if opts.is_empty() {
            return Err(GraphError::OrientationConflict);
        }
        opts.sort_by(|a, b| cmp_dna(a, b));
        let use_c = opts.into_iter().next().expect("non-empty");
        seq.push(*use_c.last().expect("k-mer"));
        cur = use_c;
    }
    Ok(seq)
}

fn pick_forward(
    canonical: &[u8],
    forward: &std::collections::HashMap<Vec<u8>, Vec<u8>>,
) -> Result<Vec<u8>, GraphError> {
    forward
        .get(canonical)
        .cloned()
        .ok_or(GraphError::OrientationConflict)
}

#[cfg(test)]
mod tests {
    use super::stitch_sequence;
    use crate::kmer::canonical_kmer;
    use std::collections::HashMap;

    #[test]
    fn stitch_sequence_dual_start_when_canonical_prefix_fails() {
        // Canonical keys ATC / AAT with read orientations GAT and ATT so overlap exists only when
        // the stitched walk starts on GAT (non-canonical read orientation of the first vertex).
        let k = 3;
        let key0 = canonical_kmer(b"GAT");
        let key1 = canonical_kmer(b"ATT");
        assert_eq!(key0, b"ATC".as_slice());
        assert_eq!(key1, b"AAT".as_slice());

        let mut forward = HashMap::new();
        forward.insert(key0.to_vec(), b"GAT".to_vec());
        forward.insert(key1.to_vec(), b"ATT".to_vec());

        let path = vec![key0.to_vec(), key1.to_vec()];
        let seq = stitch_sequence(&path, &forward, k).expect("stitch");
        assert_eq!(seq, b"GATT");
    }

    #[test]
    fn stitch_sequence_next_step_may_use_reverse_complement_of_candidate() {
        let k = 3;
        let key0 = canonical_kmer(b"AAA");
        let key1_canon = canonical_kmer(b"GTT");
        assert_eq!(key1_canon, b"AAC".as_slice());

        let mut forward = HashMap::new();
        forward.insert(key0.to_vec(), b"AAA".to_vec());
        forward.insert(key1_canon.to_vec(), b"GTT".to_vec());

        let path = vec![key0.to_vec(), key1_canon.to_vec()];
        let seq = stitch_sequence(&path, &forward, k).expect("stitch");
        assert_eq!(seq, b"AAAC");
    }

    #[test]
    fn stitch_sequence_single_vertex_is_pick_forward() {
        let k = 4;
        let key = canonical_kmer(b"ACGT");
        let mut forward = HashMap::new();
        forward.insert(key.clone(), b"ACGT".to_vec());
        let seq = stitch_sequence(&[key], &forward, k).expect("stitch");
        assert_eq!(seq, b"ACGT");
    }
}
