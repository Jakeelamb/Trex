//! **Phase-1 reference contig walker**: one **vertex-simple** path per connected component,
//! ranked by total **read-derived edge multiplicity** along arcs; ties break on **lexicographically
//! smallest** stitched sequence (**Phase-1 contig walk score** + **Phase-1 contig tie-break**).
//!
//! The search is **deterministic greedy extension** from every vertex seed in the component,
//! keeping the best scoring path found.

use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::dbg::graph::{CompactDbgGraph, DbgGraph, NodeId};
#[cfg(test)]
use crate::dbg::unitig::stitch_sequence;
use crate::error::GraphError;
use crate::kmer::{cmp_dna, reverse_complement};

/// Tie-breaking when two forward neighbors share the same **read-edge multiplicity** from `cur`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ContigWalkTieBreak {
    /// **Phase-1 contig tie-break**: lexicographically smallest neighbor *k*-mer.
    #[default]
    Phase1Lex,
    /// **Phase-2 Illumina (experimental)**: higher **trusted vertex multiplicity** on the neighbor, then lex smallest.
    Phase2DiploidNodeMul,
}

fn id_path_to_names(g: &CompactDbgGraph, path: &[NodeId]) -> Vec<Vec<u8>> {
    path.iter().map(|id| g.node_name_vec(*id)).collect()
}

fn walk_edge_score(g: &CompactDbgGraph, path: &[NodeId]) -> u64 {
    if path.len() < 2 {
        return 0;
    }
    let mut s = 0u64;
    for w in path.windows(2) {
        s += g.edge_weight(w[0], w[1]);
    }
    s
}

fn pick_best_neighbor(
    g: &CompactDbgGraph,
    cur: NodeId,
    forbidden: &HashSet<NodeId>,
    tie_break: ContigWalkTieBreak,
) -> Option<(NodeId, u64)> {
    let mut best: Option<(u64, NodeId)> = None;
    for (nb, w) in g.neighbors(cur) {
        if forbidden.contains(nb) {
            continue;
        }
        match &best {
            None => best = Some((*w, *nb)),
            Some((bw, bnb)) => {
                let take = match w.cmp(bw) {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal => match tie_break {
                        ContigWalkTieBreak::Phase1Lex => {
                            cmp_dna(g.node_name(*nb), g.node_name(*bnb)) == Ordering::Less
                        }
                        ContigWalkTieBreak::Phase2DiploidNodeMul => {
                            let mn = g.node_mul(*nb);
                            let mb = g.node_mul(*bnb);
                            mn > mb
                                || (mn == mb
                                    && cmp_dna(g.node_name(*nb), g.node_name(*bnb))
                                        == Ordering::Less)
                        }
                    },
                };
                if take {
                    best = Some((*w, *nb));
                }
            }
        }
    }
    best.map(|(w, nb)| (nb, w))
}

/// Greedy **vertex-simple** path: extend forward from `seed`, then backward from `seed` without
/// reusing vertices from the forward segment (except `seed`).
fn greedy_simple_path(
    g: &CompactDbgGraph,
    seed: NodeId,
    tie_break: ContigWalkTieBreak,
) -> (u64, Vec<NodeId>) {
    let mut forward_forbidden: HashSet<NodeId> = HashSet::new();
    forward_forbidden.insert(seed);

    let mut score = 0u64;
    let mut forward: Vec<NodeId> = vec![seed];
    let mut cur = seed;
    while let Some((nb, edge_weight)) = pick_best_neighbor(g, cur, &forward_forbidden, tie_break) {
        score += edge_weight;
        forward_forbidden.insert(nb);
        forward.push(nb);
        cur = nb;
    }

    let mut all_forbidden = forward_forbidden.clone();
    let mut back_segments: Vec<NodeId> = Vec::new();
    cur = seed;
    while let Some((nb, edge_weight)) = pick_best_neighbor(g, cur, &all_forbidden, tie_break) {
        score += edge_weight;
        all_forbidden.insert(nb);
        back_segments.push(nb);
        cur = nb;
    }

    let mut path: Vec<NodeId> = back_segments.into_iter().rev().collect();
    path.extend(forward);
    (score, path)
}

pub fn connected_components(g: &DbgGraph) -> Vec<Vec<Vec<u8>>> {
    let compact = g.compact_view();
    connected_component_ids(&compact)
        .into_iter()
        .map(|component| id_path_to_names(&compact, &component))
        .collect()
}

fn connected_component_ids(g: &CompactDbgGraph) -> Vec<Vec<NodeId>> {
    let mut unseen: BTreeSet<NodeId> = g.node_ids().filter(|id| g.degree(*id) > 0).collect();
    let mut comps: Vec<Vec<NodeId>> = Vec::new();
    while let Some(seed) = unseen.iter().next().cloned() {
        unseen.remove(&seed);
        let mut stack = vec![seed];
        let mut comp: Vec<NodeId> = Vec::new();
        while let Some(u) = stack.pop() {
            comp.push(u);
            for (v, _) in g.neighbors(u) {
                if unseen.remove(v) {
                    stack.push(*v);
                }
            }
        }
        comp.sort_unstable_by(|a, b| cmp_dna(g.node_name(*a), g.node_name(*b)));
        comps.push(comp);
    }
    comps.sort_unstable_by(|a, b| cmp_dna(g.node_name(a[0]), g.node_name(b[0])));
    comps
}

fn linear_component_path_ids(g: &CompactDbgGraph, comp: &[NodeId]) -> Option<Vec<NodeId>> {
    if comp.is_empty() || comp.iter().any(|u| g.degree(*u) > 2) {
        return None;
    }

    let start = comp
        .iter()
        .filter(|u| g.degree(**u) <= 1)
        .min_by(|a, b| cmp_dna(g.node_name(**a), g.node_name(**b)))
        .or_else(|| {
            comp.iter()
                .min_by(|a, b| cmp_dna(g.node_name(**a), g.node_name(**b)))
        })
        .copied()?;

    let mut path = Vec::with_capacity(comp.len());
    let mut seen: HashSet<NodeId> = HashSet::with_capacity(comp.len());
    let mut prev: Option<NodeId> = None;
    let mut cur = start;

    loop {
        if !seen.insert(cur) {
            break;
        }
        path.push(cur);

        let next = g
            .neighbors(cur)
            .iter()
            .map(|(nb, _)| *nb)
            .filter(|nb| prev.map(|p| *nb != p).unwrap_or(true))
            .find(|nb| !seen.contains(nb));
        let Some(nxt) = next else {
            break;
        };
        prev = Some(cur);
        cur = nxt;
    }

    if path.len() == comp.len() {
        Some(path)
    } else {
        None
    }
}

#[cfg(test)]
fn pick_best_stitchable_path<I>(
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    candidates: I,
) -> Result<Option<Vec<Vec<u8>>>, GraphError>
where
    I: IntoIterator<Item = (u64, Vec<Vec<u8>>)>,
{
    let mut best: Option<(u64, Vec<u8>, Vec<Vec<u8>>)> = None;
    for (score, path) in candidates {
        if best
            .as_ref()
            .map(|(best_score, _, _)| score < *best_score)
            .unwrap_or(false)
        {
            continue;
        }
        let seq = match stitch_sequence(&path, forward, k) {
            Ok(seq) => seq,
            Err(GraphError::OrientationConflict) => continue,
            Err(e) => return Err(e),
        };
        let replace = match &best {
            None => true,
            Some((bs, bseq, _)) => {
                score > *bs || (score == *bs && cmp_dna(&seq, bseq) == Ordering::Less)
            }
        };
        if replace {
            best = Some((score, seq, path));
        }
    }
    Ok(best.map(|(_, _, path)| path))
}

fn pick_best_stitchable_id_path<I>(
    g: &CompactDbgGraph,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    candidates: I,
) -> Result<Option<Vec<NodeId>>, GraphError>
where
    I: IntoIterator<Item = (u64, Vec<NodeId>)>,
{
    let mut best: Option<(u64, Vec<u8>, Vec<NodeId>)> = None;
    for (score, path) in candidates {
        if best
            .as_ref()
            .map(|(best_score, _, _)| score < *best_score)
            .unwrap_or(false)
        {
            continue;
        }
        let seq = match stitch_id_path(g, &path, forward, k) {
            Ok(seq) => seq,
            Err(GraphError::OrientationConflict) => continue,
            Err(e) => return Err(e),
        };
        let replace = match &best {
            None => true,
            Some((bs, bseq, _)) => {
                score > *bs || (score == *bs && cmp_dna(&seq, bseq) == Ordering::Less)
            }
        };
        if replace {
            best = Some((score, seq, path));
        }
    }
    Ok(best.map(|(_, _, path)| path))
}

fn stitch_id_path(
    g: &CompactDbgGraph,
    path: &[NodeId],
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
) -> Result<Vec<u8>, GraphError> {
    if path.is_empty() {
        return Ok(Vec::new());
    }
    if path.len() == 1 {
        return pick_forward_id(g, path[0], forward);
    }
    let f0 = pick_forward_id(g, path[0], forward)?;
    let rc0 = reverse_complement(&f0);
    let mut starters: Vec<Vec<u8>> = if f0 == rc0 { vec![f0] } else { vec![f0, rc0] };
    starters.sort_by(|a, b| cmp_dna(a, b));
    starters.dedup();

    let mut best: Option<Vec<u8>> = None;
    for s0 in starters {
        if let Ok(seq) = stitch_from_ids(g, path, forward, k, s0) {
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

fn stitch_from_ids(
    g: &CompactDbgGraph,
    path: &[NodeId],
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    mut cur: Vec<u8>,
) -> Result<Vec<u8>, GraphError> {
    let mut seq = Vec::with_capacity(k + path.len().saturating_sub(1));
    seq.extend_from_slice(&cur);
    for node in &path[1..] {
        let cand = pick_forward_id(g, *node, forward)?;
        let direct_ok = cur[1..] == cand[..k - 1];
        let rc_c = reverse_complement(&cand);
        let rc_ok = cand != rc_c && cur[1..] == rc_c[..k - 1];
        let use_c = match (direct_ok, rc_ok) {
            (true, false) => cand,
            (false, true) => rc_c,
            (true, true) => {
                if cmp_dna(&cand, &rc_c) == Ordering::Greater {
                    rc_c
                } else {
                    cand
                }
            }
            (false, false) => return Err(GraphError::OrientationConflict),
        };
        if let Some(&last) = use_c.last() {
            seq.push(last);
        }
        cur = use_c;
    }
    Ok(seq)
}

fn pick_forward_id(
    g: &CompactDbgGraph,
    id: NodeId,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
) -> Result<Vec<u8>, GraphError> {
    forward
        .get(g.node_name(id))
        .cloned()
        .ok_or(GraphError::OrientationConflict)
}

/// One **reference contig** path per connected component (**Phase-1 disconnected graph policy**).
pub fn reference_contig_paths(
    g: &DbgGraph,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    tie_break: ContigWalkTieBreak,
) -> Result<Vec<Vec<Vec<u8>>>, GraphError> {
    let compact = g.compact_view();
    let mut out: Vec<Vec<Vec<u8>>> = Vec::new();
    for comp in connected_component_ids(&compact) {
        if comp.is_empty() {
            continue;
        }
        let best = if let Some(path) = linear_component_path_ids(&compact, &comp) {
            let mut rev = path.clone();
            rev.reverse();
            let score = walk_edge_score(&compact, &path);
            pick_best_stitchable_id_path(&compact, forward, k, [(score, path), (score, rev)])?
        } else {
            let candidates = comp
                .iter()
                .map(|seed| greedy_simple_path(&compact, *seed, tie_break));
            pick_best_stitchable_id_path(&compact, forward, k, candidates)?
        };
        if let Some(path) = best {
            out.push(id_path_to_names(&compact, &path));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::BTreeMap;

    use crate::dbg::stitch_sequence;

    fn graph_with_nodes(k: usize, nodes: &[&[u8]]) -> DbgGraph {
        DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            for node in nodes {
                mm.insert((*node).to_vec(), 1);
            }
            mm
        })
    }

    #[test]
    fn linear_component_path_walks_once_from_endpoint() {
        let a = b"AAAA";
        let b = b"AAAC";
        let c = b"AACC";
        let mut g = graph_with_nodes(4, &[a, b, c]);
        g.add_undirected_edge(a, b, 1).unwrap();
        g.add_undirected_edge(b, c, 1).unwrap();
        let compact = g.compact_view();
        let comp = connected_component_ids(&compact)
            .into_iter()
            .next()
            .expect("component");

        let path = linear_component_path_ids(&compact, &comp).expect("linear path");
        assert_eq!(
            id_path_to_names(&compact, &path),
            vec![a.to_vec(), b.to_vec(), c.to_vec()]
        );
    }

    #[test]
    fn reference_paths_keep_linear_component_sequence() {
        let a = b"AAAA";
        let b = b"AAAC";
        let c = b"AACC";
        let mut g = graph_with_nodes(4, &[a, b, c]);
        g.add_undirected_edge(a, b, 1).unwrap();
        g.add_undirected_edge(b, c, 1).unwrap();
        let mut forward = HashMap::new();
        for node in [a, b, c] {
            forward.insert(node.to_vec(), node.to_vec());
        }

        let paths = reference_contig_paths(&g, &forward, 4, ContigWalkTieBreak::Phase1Lex)
            .expect("reference paths");
        assert_eq!(paths.len(), 1);
        let seq = stitch_sequence(&paths[0], &forward, 4).expect("stitched");
        assert_eq!(seq, b"AAAACC");
    }

    #[test]
    fn best_stitchable_path_ignores_unstitchable_higher_score() {
        let a = b"AAA";
        let b = b"AAC";
        let c = b"CCC";
        let mut g = graph_with_nodes(3, &[a, b, c]);
        g.add_undirected_edge(a, b, 1).unwrap();
        g.add_undirected_edge(a, c, 10).unwrap();
        let mut forward = HashMap::new();
        for node in [a, b, c] {
            forward.insert(node.to_vec(), node.to_vec());
        }

        let candidates = vec![
            (10, vec![a.to_vec(), c.to_vec()]),
            (1, vec![a.to_vec(), b.to_vec()]),
        ];
        let path = pick_best_stitchable_path(&forward, 3, candidates)
            .expect("search succeeds")
            .expect("lower score candidate is stitchable");

        assert_eq!(path, vec![a.to_vec(), b.to_vec()]);
    }
}
