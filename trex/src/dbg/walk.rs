//! **Phase-1 reference contig walker**: one **vertex-simple** path per connected component,
//! ranked by total **read-derived edge multiplicity** along arcs; ties break on **lexicographically
//! smallest** stitched sequence (**Phase-1 contig walk score** + **Phase-1 contig tie-break**).
//!
//! The search is **deterministic greedy extension** from every vertex seed in the component,
//! keeping the best scoring path found.

use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::dbg::graph::DbgGraph;
use crate::dbg::unitig::stitch_sequence;
use crate::error::GraphError;
use crate::kmer::cmp_dna;

/// Tie-breaking when two forward neighbors share the same **read-edge multiplicity** from `cur`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ContigWalkTieBreak {
    /// **Phase-1 contig tie-break**: lexicographically smallest neighbor *k*-mer.
    #[default]
    Phase1Lex,
    /// **Phase-2 Illumina (experimental)**: higher **trusted vertex multiplicity** on the neighbor, then lex smallest.
    Phase2DiploidNodeMul,
}

fn edge_weight(g: &DbgGraph, u: &[u8], v: &[u8]) -> u64 {
    g.adj.get(u).and_then(|m| m.get(v)).copied().unwrap_or(0)
}

fn walk_edge_score(g: &DbgGraph, path: &[Vec<u8>]) -> u64 {
    if path.len() < 2 {
        return 0;
    }
    let mut s = 0u64;
    for w in path.windows(2) {
        s += edge_weight(g, &w[0], &w[1]);
    }
    s
}

fn neighbor_vertex_mul(g: &DbgGraph, nb: &[u8]) -> u64 {
    g.node_mul.get(nb).copied().unwrap_or(0)
}

fn pick_best_neighbor(
    g: &DbgGraph,
    cur: &[u8],
    forbidden: &HashSet<Vec<u8>>,
    tie_break: ContigWalkTieBreak,
) -> Option<Vec<u8>> {
    let mut best: Option<(u64, Vec<u8>)> = None;
    for (nb, &w) in g.adj.get(cur).into_iter().flat_map(|m| m.iter()) {
        if forbidden.contains(nb) {
            continue;
        }
        match &best {
            None => best = Some((w, nb.clone())),
            Some((bw, bnb)) => {
                let take = match w.cmp(bw) {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal => match tie_break {
                        ContigWalkTieBreak::Phase1Lex => cmp_dna(nb, bnb) == Ordering::Less,
                        ContigWalkTieBreak::Phase2DiploidNodeMul => {
                            let mn = neighbor_vertex_mul(g, nb);
                            let mb = neighbor_vertex_mul(g, bnb);
                            mn > mb || (mn == mb && cmp_dna(nb, bnb) == Ordering::Less)
                        }
                    },
                };
                if take {
                    best = Some((w, nb.clone()));
                }
            }
        }
    }
    best.map(|(_, nb)| nb)
}

/// Greedy **vertex-simple** path: extend forward from `seed`, then backward from `seed` without
/// reusing vertices from the forward segment (except `seed`).
fn greedy_simple_path(g: &DbgGraph, seed: &[u8], tie_break: ContigWalkTieBreak) -> Vec<Vec<u8>> {
    let mut forward_forbidden: HashSet<Vec<u8>> = HashSet::new();
    forward_forbidden.insert(seed.to_vec());

    let mut forward: Vec<Vec<u8>> = vec![seed.to_vec()];
    let mut cur = seed.to_vec();
    while let Some(nb) = pick_best_neighbor(g, &cur, &forward_forbidden, tie_break) {
        forward_forbidden.insert(nb.clone());
        forward.push(nb.clone());
        cur = nb;
    }

    let mut all_forbidden = forward_forbidden.clone();
    let mut back_segments: Vec<Vec<u8>> = Vec::new();
    cur = seed.to_vec();
    while let Some(nb) = pick_best_neighbor(g, &cur, &all_forbidden, tie_break) {
        all_forbidden.insert(nb.clone());
        back_segments.push(nb.clone());
        cur = nb;
    }

    let mut path: Vec<Vec<u8>> = back_segments.into_iter().rev().collect();
    path.extend(forward);
    path
}

pub fn connected_components(g: &DbgGraph) -> Vec<Vec<Vec<u8>>> {
    let mut unseen: BTreeSet<Vec<u8>> = g.adj.keys().cloned().collect();
    let mut comps: Vec<Vec<Vec<u8>>> = Vec::new();
    while let Some(seed) = unseen.iter().next().cloned() {
        unseen.remove(&seed);
        let mut stack = vec![seed];
        let mut comp: Vec<Vec<u8>> = Vec::new();
        while let Some(u) = stack.pop() {
            comp.push(u.clone());
            if let Some(neigh) = g.adj.get(&u) {
                for v in neigh.keys() {
                    if unseen.remove(v) {
                        stack.push(v.clone());
                    }
                }
            }
        }
        comp.sort_unstable_by(|a, b| cmp_dna(a, b));
        comps.push(comp);
    }
    comps.sort_unstable_by(|a, b| cmp_dna(a[0].as_slice(), b[0].as_slice()));
    comps
}

fn linear_component_path(g: &DbgGraph, comp: &[Vec<u8>]) -> Option<Vec<Vec<u8>>> {
    if comp.is_empty() || comp.iter().any(|u| g.degree(u) > 2) {
        return None;
    }

    let start = comp
        .iter()
        .filter(|u| g.degree(u) <= 1)
        .min_by(|a, b| cmp_dna(a, b))
        .or_else(|| comp.iter().min_by(|a, b| cmp_dna(a, b)))?
        .clone();

    let mut path = Vec::with_capacity(comp.len());
    let mut seen: HashSet<Vec<u8>> = HashSet::with_capacity(comp.len());
    let mut prev: Option<Vec<u8>> = None;
    let mut cur = start;

    loop {
        if !seen.insert(cur.clone()) {
            break;
        }
        path.push(cur.clone());

        let next = g.adj.get(&cur).and_then(|neigh| {
            neigh
                .keys()
                .filter(|nb| prev.as_ref().map(|p| *nb != p).unwrap_or(true))
                .find(|nb| !seen.contains(*nb))
                .cloned()
        });
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

fn pick_best_stitchable_path(
    g: &DbgGraph,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    candidates: Vec<Vec<Vec<u8>>>,
) -> Result<Option<Vec<Vec<u8>>>, GraphError> {
    let mut best: Option<(u64, Vec<u8>, Vec<Vec<u8>>)> = None;
    for path in candidates {
        let score = walk_edge_score(g, &path);
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

/// One **reference contig** path per connected component (**Phase-1 disconnected graph policy**).
pub fn reference_contig_paths(
    g: &DbgGraph,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    tie_break: ContigWalkTieBreak,
) -> Result<Vec<Vec<Vec<u8>>>, GraphError> {
    let mut out: Vec<Vec<Vec<u8>>> = Vec::new();
    for comp in connected_components(g) {
        if comp.is_empty() {
            continue;
        }
        let best = if let Some(path) = linear_component_path(g, &comp) {
            let mut rev = path.clone();
            rev.reverse();
            pick_best_stitchable_path(g, forward, k, vec![path, rev])?
        } else {
            let candidates = comp
                .iter()
                .map(|seed| greedy_simple_path(g, seed, tie_break))
                .collect();
            pick_best_stitchable_path(g, forward, k, candidates)?
        };
        if let Some(path) = best {
            out.push(path);
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
        let comp = vec![a.to_vec(), b.to_vec(), c.to_vec()];

        let path = linear_component_path(&g, &comp).expect("linear path");
        assert_eq!(path, comp);
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
}
