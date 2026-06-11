//! **Phase-1 graph simplification**: tip clipping plus **bounded diamond bubbles** (two
//! internally vertex-disjoint length-2 paths between opposite corners), resolved using **read-derived
//! edge multiplicities** with deterministic lex tie-breaks on branch vertices (**Phase-1 bubble
//! resolution** / **Phase-1 bubble bounds**).

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::dbg::graph::DbgGraph;
use crate::kmer::cmp_dna;

#[derive(Debug, Clone)]
pub struct SimplifyParams {
    /// Maximum **sequence** length (bases) of a tip chain to clip.
    pub max_tip_bases: usize,
    /// Remove tip leaf if its **trusted** multiplicity is **≤** this floor.
    pub tip_max_multiplicity: u64,
    /// Maximum distinct vertices touched by an automatic bubble motif (including endpoints).
    pub max_bubble_vertices: usize,
    /// Conservative **sequence-span budget** (bases) for automatic bubble resolution.
    pub max_bubble_internal_bases: usize,
}

impl Default for SimplifyParams {
    fn default() -> Self {
        Self {
            max_tip_bases: 24,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 96,
        }
    }
}

impl SimplifyParams {
    pub fn for_k(k: usize) -> Self {
        Self {
            max_tip_bases: (3 * k).max(8),
            tip_max_multiplicity: 2,
            max_bubble_vertices: (2 * k).max(8).min(64),
            max_bubble_internal_bases: (8 * k).max(32),
        }
    }
}

/// Iteratively remove short low-coverage tips (degree-1 leaves).
pub fn remove_tips(graph: &mut DbgGraph, p: &SimplifyParams) {
    loop {
        let leaves: Vec<Vec<u8>> = graph
            .adj
            .keys()
            .filter(|u| graph.degree(u) == 1)
            .cloned()
            .collect();
        if leaves.is_empty() {
            break;
        }
        let mut removed_any = false;
        for leaf in leaves {
            if graph.degree(&leaf) != 1 {
                continue;
            }
            let mul = *graph.node_mul.get(&leaf).unwrap_or(&0);
            if mul > p.tip_max_multiplicity {
                continue;
            }
            let neigh: Vec<u8> = graph
                .adj
                .get(&leaf)
                .and_then(|m| m.keys().next())
                .cloned()
                .unwrap_or_default();
            if neigh.is_empty() {
                continue;
            }
            let tip_len = tip_chain_bases(graph, &leaf, &neigh, p.max_tip_bases, graph.k);
            if tip_len <= p.max_tip_bases {
                graph.remove_undirected_edge(&leaf, &neigh);
                removed_any = true;
            }
        }
        if !removed_any {
            break;
        }
    }
}

fn tip_chain_bases(graph: &DbgGraph, leaf: &[u8], nbr: &[u8], cap: usize, k: usize) -> usize {
    let mut len = k;
    let mut prev = leaf.to_vec();
    let mut cur = nbr.to_vec();
    let mut visited = BTreeSet::from([prev.clone(), cur.clone()]);
    while len < cap && graph.degree(&cur) == 2 {
        let nexts: Vec<Vec<u8>> = graph
            .adj
            .get(&cur)
            .into_iter()
            .flat_map(|m| m.keys())
            .filter(|x| *x != &prev)
            .cloned()
            .collect();
        if nexts.len() != 1 {
            break;
        }
        let nxt = nexts.into_iter().next().unwrap();
        if visited.contains(&nxt) {
            break;
        }
        visited.insert(nxt.clone());
        prev = cur;
        cur = nxt;
        len += 1;
    }
    len
}

/// **Phase-1 simplified graph invariants:** forbid self-adjacency.
pub fn assert_no_self_loops(graph: &DbgGraph) -> Result<(), crate::error::GraphError> {
    for (u, neigh) in &graph.adj {
        if neigh.contains_key(u) {
            return Err(crate::error::GraphError::SimplifiedSelfLoop);
        }
    }
    Ok(())
}

fn edge_weight(g: &DbgGraph, u: &[u8], v: &[u8]) -> u64 {
    g.adj
        .get(u)
        .and_then(|m| m.get(v))
        .copied()
        .unwrap_or(0)
}

fn branch_score(g: &DbgGraph, u: &[u8], x: &[u8], m: &[u8]) -> u64 {
    edge_weight(g, u, x).saturating_add(edge_weight(g, x, m))
}

/// Marker for **Phase-2 Illumina diploid** diamond handling: balanced branches are left intact.
#[derive(Debug, Clone, Copy)]
pub struct DiploidSimplifyMode;

fn branch_scores_nearly_balanced(s_a: u64, s_b: u64, max_relative_diff_percent: u64) -> bool {
    if s_a == s_b {
        return true;
    }
    let hi = s_a.max(s_b);
    let lo = s_a.min(s_b);
    if hi == 0 {
        return true;
    }
    (hi - lo).saturating_mul(100) <= hi.saturating_mul(max_relative_diff_percent)
}

/// Resolve **diamond** bubbles `u–a–m` vs `u–b–m` when both `a` and `b` are degree-2 junctions and
/// the motif fits **Phase-1 bubble bounds**. Lower-scoring branch (read-edge support) is removed.
///
/// When `diploid` is **`Some`**, branches with **equal** or **near-equal** read-edge support (within
/// 5% of the stronger branch) are **not** collapsed so heterozygous structure can remain.
pub fn remove_diamond_bubbles_ext(
    graph: &mut DbgGraph,
    p: &SimplifyParams,
    diploid: Option<DiploidSimplifyMode>,
) {
    if p.max_bubble_vertices < 4 {
        return;
    }
    let est_bases = graph.k.saturating_mul(3);
    if est_bases > p.max_bubble_internal_bases {
        return;
    }

    loop {
        let verts: Vec<Vec<u8>> = graph.adj.keys().cloned().collect();
        let mut action: Option<(Vec<u8>, Vec<u8>)> = None;

        'outer: for a in &verts {
            if graph.degree(a) != 2 {
                continue;
            }
            let nbrs: Vec<Vec<u8>> = graph.adj[a].keys().cloned().collect();
            if nbrs.len() != 2 {
                continue;
            }
            let u = nbrs[0].clone();
            let m = nbrs[1].clone();

            for b in graph.adj.get(&m).into_iter().flat_map(|mp| mp.keys()) {
                if b == a || b == &u {
                    continue;
                }
                if graph.degree(b) != 2 {
                    continue;
                }
                let bn: Vec<Vec<u8>> = graph.adj[b].keys().cloned().collect();
                if bn.len() != 2 {
                    continue;
                }
                let has_u = bn.iter().any(|x| x == &u);
                let has_m = bn.iter().any(|x| x == &m);
                if !has_u || !has_m {
                    continue;
                }

                let mut nodes = BTreeSet::new();
                nodes.insert(a.clone());
                nodes.insert(b.clone());
                nodes.insert(u.clone());
                nodes.insert(m.clone());
                if nodes.len() != 4 {
                    continue;
                }
                if graph.degree(a) == 2
                    && graph.degree(b) == 2
                    && graph.degree(&u) == 2
                    && graph.degree(&m) == 2
                {
                    // Pure **K₂,₂**: two valid opposite-corner labelings; skip automatic resolution.
                    continue;
                }
                if nodes.len() > p.max_bubble_vertices {
                    continue;
                }

                let s_a = branch_score(graph, &u, a, &m);
                let s_b = branch_score(graph, &u, b, &m);
                if diploid.is_some() && branch_scores_nearly_balanced(s_a, s_b, 5) {
                    continue;
                }
                let remove_b = match s_a.cmp(&s_b) {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal => cmp_dna(a, b) == Ordering::Less,
                };

                let (drop_u, drop_x) = if remove_b {
                    (u.clone(), b.clone())
                } else {
                    (u.clone(), a.clone())
                };

                action = Some((drop_u, drop_x));
                break 'outer;
            }
        }

        let Some((drop_u, drop_x)) = action else {
            break;
        };
        graph.remove_undirected_edge(&drop_u, &drop_x);
        graph.remove_vertex_from_adj(&drop_x);
    }
}

/// Resolve diamond bubbles using **Phase-1** rules only (collapse every resolvable diamond).
pub fn remove_diamond_bubbles(graph: &mut DbgGraph, p: &SimplifyParams) {
    remove_diamond_bubbles_ext(graph, p, None);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn diamond_removes_lower_branch() {
        let k = 4usize;
        let u = b"AAAA".to_vec();
        let a = b"AAAT".to_vec();
        let m = b"AATT".to_vec();
        let b = b"AATA".to_vec();
        // Leaf on `u` so `u` is not degree-2; avoids ambiguous K₂,₂ corner labeling when scanning `a`.
        let leaf = b"AAAC".to_vec();

        let mut g = DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            for v in [&u, &a, &m, &b, &leaf] {
                mm.insert(v.clone(), 10);
            }
            mm
        });

        g.add_undirected_edge(&u, &a, 5).unwrap();
        g.add_undirected_edge(&a, &m, 5).unwrap();
        g.add_undirected_edge(&u, &b, 1).unwrap();
        g.add_undirected_edge(&b, &m, 1).unwrap();
        g.add_undirected_edge(&u, &leaf, 1).unwrap();

        let p = SimplifyParams {
            max_tip_bases: 8,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 1000,
        };
        remove_diamond_bubbles(&mut g, &p);

        assert!(g.adj.get(&b).is_none() || g.degree(&b) == 0);
        assert!(g.adj.contains_key(&a));
    }

    #[test]
    fn diamond_diploid_retains_near_balanced_branches() {
        let k = 4usize;
        let u = b"AAAA".to_vec();
        let a = b"AAAT".to_vec();
        let m = b"AATT".to_vec();
        let b = b"AATA".to_vec();
        let leaf = b"AAAC".to_vec();

        let mut g = DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            for v in [&u, &a, &m, &b, &leaf] {
                mm.insert(v.clone(), 10);
            }
            mm
        });

        g.add_undirected_edge(&u, &a, 10).unwrap();
        g.add_undirected_edge(&a, &m, 10).unwrap();
        g.add_undirected_edge(&u, &b, 10).unwrap();
        g.add_undirected_edge(&b, &m, 9).unwrap();
        g.add_undirected_edge(&u, &leaf, 1).unwrap();

        let p = SimplifyParams {
            max_tip_bases: 8,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 1000,
        };
        remove_diamond_bubbles_ext(&mut g, &p, Some(DiploidSimplifyMode));

        assert!(g.adj.contains_key(&b));
        assert!(g.adj.contains_key(&a));
    }
}
