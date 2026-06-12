//! **Phase-1 graph simplification**: tip clipping plus **bounded diamond bubbles** (two
//! internally vertex-disjoint length-2 paths between opposite corners), resolved using **read-derived
//! edge multiplicities** with deterministic lex tie-breaks on branch vertices (**Phase-1 bubble
//! resolution** / **Phase-1 bubble bounds**).

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::dbg::graph::DbgGraph;
use crate::kmer::cmp_dna;

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SimplifyDecisionAction {
    RemoveTipEdge,
    RemoveDiamondBranch,
    RetainDiploidDiamond,
    SkipAmbiguousK22Diamond,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SimplifyDecision {
    pub action: SimplifyDecisionAction,
    pub reason: String,
    pub nodes: Vec<String>,
    pub removed_node: Option<String>,
    pub removed_edge: Option<[String; 2]>,
    pub score_a: Option<u64>,
    pub score_b: Option<u64>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SimplifyDecisionLog {
    pub tips: Vec<SimplifyDecision>,
    pub diamonds: Vec<SimplifyDecision>,
}

/// Decision counters from automatic graph simplification.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SimplifyStats {
    pub tips_removed: usize,
    pub diamond_bubbles_resolved: usize,
    pub diploid_diamonds_retained: usize,
    pub ambiguous_k22_diamonds_skipped: usize,
}

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
            max_bubble_vertices: (2 * k).clamp(8, 64),
            max_bubble_internal_bases: (8 * k).max(32),
        }
    }
}

/// Iteratively remove short low-coverage tips (degree-1 leaves).
pub fn remove_tips(graph: &mut DbgGraph, p: &SimplifyParams) -> usize {
    remove_tips_with_decisions(graph, p).0
}

pub fn remove_tips_with_decisions(
    graph: &mut DbgGraph,
    p: &SimplifyParams,
) -> (usize, Vec<SimplifyDecision>) {
    let mut removed = 0usize;
    let mut decisions = Vec::new();
    loop {
        let planned = plan_tip_clips(graph, p);
        if planned.is_empty() {
            break;
        }
        let mut removed_any = false;
        for decision in planned {
            let Some(edge) = decision.removed_edge.as_ref() else {
                continue;
            };
            let leaf = edge[0].as_bytes();
            let neigh = edge[1].as_bytes();
            if graph.degree(leaf) == 1 && graph.adj.get(leaf).is_some_and(|m| m.contains_key(neigh))
            {
                graph.remove_undirected_edge(leaf, neigh);
                removed += 1;
                removed_any = true;
                decisions.push(decision);
            }
        }
        if !removed_any {
            break;
        }
    }
    (removed, decisions)
}

pub fn plan_tip_clips(graph: &DbgGraph, p: &SimplifyParams) -> Vec<SimplifyDecision> {
    let leaves: Vec<Vec<u8>> = graph
        .adj
        .keys()
        .filter(|u| graph.degree(u) == 1)
        .cloned()
        .collect();
    let mut decisions = Vec::new();
    for leaf in leaves {
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
        if let Some(chain_bases) =
            tip_chain_bases_to_junction(graph, &leaf, &neigh, p.max_tip_bases, graph.k)
        {
            decisions.push(SimplifyDecision {
                action: SimplifyDecisionAction::RemoveTipEdge,
                reason: format!(
                    "tip leaf multiplicity {mul} <= {} and chain span {chain_bases} <= {}",
                    p.tip_max_multiplicity, p.max_tip_bases
                ),
                nodes: vec![node_label(&leaf), node_label(&neigh)],
                removed_node: Some(node_label(&leaf)),
                removed_edge: Some([node_label(&leaf), node_label(&neigh)]),
                score_a: None,
                score_b: None,
            });
        }
    }
    decisions
}

fn tip_chain_bases_to_junction(
    graph: &DbgGraph,
    leaf: &[u8],
    nbr: &[u8],
    cap: usize,
    k: usize,
) -> Option<usize> {
    let mut len = k;
    let mut prev = leaf.to_vec();
    let mut cur = nbr.to_vec();
    let mut visited = BTreeSet::from([prev.clone(), cur.clone()]);
    loop {
        match graph.degree(&cur) {
            0 | 1 => return None,
            2 => {}
            _ => return Some(len),
        }
        if len >= cap {
            return None;
        }
        let nexts: Vec<Vec<u8>> = graph
            .adj
            .get(&cur)
            .into_iter()
            .flat_map(|m| m.keys())
            .filter(|x| *x != &prev)
            .cloned()
            .collect();
        if nexts.len() != 1 {
            return None;
        }
        let nxt = nexts.into_iter().next().unwrap();
        if visited.contains(&nxt) {
            return None;
        }
        visited.insert(nxt.clone());
        prev = cur;
        cur = nxt;
        len += 1;
    }
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
    g.adj.get(u).and_then(|m| m.get(v)).copied().unwrap_or(0)
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
) -> SimplifyStats {
    remove_diamond_bubbles_ext_with_decisions(graph, p, diploid).0
}

pub fn remove_diamond_bubbles_ext_with_decisions(
    graph: &mut DbgGraph,
    p: &SimplifyParams,
    diploid: Option<DiploidSimplifyMode>,
) -> (SimplifyStats, Vec<SimplifyDecision>) {
    let mut stats = SimplifyStats::default();
    let mut retained_diploid_motifs: BTreeSet<Vec<Vec<u8>>> = BTreeSet::new();
    let mut skipped_k22_motifs: BTreeSet<Vec<Vec<u8>>> = BTreeSet::new();
    let mut decisions = Vec::new();
    if p.max_bubble_vertices < 4 {
        return (stats, decisions);
    }
    let est_bases = graph.k.saturating_mul(3);
    if est_bases > p.max_bubble_internal_bases {
        return (stats, decisions);
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
                    let mut motif = vec![u.clone(), a.clone(), b.clone(), m.clone()];
                    motif.sort_by(|x, y| cmp_dna(x, y));
                    if skipped_k22_motifs.insert(motif) {
                        stats.ambiguous_k22_diamonds_skipped += 1;
                        decisions.push(SimplifyDecision {
                            action: SimplifyDecisionAction::SkipAmbiguousK22Diamond,
                            reason: "pure K2,2 diamond has two valid opposite-corner labelings"
                                .to_string(),
                            nodes: sorted_node_labels([&u, a, b, &m]),
                            removed_node: None,
                            removed_edge: None,
                            score_a: None,
                            score_b: None,
                        });
                    }
                    continue;
                }
                if nodes.len() > p.max_bubble_vertices {
                    continue;
                }

                let s_a = branch_score(graph, &u, a, &m);
                let s_b = branch_score(graph, &u, b, &m);
                if diploid.is_some() && branch_scores_nearly_balanced(s_a, s_b, 5) {
                    let mut motif = vec![u.clone(), a.clone(), b.clone(), m.clone()];
                    motif.sort_by(|x, y| cmp_dna(x, y));
                    if retained_diploid_motifs.insert(motif) {
                        stats.diploid_diamonds_retained += 1;
                        decisions.push(SimplifyDecision {
                            action: SimplifyDecisionAction::RetainDiploidDiamond,
                            reason: "diploid mode retained near-balanced diamond branches"
                                .to_string(),
                            nodes: sorted_node_labels([&u, a, b, &m]),
                            removed_node: None,
                            removed_edge: None,
                            score_a: Some(s_a),
                            score_b: Some(s_b),
                        });
                    }
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

                let kept = if remove_b { a.clone() } else { b.clone() };
                decisions.push(SimplifyDecision {
                    action: SimplifyDecisionAction::RemoveDiamondBranch,
                    reason: "lower edge-support diamond branch removed".to_string(),
                    nodes: sorted_node_labels([&u, a, b, &m]),
                    removed_node: Some(node_label(&drop_x)),
                    removed_edge: Some([node_label(&drop_u), node_label(&drop_x)]),
                    score_a: Some(s_a),
                    score_b: Some(s_b),
                });
                debug_assert_ne!(kept, drop_x);
                action = Some((drop_u, drop_x));
                break 'outer;
            }
        }

        let Some((drop_u, drop_x)) = action else {
            break;
        };
        graph.remove_undirected_edge(&drop_u, &drop_x);
        graph.remove_vertex_from_adj(&drop_x);
        stats.diamond_bubbles_resolved += 1;
    }
    (stats, decisions)
}

/// Resolve diamond bubbles using **Phase-1** rules only (collapse every resolvable diamond).
pub fn remove_diamond_bubbles(graph: &mut DbgGraph, p: &SimplifyParams) -> SimplifyStats {
    remove_diamond_bubbles_ext(graph, p, None)
}

fn node_label(node: &[u8]) -> String {
    String::from_utf8_lossy(node).into_owned()
}

fn sorted_node_labels<const N: usize>(nodes: [&Vec<u8>; N]) -> Vec<String> {
    let mut out: Vec<Vec<u8>> = nodes.into_iter().cloned().collect();
    out.sort_by(|x, y| cmp_dna(x, y));
    out.into_iter().map(|node| node_label(&node)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn graph_with_nodes(k: usize, nodes: &[&[u8]], mul: u64) -> DbgGraph {
        DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            for v in nodes {
                mm.insert((*v).to_vec(), mul);
            }
            mm
        })
    }

    #[test]
    fn tip_clipping_retains_isolated_linear_component() {
        let k = 4usize;
        let a = b"AAAA";
        let b = b"AAAC";
        let c = b"AACC";
        let mut g = graph_with_nodes(k, &[a, b, c], 1);
        g.add_undirected_edge(a, b, 1).unwrap();
        g.add_undirected_edge(b, c, 1).unwrap();

        let p = SimplifyParams {
            max_tip_bases: 12,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 1000,
        };
        let removed = remove_tips(&mut g, &p);

        assert_eq!(removed, 0);
        assert_eq!(g.degree(a), 1);
        assert_eq!(g.degree(b), 2);
        assert_eq!(g.degree(c), 1);
    }

    #[test]
    fn tip_clipping_removes_short_leaf_attached_to_junction() {
        let k = 4usize;
        let tip = b"AAAT";
        let hub = b"AAAC";
        let keep_a = b"AACC";
        let keep_b = b"ACCC";
        let mut g = DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            mm.insert(tip.to_vec(), 1);
            mm.insert(hub.to_vec(), 10);
            mm.insert(keep_a.to_vec(), 10);
            mm.insert(keep_b.to_vec(), 10);
            mm
        });
        g.add_undirected_edge(tip, hub, 1).unwrap();
        g.add_undirected_edge(hub, keep_a, 10).unwrap();
        g.add_undirected_edge(hub, keep_b, 10).unwrap();

        let p = SimplifyParams {
            max_tip_bases: 12,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 1000,
        };
        let removed = remove_tips(&mut g, &p);

        assert_eq!(removed, 1);
        assert_eq!(g.degree(tip), 0);
        assert_eq!(g.degree(hub), 2);
        assert!(g.adj[hub.as_slice()].contains_key(keep_a.as_slice()));
        assert!(g.adj[hub.as_slice()].contains_key(keep_b.as_slice()));
    }

    #[test]
    fn tip_decision_matches_tip_mutation() {
        let k = 4usize;
        let tip = b"AAAT";
        let hub = b"AAAC";
        let keep_a = b"AACC";
        let keep_b = b"ACCC";
        let mut g = DbgGraph::new(k, {
            let mut mm = BTreeMap::new();
            mm.insert(tip.to_vec(), 1);
            mm.insert(hub.to_vec(), 10);
            mm.insert(keep_a.to_vec(), 10);
            mm.insert(keep_b.to_vec(), 10);
            mm
        });
        g.add_undirected_edge(tip, hub, 1).unwrap();
        g.add_undirected_edge(hub, keep_a, 10).unwrap();
        g.add_undirected_edge(hub, keep_b, 10).unwrap();

        let p = SimplifyParams {
            max_tip_bases: 12,
            tip_max_multiplicity: 2,
            max_bubble_vertices: 16,
            max_bubble_internal_bases: 1000,
        };
        let planned = plan_tip_clips(&g, &p);
        let (removed, decisions) = remove_tips_with_decisions(&mut g, &p);

        assert_eq!(planned.len(), 1);
        assert_eq!(removed, 1);
        assert_eq!(decisions, planned);
        assert_eq!(decisions[0].action, SimplifyDecisionAction::RemoveTipEdge);
        assert_eq!(decisions[0].removed_node.as_deref(), Some("AAAT"));
        assert_eq!(g.degree(tip), 0);
    }

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
        let stats = remove_diamond_bubbles(&mut g, &p);

        assert_eq!(stats.diamond_bubbles_resolved, 1);
        assert!(!g.adj.contains_key(&b) || g.degree(&b) == 0);
        assert!(g.adj.contains_key(&a));
    }

    #[test]
    fn diamond_decision_matches_branch_mutation() {
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
        let (stats, decisions) = remove_diamond_bubbles_ext_with_decisions(&mut g, &p, None);

        assert_eq!(stats.diamond_bubbles_resolved, 1);
        assert_eq!(decisions.len(), 1);
        assert_eq!(
            decisions[0].action,
            SimplifyDecisionAction::RemoveDiamondBranch
        );
        assert_eq!(decisions[0].removed_node.as_deref(), Some("AATA"));
        assert_eq!(decisions[0].score_a, Some(10));
        assert_eq!(decisions[0].score_b, Some(2));
        assert!(!g.adj.contains_key(&b) || g.degree(&b) == 0);
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
        let stats = remove_diamond_bubbles_ext(&mut g, &p, Some(DiploidSimplifyMode));

        assert_eq!(stats.diploid_diamonds_retained, 1);
        assert_eq!(stats.diamond_bubbles_resolved, 0);
        assert!(g.adj.contains_key(&b));
        assert!(g.adj.contains_key(&a));
    }

    #[test]
    fn diploid_retention_has_decision_reason() {
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
        let (stats, decisions) =
            remove_diamond_bubbles_ext_with_decisions(&mut g, &p, Some(DiploidSimplifyMode));

        assert_eq!(stats.diploid_diamonds_retained, 1);
        assert_eq!(decisions.len(), 1);
        assert_eq!(
            decisions[0].action,
            SimplifyDecisionAction::RetainDiploidDiamond
        );
        assert!(decisions[0].removed_node.is_none());
        assert!(g.adj.contains_key(&b));
        assert!(g.adj.contains_key(&a));
    }
}
