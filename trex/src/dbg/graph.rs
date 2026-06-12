//! Build a **trusted** de Bruijn graph: nodes are **canonical** *k*-mers; edges follow read-adjacent
//! forward windows whose overlap is **k − 1** (**Phase-1 k-mer identity** + counting orientation).

use std::collections::BTreeMap;
use std::collections::HashMap;

use crate::error::GraphError;
use crate::illumina::preprocess::n_free_acgt_segments;
use crate::illumina::read::Read;
use crate::kmer::canonical_kmer;

/// Dense internal vertex id for read-only graph walks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct NodeId(pub(crate) usize);

/// Read-only compact view over a [`DbgGraph`], keyed by dense [`NodeId`] values.
///
/// Public graph names remain canonical *k*-mer byte strings; this view is a hot-path adapter for
/// traversal code that would otherwise repeat ordered-map lookups and clone keys during walks.
#[derive(Debug, Clone)]
pub(crate) struct CompactDbgGraph {
    node_names: Vec<Vec<u8>>,
    node_mul: Vec<u64>,
    adj: Vec<Vec<(NodeId, u64)>>,
}

impl CompactDbgGraph {
    pub(crate) fn node_ids(&self) -> impl Iterator<Item = NodeId> + '_ {
        (0..self.node_names.len()).map(NodeId)
    }

    pub(crate) fn node_name(&self, id: NodeId) -> &[u8] {
        &self.node_names[id.0]
    }

    pub(crate) fn node_name_vec(&self, id: NodeId) -> Vec<u8> {
        self.node_names[id.0].clone()
    }

    pub(crate) fn node_mul(&self, id: NodeId) -> u64 {
        self.node_mul[id.0]
    }

    pub(crate) fn degree(&self, id: NodeId) -> usize {
        self.adj[id.0].len()
    }

    pub(crate) fn neighbors(&self, id: NodeId) -> &[(NodeId, u64)] {
        &self.adj[id.0]
    }

    pub(crate) fn edge_weight(&self, a: NodeId, b: NodeId) -> u64 {
        self.adj[a.0]
            .iter()
            .find_map(|(nb, weight)| (*nb == b).then_some(*weight))
            .unwrap_or(0)
    }
}

/// Undirected multigraph on canonical *k*-mers (adjacency lists with parallel edge multiplicity).
#[derive(Debug, Clone)]
pub struct DbgGraph {
    pub k: usize,
    /// Undirected adjacency: `adj[u][v]` counts traversals u–v (stored for both endpoints).
    pub adj: BTreeMap<Vec<u8>, BTreeMap<Vec<u8>, u64>>,
    /// Multiplicity from the trusted *k*-mer table (vertex weight).
    pub node_mul: BTreeMap<Vec<u8>, u64>,
}

impl DbgGraph {
    pub fn new(k: usize, node_mul: BTreeMap<Vec<u8>, u64>) -> Self {
        Self {
            k,
            adj: BTreeMap::new(),
            node_mul,
        }
    }

    pub fn degree(&self, u: &[u8]) -> usize {
        self.adj.get(u).map(|m| m.len()).unwrap_or(0)
    }

    pub(crate) fn compact_view(&self) -> CompactDbgGraph {
        let mut name_to_id = BTreeMap::new();
        for node in self.node_mul.keys() {
            let next = NodeId(name_to_id.len());
            name_to_id.entry(node.clone()).or_insert(next);
        }
        for (node, neighbors) in &self.adj {
            let next = NodeId(name_to_id.len());
            name_to_id.entry(node.clone()).or_insert(next);
            for neighbor in neighbors.keys() {
                let next = NodeId(name_to_id.len());
                name_to_id.entry(neighbor.clone()).or_insert(next);
            }
        }

        let mut node_names = vec![Vec::new(); name_to_id.len()];
        let mut node_mul = vec![0u64; name_to_id.len()];
        for (node, id) in &name_to_id {
            node_names[id.0] = node.clone();
            node_mul[id.0] = self.node_mul.get(node).copied().unwrap_or(0);
        }

        let mut adj = vec![Vec::new(); name_to_id.len()];
        for (node, neighbors) in &self.adj {
            let Some(from) = name_to_id.get(node).copied() else {
                continue;
            };
            adj[from.0].reserve(neighbors.len());
            for (neighbor, weight) in neighbors {
                if let Some(to) = name_to_id.get(neighbor).copied() {
                    adj[from.0].push((to, *weight));
                }
            }
        }

        CompactDbgGraph {
            node_names,
            node_mul,
            adj,
        }
    }

    pub(crate) fn add_undirected_edge(
        &mut self,
        a: &[u8],
        b: &[u8],
        w: u64,
    ) -> Result<(), GraphError> {
        if a == b {
            return Err(GraphError::SelfLoop);
        }
        if w == 0 {
            return Ok(());
        }
        *self
            .adj
            .entry(a.to_vec())
            .or_default()
            .entry(b.to_vec())
            .or_insert(0) += w;
        *self
            .adj
            .entry(b.to_vec())
            .or_default()
            .entry(a.to_vec())
            .or_insert(0) += w;
        Ok(())
    }

    /// Restore a simplified **DBG** from checkpoint parts: vertices from `node_mul`, then undirected edges.
    pub fn from_checkpoint_parts(
        k: usize,
        node_mul: BTreeMap<Vec<u8>, u64>,
        edges: Vec<(Vec<u8>, Vec<u8>, u64)>,
    ) -> Result<Self, GraphError> {
        let mut g = DbgGraph::new(k, node_mul);
        for (u, v, w) in edges {
            g.add_undirected_edge(&u, &v, w)?;
        }
        Ok(g)
    }

    /// Remove an undirected edge entirely (both directions).
    pub fn remove_undirected_edge(&mut self, a: &[u8], b: &[u8]) {
        if let Some(m) = self.adj.get_mut(a) {
            m.remove(b);
            if m.is_empty() {
                self.adj.remove(a);
            }
        }
        if let Some(m) = self.adj.get_mut(b) {
            m.remove(a);
            if m.is_empty() {
                self.adj.remove(b);
            }
        }
    }

    /// Remove vertex and all incident edges from `adj`.
    pub fn remove_vertex_from_adj(&mut self, u: &[u8]) {
        if let Some(neigh) = self.adj.remove(u) {
            for v in neigh.keys() {
                if let Some(m) = self.adj.get_mut(v.as_slice()) {
                    m.remove(u);
                    if m.is_empty() {
                        self.adj.remove(v);
                    }
                }
            }
        }
    }
}

/// Build DBG from reads using **forward** consecutive windows for overlap correctness.
pub fn build_dbg(
    reads: &[Read],
    k: usize,
    trusted: &[(Vec<u8>, u64)],
) -> Result<DbgGraph, GraphError> {
    let mut mul: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
    for (key, c) in trusted {
        mul.insert(key.clone(), *c);
    }
    let mut g = DbgGraph::new(k, mul);
    let trusted_set: HashMap<Vec<u8>, u64> = trusted.iter().cloned().collect();

    for read in reads {
        for seg in n_free_acgt_segments(&read.sequence) {
            if seg.len() <= k {
                continue;
            }
            for i in 0..(seg.len() - k) {
                let f1 = &seg[i..i + k];
                let f2 = &seg[i + 1..i + k + 1];
                debug_assert_eq!(&f1[1..], &f2[..k - 1]);
                let u = canonical_kmer(f1);
                let v = canonical_kmer(f2);
                // Adjacent forward windows can **canonicalize to the same** *k*-mer (e.g. homopolymers);
                // an undirected self-edge is forbidden (**Phase-1 simplified graph invariants**), so skip.
                if u == v {
                    continue;
                }
                if trusted_set.contains_key(&u) && trusted_set.contains_key(&v) {
                    g.add_undirected_edge(&u, &v, 1)?;
                }
            }
        }
    }
    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::illumina::read::Read;

    #[test]
    fn homopolymer_consecutive_windows_skip_self_edge() {
        let reads = vec![Read {
            id: "h1".into(),
            sequence: b"AAAAAAAAAAAA".to_vec(),
        }];
        let trusted = vec![(b"AAAA".to_vec(), 10u64)];
        let g = build_dbg(&reads, 4, &trusted).expect("build_dbg");
        for (u, mp) in &g.adj {
            assert!(!mp.contains_key(u), "unexpected self-loop at {u:?}");
        }
    }

    #[test]
    fn compact_view_preserves_sorted_names_multiplicity_degree_and_weights() {
        let mut nodes = BTreeMap::new();
        nodes.insert(b"AAA".to_vec(), 3);
        nodes.insert(b"AAC".to_vec(), 5);
        nodes.insert(b"ACC".to_vec(), 7);
        let mut g = DbgGraph::new(3, nodes);
        g.add_undirected_edge(b"AAA", b"AAC", 2).expect("edge");
        g.add_undirected_edge(b"AAC", b"ACC", 4).expect("edge");

        let compact = g.compact_view();
        let ids: Vec<_> = compact.node_ids().collect();
        let a = ids[0];
        let b = ids[1];
        let c = ids[2];

        assert_eq!(ids.len(), 3);
        assert_eq!(compact.node_name(a), b"AAA");
        assert_eq!(compact.node_mul(b), 5);
        assert_eq!(compact.degree(b), 2);
        assert_eq!(compact.edge_weight(a, b), 2);
        assert_eq!(compact.edge_weight(b, a), 2);
        assert_eq!(compact.edge_weight(b, c), 4);
        assert_eq!(compact.edge_weight(a, c), 0);
    }
}
