//! Copy-number and repeat-risk annotations over a simplified DBG without mutating topology.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::dbg::graph::DbgGraph;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NodeDepthClass {
    LowCopy,
    SingleCopy,
    HighCopy,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NodeAnnotation {
    pub multiplicity: u64,
    pub degree: usize,
    pub depth_class: NodeDepthClass,
    pub repeat_suspected: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnitigAnnotation {
    pub unitig_index: usize,
    pub node_count: usize,
    pub min_multiplicity: u64,
    pub max_multiplicity: u64,
    pub mean_multiplicity: f64,
    pub repeat_suspected: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GraphAnnotationSummary {
    pub baseline_multiplicity: u64,
    pub node_count: usize,
    pub low_copy_nodes: usize,
    pub single_copy_nodes: usize,
    pub high_copy_nodes: usize,
    pub repeat_suspected_nodes: usize,
    pub unitig_count: usize,
    pub repeat_suspected_unitigs: usize,
    pub unitigs: Vec<UnitigAnnotation>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GraphAnnotations {
    pub nodes: BTreeMap<String, NodeAnnotation>,
    pub summary: GraphAnnotationSummary,
}

pub fn annotate_graph(graph: &DbgGraph, unitig_paths: &[Vec<Vec<u8>>]) -> GraphAnnotations {
    let active_nodes = active_graph_nodes(graph, unitig_paths);
    let baseline = median_multiplicity(graph, &active_nodes);
    let mut nodes = BTreeMap::new();
    let mut low_copy_nodes = 0usize;
    let mut single_copy_nodes = 0usize;
    let mut high_copy_nodes = 0usize;
    let mut repeat_suspected_nodes = 0usize;

    for node in &active_nodes {
        let multiplicity = graph.node_mul.get(node).copied().unwrap_or(0);
        let degree = graph.degree(node);
        let depth_class = classify_depth(multiplicity, baseline);
        match depth_class {
            NodeDepthClass::LowCopy => low_copy_nodes += 1,
            NodeDepthClass::SingleCopy => single_copy_nodes += 1,
            NodeDepthClass::HighCopy => high_copy_nodes += 1,
        }
        let repeat_suspected = depth_class == NodeDepthClass::HighCopy || degree > 2;
        if repeat_suspected {
            repeat_suspected_nodes += 1;
        }
        nodes.insert(
            String::from_utf8_lossy(node).into_owned(),
            NodeAnnotation {
                multiplicity,
                degree,
                depth_class,
                repeat_suspected,
            },
        );
    }

    let unitigs: Vec<UnitigAnnotation> = unitig_paths
        .iter()
        .enumerate()
        .map(|(idx, path)| annotate_unitig(graph, idx, path, &nodes))
        .collect();
    let repeat_suspected_unitigs = unitigs
        .iter()
        .filter(|unitig| unitig.repeat_suspected)
        .count();

    GraphAnnotations {
        nodes,
        summary: GraphAnnotationSummary {
            baseline_multiplicity: baseline,
            node_count: active_nodes.len(),
            low_copy_nodes,
            single_copy_nodes,
            high_copy_nodes,
            repeat_suspected_nodes,
            unitig_count: unitigs.len(),
            repeat_suspected_unitigs,
            unitigs,
        },
    }
}

fn active_graph_nodes(graph: &DbgGraph, unitig_paths: &[Vec<Vec<u8>>]) -> BTreeSet<Vec<u8>> {
    let mut nodes: BTreeSet<Vec<u8>> = graph.adj.keys().cloned().collect();
    for path in unitig_paths {
        for node in path {
            nodes.insert(node.clone());
        }
    }
    nodes
}

fn median_multiplicity(graph: &DbgGraph, active_nodes: &BTreeSet<Vec<u8>>) -> u64 {
    let mut values: Vec<u64> = active_nodes
        .iter()
        .filter_map(|node| graph.node_mul.get(node).copied())
        .filter(|value| *value > 0)
        .collect();
    if values.is_empty() {
        return 0;
    }
    values.sort_unstable();
    values[values.len() / 2]
}

fn classify_depth(multiplicity: u64, baseline: u64) -> NodeDepthClass {
    if baseline == 0 {
        return NodeDepthClass::SingleCopy;
    }
    if multiplicity.saturating_mul(2) < baseline {
        NodeDepthClass::LowCopy
    } else if multiplicity >= baseline.saturating_mul(2) && multiplicity > baseline {
        NodeDepthClass::HighCopy
    } else {
        NodeDepthClass::SingleCopy
    }
}

fn annotate_unitig(
    graph: &DbgGraph,
    unitig_index: usize,
    path: &[Vec<u8>],
    nodes: &BTreeMap<String, NodeAnnotation>,
) -> UnitigAnnotation {
    let mut min_multiplicity = u64::MAX;
    let mut max_multiplicity = 0u64;
    let mut total = 0u64;
    let mut repeat_suspected = false;
    let mut observed = 0usize;

    for node in path {
        let multiplicity = graph.node_mul.get(node).copied().unwrap_or(0);
        min_multiplicity = min_multiplicity.min(multiplicity);
        max_multiplicity = max_multiplicity.max(multiplicity);
        total = total.saturating_add(multiplicity);
        observed += 1;
        if nodes
            .get(&String::from_utf8_lossy(node).into_owned())
            .map(|annotation| annotation.repeat_suspected)
            .unwrap_or(false)
        {
            repeat_suspected = true;
        }
    }

    if observed == 0 {
        min_multiplicity = 0;
    }

    UnitigAnnotation {
        unitig_index,
        node_count: observed,
        min_multiplicity,
        max_multiplicity,
        mean_multiplicity: if observed == 0 {
            0.0
        } else {
            total as f64 / observed as f64
        },
        repeat_suspected,
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{annotate_graph, NodeDepthClass};
    use crate::dbg::graph::DbgGraph;

    fn graph_with_multiplicities(rows: &[(&[u8], u64)]) -> DbgGraph {
        let mut node_mul = BTreeMap::new();
        for (node, multiplicity) in rows {
            node_mul.insert(node.to_vec(), *multiplicity);
        }
        DbgGraph::new(3, node_mul)
    }

    #[test]
    fn classifies_low_single_and_high_copy_nodes() {
        let graph = graph_with_multiplicities(&[(b"AAA", 4), (b"AAC", 10), (b"ACC", 25)]);
        let annotations = annotate_graph(
            &graph,
            &[vec![b"AAA".to_vec(), b"AAC".to_vec(), b"ACC".to_vec()]],
        );

        assert_eq!(annotations.summary.baseline_multiplicity, 10);
        assert_eq!(annotations.summary.low_copy_nodes, 1);
        assert_eq!(annotations.summary.single_copy_nodes, 1);
        assert_eq!(annotations.summary.high_copy_nodes, 1);
        assert_eq!(
            annotations.nodes["AAA"].depth_class,
            NodeDepthClass::LowCopy
        );
        assert_eq!(
            annotations.nodes["AAC"].depth_class,
            NodeDepthClass::SingleCopy
        );
        assert_eq!(
            annotations.nodes["ACC"].depth_class,
            NodeDepthClass::HighCopy
        );
    }

    #[test]
    fn flags_branching_and_high_copy_repeat_suspicion() {
        let graph = DbgGraph::from_checkpoint_parts(
            3,
            BTreeMap::from([
                (b"AAA".to_vec(), 10),
                (b"AAC".to_vec(), 10),
                (b"ACC".to_vec(), 25),
                (b"CCC".to_vec(), 10),
                (b"CCG".to_vec(), 10),
            ]),
            vec![
                (b"AAA".to_vec(), b"AAC".to_vec(), 1),
                (b"AAA".to_vec(), b"ACC".to_vec(), 1),
                (b"AAA".to_vec(), b"CCC".to_vec(), 1),
                (b"CCC".to_vec(), b"CCG".to_vec(), 1),
            ],
        )
        .expect("graph");

        let annotations = annotate_graph(&graph, &[vec![b"AAA".to_vec(), b"AAC".to_vec()]]);

        assert!(annotations.nodes["AAA"].repeat_suspected);
        assert!(annotations.nodes["ACC"].repeat_suspected);
        assert_eq!(annotations.summary.repeat_suspected_nodes, 2);
        assert_eq!(annotations.summary.repeat_suspected_unitigs, 1);
    }

    #[test]
    fn computes_unitig_min_mean_max_multiplicity() {
        let graph = graph_with_multiplicities(&[(b"AAA", 5), (b"AAC", 10), (b"ACC", 15)]);
        let annotations = annotate_graph(
            &graph,
            &[vec![b"AAA".to_vec(), b"AAC".to_vec(), b"ACC".to_vec()]],
        );
        let unitig = &annotations.summary.unitigs[0];

        assert_eq!(unitig.unitig_index, 0);
        assert_eq!(unitig.node_count, 3);
        assert_eq!(unitig.min_multiplicity, 5);
        assert_eq!(unitig.max_multiplicity, 15);
        assert_eq!(unitig.mean_multiplicity, 10.0);
    }

    #[test]
    fn ignores_stale_node_multiplicity_not_present_in_active_graph() {
        let graph = graph_with_multiplicities(&[(b"AAA", 10), (b"AAC", 11), (b"ZZZ", 1000)]);
        let annotations = annotate_graph(&graph, &[vec![b"AAA".to_vec(), b"AAC".to_vec()]]);

        assert_eq!(annotations.summary.node_count, 2);
        assert_eq!(annotations.summary.baseline_multiplicity, 11);
        assert!(!annotations.nodes.contains_key("ZZZ"));
    }
}
