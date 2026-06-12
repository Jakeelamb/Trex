//! Report-only contig endpoint diagnosis for fragmentation triage.

use serde::{Deserialize, Serialize};

use crate::dbg::{DbgGraph, GraphAnnotations};

pub const FRAGMENTATION_REPORT_SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FragmentStopReason {
    EmptyContigPath,
    GraphDeadEnd,
    BranchTangle,
    RepeatSuspected,
    LinearEndpoint,
    CircularOrInternal,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FragmentEndpointReport {
    pub contig: String,
    pub side: String,
    pub node: Option<String>,
    pub degree: usize,
    pub multiplicity: u64,
    pub repeat_suspected: bool,
    pub reason: FragmentStopReason,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct FragmentationSummary {
    pub contigs: usize,
    pub contig_path_nodes: usize,
    pub endpoints: usize,
    pub graph_dead_end_endpoints: usize,
    pub branch_tangle_endpoints: usize,
    pub repeat_suspected_endpoints: usize,
    pub linear_endpoint_endpoints: usize,
    pub circular_or_internal_endpoints: usize,
    pub empty_contig_paths: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FragmentationReport {
    pub schema_version: u64,
    pub summary: FragmentationSummary,
    pub endpoints: Vec<FragmentEndpointReport>,
}

pub fn diagnose_fragmentation(
    graph: &DbgGraph,
    contig_paths: &[Vec<Vec<u8>>],
    annotations: &GraphAnnotations,
) -> FragmentationReport {
    let mut summary = FragmentationSummary {
        contigs: contig_paths.len(),
        contig_path_nodes: contig_paths.iter().map(Vec::len).sum(),
        ..FragmentationSummary::default()
    };
    let mut endpoints = Vec::new();

    for (idx, path) in contig_paths.iter().enumerate() {
        if path.is_empty() {
            summary.empty_contig_paths += 1;
            continue;
        }
        let contig = format!("ctg{}", idx + 1);
        endpoints.push(classify_endpoint(
            graph,
            annotations,
            contig.clone(),
            "left",
            &path[0],
            path.len(),
        ));
        endpoints.push(classify_endpoint(
            graph,
            annotations,
            contig,
            "right",
            &path[path.len() - 1],
            path.len(),
        ));
    }

    summary.endpoints = endpoints.len();
    for endpoint in &endpoints {
        match endpoint.reason {
            FragmentStopReason::EmptyContigPath => summary.empty_contig_paths += 1,
            FragmentStopReason::GraphDeadEnd => summary.graph_dead_end_endpoints += 1,
            FragmentStopReason::BranchTangle => summary.branch_tangle_endpoints += 1,
            FragmentStopReason::RepeatSuspected => summary.repeat_suspected_endpoints += 1,
            FragmentStopReason::LinearEndpoint => summary.linear_endpoint_endpoints += 1,
            FragmentStopReason::CircularOrInternal => summary.circular_or_internal_endpoints += 1,
        }
    }

    FragmentationReport {
        schema_version: FRAGMENTATION_REPORT_SCHEMA_VERSION,
        summary,
        endpoints,
    }
}

fn classify_endpoint(
    graph: &DbgGraph,
    annotations: &GraphAnnotations,
    contig: String,
    side: &str,
    node: &[u8],
    path_len: usize,
) -> FragmentEndpointReport {
    let node_name = String::from_utf8_lossy(node).into_owned();
    let degree = graph.degree(node);
    let multiplicity = graph.node_mul.get(node).copied().unwrap_or(0);
    let repeat_suspected = annotations
        .nodes
        .get(&node_name)
        .map(|annotation| annotation.repeat_suspected)
        .unwrap_or(false);
    let reason = if degree <= 1 {
        FragmentStopReason::GraphDeadEnd
    } else if degree > 2 {
        FragmentStopReason::BranchTangle
    } else if repeat_suspected {
        FragmentStopReason::RepeatSuspected
    } else if path_len > 1 && degree == 2 {
        FragmentStopReason::LinearEndpoint
    } else {
        FragmentStopReason::CircularOrInternal
    };

    FragmentEndpointReport {
        contig,
        side: side.to_string(),
        node: Some(node_name),
        degree,
        multiplicity,
        repeat_suspected,
        reason,
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{diagnose_fragmentation, FragmentStopReason};
    use crate::dbg::{annotate_graph, DbgGraph};

    #[test]
    fn classifies_dead_end_branch_repeat_and_linear_endpoints() {
        let graph = DbgGraph::from_checkpoint_parts(
            3,
            BTreeMap::from([
                (b"AAA".to_vec(), 10),
                (b"AAC".to_vec(), 10),
                (b"ACC".to_vec(), 10),
                (b"CCC".to_vec(), 10),
                (b"CCG".to_vec(), 10),
                (b"CGG".to_vec(), 30),
                (b"GGT".to_vec(), 10),
            ]),
            vec![
                (b"AAA".to_vec(), b"AAC".to_vec(), 1),
                (b"AAC".to_vec(), b"ACC".to_vec(), 1),
                (b"ACC".to_vec(), b"CCC".to_vec(), 1),
                (b"ACC".to_vec(), b"CCG".to_vec(), 1),
                (b"ACC".to_vec(), b"CGG".to_vec(), 1),
                (b"CGG".to_vec(), b"GGT".to_vec(), 1),
            ],
        )
        .expect("graph");
        let contig_paths = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
            vec![b"CGG".to_vec(), b"GGT".to_vec()],
        ];
        let annotations = annotate_graph(&graph, &contig_paths);

        let report = diagnose_fragmentation(&graph, &contig_paths, &annotations);

        assert_eq!(report.summary.contigs, 3);
        assert_eq!(report.summary.endpoints, 6);
        assert_eq!(report.summary.graph_dead_end_endpoints, 3);
        assert_eq!(report.summary.linear_endpoint_endpoints, 1);
        assert_eq!(report.summary.branch_tangle_endpoints, 1);
        assert_eq!(report.summary.repeat_suspected_endpoints, 1);
        assert!(report
            .endpoints
            .iter()
            .any(|endpoint| endpoint.node.as_deref() == Some("CGG")
                && endpoint.reason == FragmentStopReason::RepeatSuspected));
    }

    #[test]
    fn counts_empty_contig_paths_without_endpoint_rows() {
        let graph = DbgGraph::new(3, BTreeMap::new());
        let annotations = annotate_graph(&graph, &[Vec::new()]);

        let report = diagnose_fragmentation(&graph, &[Vec::new()], &annotations);

        assert_eq!(report.summary.contigs, 1);
        assert_eq!(report.summary.empty_contig_paths, 1);
        assert_eq!(report.summary.endpoints, 0);
        assert!(report.endpoints.is_empty());
    }
}
