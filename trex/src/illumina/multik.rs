//! Explicit multi-k candidate scoring and one-graph selection.

use serde::{Deserialize, Serialize};

use crate::dbg::{annotate_graph, build_dbg, extract_unitigs};
use crate::error::{IngestError, TrexError};
use crate::illumina::counts::enumerate_sorted_counts;
use crate::illumina::read::Read;
use crate::kmer::apply_trusted_threshold;

pub const MULTI_K_SELECTION_SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone, Default)]
pub struct MultiKParams {
    pub ladder: Vec<usize>,
}

impl MultiKParams {
    pub fn enabled(&self) -> bool {
        !self.ladder.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MultiKCandidateReport {
    pub k: usize,
    pub feasible: bool,
    pub reason: Option<String>,
    pub unique_kmers: usize,
    pub trusted_kmers: usize,
    pub graph_nodes: usize,
    pub graph_edges: usize,
    pub unitigs: usize,
    pub branch_nodes: usize,
    pub repeat_suspected_nodes: usize,
    pub read_kmer_completeness: f64,
    pub contiguity_score: f64,
    pub branchiness_score: f64,
    pub repeat_risk_score: f64,
    pub total_score: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MultiKSelectionReport {
    pub schema_version: u64,
    pub enabled: bool,
    pub requested_k: usize,
    pub selected_k: usize,
    pub candidates: Vec<MultiKCandidateReport>,
}

impl MultiKSelectionReport {
    pub fn single_k(k: usize) -> Self {
        Self {
            schema_version: MULTI_K_SELECTION_SCHEMA_VERSION,
            enabled: false,
            requested_k: k,
            selected_k: k,
            candidates: Vec::new(),
        }
    }
}

pub fn select_k(
    reads: &[Read],
    requested_k: usize,
    trusted_threshold: u64,
    multi_k: &MultiKParams,
) -> Result<MultiKSelectionReport, TrexError> {
    if !multi_k.enabled() {
        return Ok(MultiKSelectionReport::single_k(requested_k));
    }

    let shortest = Read::shortest_length(reads).ok_or_else(|| {
        TrexError::Ingest(IngestError::FastqFormat("no reads ingested".to_string()))
    })?;
    let mut candidates = Vec::new();
    let mut best: Option<MultiKCandidateReport> = None;
    for k in normalized_ladder(&multi_k.ladder) {
        let report = score_candidate(reads, k, shortest, trusted_threshold)?;
        if report.feasible {
            let replace = match &best {
                None => true,
                Some(prev) => better_candidate(&report, prev),
            };
            if replace {
                best = Some(report.clone());
            }
        }
        candidates.push(report);
    }

    let Some(best) = best else {
        return Err(IngestError::KTooLarge {
            k: *multi_k.ladder.iter().min().unwrap_or(&requested_k),
            shortest,
        }
        .into());
    };

    Ok(MultiKSelectionReport {
        schema_version: MULTI_K_SELECTION_SCHEMA_VERSION,
        enabled: true,
        requested_k,
        selected_k: best.k,
        candidates,
    })
}

fn normalized_ladder(ladder: &[usize]) -> Vec<usize> {
    let mut out: Vec<usize> = ladder.iter().copied().filter(|k| *k > 0).collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn score_candidate(
    reads: &[Read],
    k: usize,
    shortest: usize,
    trusted_threshold: u64,
) -> Result<MultiKCandidateReport, TrexError> {
    if k > shortest {
        return Ok(infeasible_candidate(
            k,
            format!("k ({k}) exceeds shortest post-preprocess read length ({shortest})"),
        ));
    }
    let merged = enumerate_sorted_counts(reads, k)?;
    let unique_kmers = merged.len();
    let trusted = apply_trusted_threshold(merged, trusted_threshold);
    let graph = build_dbg(reads, k, &trusted)?;
    let unitig_paths = extract_unitigs(&graph);
    let annotations = annotate_graph(&graph, &unitig_paths);
    let graph_edges = graph
        .adj
        .values()
        .map(|neighbors| neighbors.len())
        .sum::<usize>()
        / 2;
    let branch_nodes = graph
        .adj
        .keys()
        .filter(|node| graph.degree(node) > 2)
        .count();
    let max_unitig_nodes = unitig_paths.iter().map(Vec::len).max().unwrap_or(0);
    let read_kmer_completeness = fraction(trusted.len(), unique_kmers);
    let contiguity_score = max_unitig_nodes as f64;
    let branchiness_score = fraction(branch_nodes, graph.node_mul.len());
    let repeat_risk_score = fraction(
        annotations.summary.repeat_suspected_nodes,
        annotations.summary.node_count,
    );
    let total_score = read_kmer_completeness.mul_add(
        1000.0,
        contiguity_score - branchiness_score * 100.0 - repeat_risk_score * 50.0,
    );

    Ok(MultiKCandidateReport {
        k,
        feasible: true,
        reason: None,
        unique_kmers,
        trusted_kmers: trusted.len(),
        graph_nodes: graph.node_mul.len(),
        graph_edges,
        unitigs: unitig_paths.len(),
        branch_nodes,
        repeat_suspected_nodes: annotations.summary.repeat_suspected_nodes,
        read_kmer_completeness,
        contiguity_score,
        branchiness_score,
        repeat_risk_score,
        total_score,
    })
}

fn infeasible_candidate(k: usize, reason: String) -> MultiKCandidateReport {
    MultiKCandidateReport {
        k,
        feasible: false,
        reason: Some(reason),
        unique_kmers: 0,
        trusted_kmers: 0,
        graph_nodes: 0,
        graph_edges: 0,
        unitigs: 0,
        branch_nodes: 0,
        repeat_suspected_nodes: 0,
        read_kmer_completeness: 0.0,
        contiguity_score: 0.0,
        branchiness_score: 1.0,
        repeat_risk_score: 1.0,
        total_score: f64::NEG_INFINITY,
    }
}

fn better_candidate(candidate: &MultiKCandidateReport, previous: &MultiKCandidateReport) -> bool {
    candidate
        .total_score
        .partial_cmp(&previous.total_score)
        .map(|ord| {
            ord.is_gt()
                || (ord.is_eq()
                    && (candidate.k > previous.k
                        || (candidate.k == previous.k
                            && candidate.graph_edges < previous.graph_edges)))
        })
        .unwrap_or(false)
}

fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[cfg(test)]
mod tests {
    use super::{select_k, MultiKParams};
    use crate::illumina::read::Read;

    #[test]
    fn empty_ladder_keeps_requested_single_k() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        let report = select_k(&reads, 4, 1, &MultiKParams::default()).expect("select");

        assert!(!report.enabled);
        assert_eq!(report.requested_k, 4);
        assert_eq!(report.selected_k, 4);
        assert!(report.candidates.is_empty());
    }

    #[test]
    fn ladder_scores_feasible_candidates_and_selects_one() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
        }];
        let report = select_k(
            &reads,
            4,
            1,
            &MultiKParams {
                ladder: vec![3, 4, 5],
            },
        )
        .expect("select");

        assert!(report.enabled);
        assert!([3, 4, 5].contains(&report.selected_k));
        assert_eq!(report.candidates.len(), 3);
        assert!(report.candidates.iter().all(|candidate| candidate.feasible));
    }

    #[test]
    fn ladder_records_infeasible_candidates() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGT".to_vec(),
        }];
        let report = select_k(&reads, 4, 1, &MultiKParams { ladder: vec![3, 9] }).expect("select");

        assert_eq!(report.selected_k, 3);
        assert_eq!(report.candidates.len(), 2);
        assert!(!report.candidates[1].feasible);
        assert!(report.candidates[1]
            .reason
            .as_deref()
            .unwrap_or_default()
            .contains("exceeds shortest"));
    }
}
