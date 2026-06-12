//! Explicit multi-k candidate scoring and one-graph selection.

use serde::{Deserialize, Serialize};

use crate::dbg::{annotate_graph, build_dbg, extract_unitigs};
use crate::error::{IngestError, TrexError};
use crate::illumina::counts::enumerate_sorted_counts;
use crate::illumina::read::Read;
use crate::kmer::apply_trusted_threshold;

pub const MULTI_K_SELECTION_SCHEMA_VERSION: u64 = 2;

#[derive(Debug, Clone, Default)]
pub struct MultiKParams {
    pub auto: bool,
    pub ladder: Vec<usize>,
}

impl MultiKParams {
    pub fn enabled(&self) -> bool {
        self.auto || self.ladder.iter().any(|k| *k > 0)
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
    pub graph_density: f64,
    pub unitigs: usize,
    pub max_unitig_nodes: usize,
    pub unitig_n50_nodes: usize,
    pub dead_end_nodes: usize,
    pub branch_nodes: usize,
    pub tangle_nodes: usize,
    pub repeat_suspected_nodes: usize,
    pub read_kmer_completeness: f64,
    pub contiguity_score: f64,
    pub branchiness_score: f64,
    pub dead_end_score: f64,
    pub tangle_score: f64,
    pub repeat_risk_score: f64,
    pub graph_density_score: f64,
    pub score_terms: MultiKScoreTerms,
    pub total_score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MultiKScoreTerms {
    pub read_kmer_completeness: f64,
    pub contiguity: f64,
    pub branchiness_penalty: f64,
    pub dead_end_penalty: f64,
    pub tangle_penalty: f64,
    pub repeat_risk_penalty: f64,
    pub graph_density_penalty: f64,
}

impl MultiKScoreTerms {
    fn total(self) -> f64 {
        self.read_kmer_completeness + self.contiguity
            - self.branchiness_penalty
            - self.dead_end_penalty
            - self.tangle_penalty
            - self.repeat_risk_penalty
            - self.graph_density_penalty
    }
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
    let ladder = resolved_ladder(multi_k, shortest);
    for k in ladder.iter().copied() {
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
            k: *ladder.iter().min().unwrap_or(&requested_k),
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

fn resolved_ladder(multi_k: &MultiKParams, shortest: usize) -> Vec<usize> {
    let mut out = normalized_ladder(&multi_k.ladder);
    if multi_k.auto {
        out.extend(auto_ladder_for_shortest_read(shortest));
        out.sort_unstable();
        out.dedup();
    }
    out
}

pub fn auto_ladder_for_shortest_read(shortest: usize) -> Vec<usize> {
    let max_k = largest_odd_at_most(shortest.min(127));
    if max_k < 3 {
        return Vec::new();
    }
    let candidates: &[usize] = if max_k <= 37 {
        &[25, 29, 31, 33, 35]
    } else if max_k <= 75 {
        &[21, 33, 45, 55, 65, 75]
    } else {
        &[21, 33, 55, 77, 99, 127]
    };
    let mut out: Vec<usize> = candidates.iter().copied().filter(|k| *k <= max_k).collect();
    if out.is_empty() || (!out.contains(&max_k) && max_k >= 21) {
        out.push(max_k);
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn largest_odd_at_most(value: usize) -> usize {
    if value % 2 == 1 {
        value
    } else {
        value.saturating_sub(1)
    }
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
    let graph_nodes = graph.node_mul.len();
    let graph_density = fraction(graph_edges, graph_nodes);
    let dead_end_nodes = graph
        .node_mul
        .keys()
        .filter(|node| graph.degree(node) <= 1)
        .count();
    let branch_nodes = graph
        .node_mul
        .keys()
        .filter(|node| graph.degree(node) > 2)
        .count();
    let tangle_nodes = graph
        .node_mul
        .keys()
        .filter(|node| graph.degree(node) > 3)
        .count();
    let max_unitig_nodes = unitig_paths.iter().map(Vec::len).max().unwrap_or(0);
    let unitig_n50_nodes = n50_unitig_nodes(&unitig_paths);
    let read_kmer_completeness = fraction(trusted.len(), unique_kmers);
    let contiguity_score = max_unitig_nodes as f64;
    let branchiness_score = fraction(branch_nodes, graph_nodes);
    let dead_end_score = fraction(dead_end_nodes, graph_nodes);
    let tangle_score = fraction(tangle_nodes, graph_nodes);
    let repeat_risk_score = fraction(
        annotations.summary.repeat_suspected_nodes,
        annotations.summary.node_count,
    );
    let graph_density_score = (graph_density - 1.25).max(0.0);
    let score_terms = MultiKScoreTerms {
        read_kmer_completeness: read_kmer_completeness * 1000.0,
        contiguity: contiguity_score,
        branchiness_penalty: branchiness_score * 100.0,
        dead_end_penalty: dead_end_score * 50.0,
        tangle_penalty: tangle_score * 150.0,
        repeat_risk_penalty: repeat_risk_score * 50.0,
        graph_density_penalty: graph_density_score * 25.0,
    };
    let total_score = score_terms.total();

    Ok(MultiKCandidateReport {
        k,
        feasible: true,
        reason: None,
        unique_kmers,
        trusted_kmers: trusted.len(),
        graph_nodes,
        graph_edges,
        graph_density,
        unitigs: unitig_paths.len(),
        max_unitig_nodes,
        unitig_n50_nodes,
        dead_end_nodes,
        branch_nodes,
        tangle_nodes,
        repeat_suspected_nodes: annotations.summary.repeat_suspected_nodes,
        read_kmer_completeness,
        contiguity_score,
        branchiness_score,
        dead_end_score,
        tangle_score,
        repeat_risk_score,
        graph_density_score,
        score_terms,
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
        graph_density: 0.0,
        unitigs: 0,
        max_unitig_nodes: 0,
        unitig_n50_nodes: 0,
        dead_end_nodes: 0,
        branch_nodes: 0,
        tangle_nodes: 0,
        repeat_suspected_nodes: 0,
        read_kmer_completeness: 0.0,
        contiguity_score: 0.0,
        branchiness_score: 1.0,
        dead_end_score: 1.0,
        tangle_score: 1.0,
        repeat_risk_score: 1.0,
        graph_density_score: 1.0,
        score_terms: MultiKScoreTerms {
            read_kmer_completeness: 0.0,
            contiguity: 0.0,
            branchiness_penalty: 100.0,
            dead_end_penalty: 50.0,
            tangle_penalty: 150.0,
            repeat_risk_penalty: 50.0,
            graph_density_penalty: 25.0,
        },
        total_score: f64::NEG_INFINITY,
    }
}

fn n50_unitig_nodes(unitig_paths: &[Vec<Vec<u8>>]) -> usize {
    let mut lengths: Vec<usize> = unitig_paths
        .iter()
        .map(Vec::len)
        .filter(|len| *len > 0)
        .collect();
    if lengths.is_empty() {
        return 0;
    }
    let total: usize = lengths.iter().sum();
    let midpoint = total.div_ceil(2);
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let mut running = 0usize;
    for len in lengths {
        running += len;
        if running >= midpoint {
            return len;
        }
    }
    0
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
    use super::{auto_ladder_for_shortest_read, select_k, MultiKParams};
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
    fn zero_only_ladder_keeps_requested_single_k() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        let report = select_k(
            &reads,
            4,
            1,
            &MultiKParams {
                auto: false,
                ladder: vec![0],
            },
        )
        .expect("select");

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
                auto: false,
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
    fn candidate_report_includes_spades_style_graph_pressure_terms() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
        }];
        let report = select_k(
            &reads,
            4,
            1,
            &MultiKParams {
                auto: false,
                ladder: vec![4],
            },
        )
        .expect("select");

        let candidate = &report.candidates[0];
        assert!(candidate.feasible);
        assert!(candidate.graph_density > 0.0);
        assert_eq!(candidate.max_unitig_nodes, candidate.unitig_n50_nodes);
        assert!(candidate.dead_end_nodes > 0);
        assert_eq!(
            candidate.dead_end_score,
            candidate.dead_end_nodes as f64 / candidate.graph_nodes as f64
        );
        assert_eq!(
            candidate.total_score,
            candidate.score_terms.read_kmer_completeness + candidate.score_terms.contiguity
                - candidate.score_terms.branchiness_penalty
                - candidate.score_terms.dead_end_penalty
                - candidate.score_terms.tangle_penalty
                - candidate.score_terms.repeat_risk_penalty
                - candidate.score_terms.graph_density_penalty
        );
    }

    #[test]
    fn ladder_records_infeasible_candidates() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGT".to_vec(),
        }];
        let report = select_k(
            &reads,
            4,
            1,
            &MultiKParams {
                auto: false,
                ladder: vec![3, 9],
            },
        )
        .expect("select");

        assert_eq!(report.selected_k, 3);
        assert_eq!(report.candidates.len(), 2);
        assert!(!report.candidates[1].feasible);
        assert!(report.candidates[1]
            .reason
            .as_deref()
            .unwrap_or_default()
            .contains("exceeds shortest"));
        assert!(report.candidates[1].score_terms.dead_end_penalty > 0.0);
        assert!(report.candidates[1].score_terms.tangle_penalty > 0.0);
    }

    #[test]
    fn auto_ladder_for_short_short_reads_includes_short_read_candidates() {
        assert_eq!(auto_ladder_for_shortest_read(36), vec![25, 29, 31, 33, 35]);
    }

    #[test]
    fn auto_ladder_caps_long_reads_at_127() {
        assert_eq!(
            auto_ladder_for_shortest_read(151),
            vec![21, 33, 55, 77, 99, 127]
        );
    }

    #[test]
    fn auto_ladder_handles_very_short_reads() {
        assert_eq!(auto_ladder_for_shortest_read(12), vec![11]);
    }

    #[test]
    fn auto_k_scores_derived_candidates() {
        let reads = vec![Read {
            id: "r1".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
        }];
        let report = select_k(
            &reads,
            21,
            1,
            &MultiKParams {
                auto: true,
                ladder: Vec::new(),
            },
        )
        .expect("select");

        assert!(report.enabled);
        assert_eq!(report.requested_k, 21);
        assert_eq!(report.candidates.len(), 1);
        assert_eq!(report.candidates[0].k, 11);
        assert_eq!(report.selected_k, 11);
    }
}
