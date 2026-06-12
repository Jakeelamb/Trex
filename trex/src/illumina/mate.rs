//! **Phase-2 Illumina mate usage**: strengthen existing **DBG** edges using paired-end layout
//! (does not alter **Phase-1 k-mer count representation** or trusted vertex weights).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::dbg::graph::DbgGraph;
use crate::evidence::{EvidenceKind, EvidenceRecord, EvidenceSourceStage, SupportCounts};
use crate::illumina::preprocess::n_free_acgt_segments;
use crate::illumina::read::Read;
use crate::kmer::canonical_kmer;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatePairOrientation {
    R1TailToR2Head,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct MateDistanceEvidence {
    pub insert_mean_bp: u64,
    pub insert_stddev_bp: Option<u64>,
    pub read_bases_bp: u64,
    pub estimated_gap_bp: i64,
    pub confidence: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MateBridgeCandidate {
    pub from_node: String,
    pub to_node: String,
    pub orientation: MatePairOrientation,
    pub distance: Option<MateDistanceEvidence>,
    pub support_pairs: usize,
    pub conflict_pairs: usize,
    pub score: u64,
    pub existing_dbg_edge: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MateBridgeStats {
    pub pairs_seen: usize,
    pub pairs_with_endpoint_kmers: usize,
    pub trusted_endpoint_pairs: usize,
    pub existing_edge_pairs: usize,
    pub report_only_pairs: usize,
    pub boosted_edges: usize,
    pub candidates: Vec<MateBridgeCandidate>,
}

type MateBridgeCandidateKey = (
    Vec<u8>,
    Vec<u8>,
    MatePairOrientation,
    Option<MateDistanceEvidence>,
    bool,
);

impl MateBridgeStats {
    pub fn evidence_record(&self) -> EvidenceRecord {
        EvidenceRecord::new(
            EvidenceKind::MateBridgeExistingEdge,
            EvidenceSourceStage::Phase2MateBridge,
            SupportCounts {
                observed: self.pairs_seen as u64,
                eligible: self.trusted_endpoint_pairs as u64,
                supporting: self.existing_edge_pairs as u64,
                applied: self.boosted_edges as u64,
            },
        )
        .with_counter("pairs_seen", self.pairs_seen as u64)
        .with_counter(
            "pairs_with_endpoint_kmers",
            self.pairs_with_endpoint_kmers as u64,
        )
        .with_counter("trusted_endpoint_pairs", self.trusted_endpoint_pairs as u64)
        .with_counter("existing_edge_pairs", self.existing_edge_pairs as u64)
        .with_counter("report_only_pairs", self.report_only_pairs as u64)
        .with_counter("boosted_edges", self.boosted_edges as u64)
    }
}

fn node_label(node: &[u8]) -> String {
    String::from_utf8_lossy(node).into_owned()
}

fn first_canonical_kmer_forward(seq: &[u8], k: usize) -> Option<Vec<u8>> {
    for seg in n_free_acgt_segments(seq) {
        if seg.len() >= k {
            return Some(canonical_kmer(&seg[..k]));
        }
    }
    None
}

fn last_canonical_kmer_forward(seq: &[u8], k: usize) -> Option<Vec<u8>> {
    let mut last: Option<Vec<u8>> = None;
    for seg in n_free_acgt_segments(seq) {
        if seg.len() >= k {
            let w = &seg[seg.len() - k..];
            last = Some(canonical_kmer(w));
        }
    }
    last
}

fn distance_evidence(
    insert_mean_bp: Option<u64>,
    insert_stddev_bp: Option<u64>,
    r1_len: usize,
    r2_len: usize,
) -> Option<MateDistanceEvidence> {
    let insert_mean_bp = insert_mean_bp?;
    let read_bases_bp = r1_len.saturating_add(r2_len) as u64;
    let estimated_gap_bp = insert_mean_bp as i64 - read_bases_bp as i64;
    let confidence = if insert_mean_bp == 0 {
        0
    } else if let Some(stddev) = insert_stddev_bp {
        let penalty = stddev.saturating_mul(100) / insert_mean_bp;
        100u64.saturating_sub(penalty.min(100))
    } else {
        50
    };
    Some(MateDistanceEvidence {
        insert_mean_bp,
        insert_stddev_bp,
        read_bases_bp,
        estimated_gap_bp,
        confidence,
    })
}

fn score_candidate(
    support_pairs: usize,
    conflict_pairs: usize,
    distance: Option<&MateDistanceEvidence>,
) -> u64 {
    let distance_confidence = distance.map(|distance| distance.confidence).unwrap_or(0);
    (support_pairs as u64)
        .saturating_mul(100)
        .saturating_add(distance_confidence)
        .saturating_sub((conflict_pairs as u64).saturating_mul(100))
}

/// For each mate pair, if **R1**'s last forward *k*-mer and **R2**'s first forward *k*-mer are both
/// trusted vertices and already adjacent in the **DBG**, increment that undirected edge weight by **1**.
///
/// This is a conservative **Phase-2 Illumina** bridge: no new edges, no changes to `node_mul`.
/// It runs only when the operator supplied an **insert-size prior** (see `assemble_illumina` caller).
pub fn boost_mate_pairs_on_existing_dbg_edges(
    graph: &mut DbgGraph,
    reads: &[Read],
    r1_count: usize,
    k: usize,
    insert_mean_bp: Option<u64>,
    insert_stddev_bp: Option<u64>,
) -> MateBridgeStats {
    if r1_count == 0 || reads.len() < r1_count * 2 {
        return MateBridgeStats::default();
    }
    let mut stats = MateBridgeStats::default();
    let mut candidates: BTreeMap<MateBridgeCandidateKey, usize> = BTreeMap::new();
    for i in 0..r1_count {
        stats.pairs_seen += 1;
        let r1 = &reads[i];
        let r2 = &reads[i + r1_count];
        let Some(a) = last_canonical_kmer_forward(&r1.sequence, k) else {
            continue;
        };
        let Some(b) = first_canonical_kmer_forward(&r2.sequence, k) else {
            continue;
        };
        stats.pairs_with_endpoint_kmers += 1;
        if !graph.node_mul.contains_key(&a) || !graph.node_mul.contains_key(&b) {
            continue;
        }
        stats.trusted_endpoint_pairs += 1;
        if a == b {
            continue;
        }
        let w0 = graph
            .adj
            .get(&a)
            .and_then(|m| m.get(&b))
            .copied()
            .unwrap_or(0);
        let existing_dbg_edge = w0 > 0;
        let orientation = MatePairOrientation::R1TailToR2Head;
        let distance = distance_evidence(
            insert_mean_bp,
            insert_stddev_bp,
            r1.sequence.len(),
            r2.sequence.len(),
        );
        *candidates
            .entry((
                a.clone(),
                b.clone(),
                orientation,
                distance,
                existing_dbg_edge,
            ))
            .or_insert(0) += 1;
        if existing_dbg_edge {
            stats.existing_edge_pairs += 1;
            if graph.add_undirected_edge(&a, &b, 1).is_ok() {
                stats.boosted_edges += 1;
            }
        } else {
            stats.report_only_pairs += 1;
        }
    }
    let mut from_support: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut to_support: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    let mut pair_support: BTreeMap<(Vec<u8>, Vec<u8>), usize> = BTreeMap::new();
    for ((from, to, _, _, _), support_pairs) in &candidates {
        *from_support.entry(from.clone()).or_insert(0) += *support_pairs;
        *to_support.entry(to.clone()).or_insert(0) += *support_pairs;
        *pair_support.entry((from.clone(), to.clone())).or_insert(0) += *support_pairs;
    }
    stats.candidates = candidates
        .into_iter()
        .map(
            |((from, to, orientation, distance, existing_dbg_edge), support_pairs)| {
                let same_from = from_support.get(&from).copied().unwrap_or(0);
                let same_to = to_support.get(&to).copied().unwrap_or(0);
                let same_pair = pair_support
                    .get(&(from.clone(), to.clone()))
                    .copied()
                    .unwrap_or(0);
                let conflict_pairs = same_from
                    .saturating_add(same_to)
                    .saturating_sub(same_pair * 2);
                let score = score_candidate(support_pairs, conflict_pairs, distance.as_ref());
                MateBridgeCandidate {
                    from_node: node_label(&from),
                    to_node: node_label(&to),
                    orientation,
                    distance,
                    support_pairs,
                    conflict_pairs,
                    score,
                    existing_dbg_edge,
                }
            },
        )
        .collect();
    stats
}

#[cfg(test)]
mod tests {
    use super::{boost_mate_pairs_on_existing_dbg_edges, MateBridgeStats, MatePairOrientation};
    use crate::dbg::graph::DbgGraph;
    use crate::evidence::{EvidenceKind, EvidenceSourceStage};
    use crate::illumina::read::Read;
    use std::collections::BTreeMap;

    #[test]
    fn mate_bridge_stats_builds_evidence_record() {
        let stats = MateBridgeStats {
            pairs_seen: 6,
            pairs_with_endpoint_kmers: 5,
            trusted_endpoint_pairs: 4,
            existing_edge_pairs: 3,
            report_only_pairs: 1,
            boosted_edges: 2,
            candidates: Vec::new(),
        };

        let record = stats.evidence_record();
        assert_eq!(record.kind, EvidenceKind::MateBridgeExistingEdge);
        assert_eq!(record.source_stage, EvidenceSourceStage::Phase2MateBridge);
        assert_eq!(record.support.observed, 6);
        assert_eq!(record.support.eligible, 4);
        assert_eq!(record.support.supporting, 3);
        assert_eq!(record.support.applied, 2);
        assert_eq!(record.counters["trusted_endpoint_pairs"], 4);
        assert_eq!(record.counters["report_only_pairs"], 1);
    }

    #[test]
    fn mate_bridge_stats_include_aggregated_candidates() {
        let mut graph = DbgGraph::new(
            4,
            BTreeMap::from([(b"ACGT".to_vec(), 2), (b"CGTA".to_vec(), 2)]),
        );
        graph
            .add_undirected_edge(b"ACGT", b"CGTA", 1)
            .expect("edge");
        let reads = vec![
            Read {
                id: "r1/1".into(),
                sequence: b"ACGT".to_vec(),
            },
            Read {
                id: "r2/1".into(),
                sequence: b"CGTA".to_vec(),
            },
        ];

        let stats =
            boost_mate_pairs_on_existing_dbg_edges(&mut graph, &reads, 1, 4, Some(12), Some(2));

        assert_eq!(stats.existing_edge_pairs, 1);
        assert_eq!(stats.report_only_pairs, 0);
        assert_eq!(stats.boosted_edges, 1);
        assert_eq!(stats.candidates.len(), 1);
        assert_eq!(stats.candidates[0].from_node, "ACGT");
        assert_eq!(stats.candidates[0].to_node, "CGTA");
        assert_eq!(
            stats.candidates[0].orientation,
            MatePairOrientation::R1TailToR2Head
        );
        assert_eq!(
            stats.candidates[0]
                .distance
                .as_ref()
                .map(|distance| distance.estimated_gap_bp),
            Some(4)
        );
        assert_eq!(
            stats.candidates[0]
                .distance
                .as_ref()
                .map(|distance| distance.confidence),
            Some(84)
        );
        assert_eq!(stats.candidates[0].support_pairs, 1);
        assert_eq!(stats.candidates[0].conflict_pairs, 0);
        assert_eq!(stats.candidates[0].score, 184);
        assert!(stats.candidates[0].existing_dbg_edge);
    }

    #[test]
    fn mate_bridge_keeps_absent_edges_as_report_only_candidates() {
        let mut graph = DbgGraph::new(
            4,
            BTreeMap::from([(b"ACGT".to_vec(), 2), (b"AAAA".to_vec(), 2)]),
        );
        let reads = vec![
            Read {
                id: "r1/1".into(),
                sequence: b"ACGT".to_vec(),
            },
            Read {
                id: "r2/1".into(),
                sequence: b"TTTT".to_vec(),
            },
        ];

        let stats =
            boost_mate_pairs_on_existing_dbg_edges(&mut graph, &reads, 1, 4, Some(10), None);

        assert_eq!(stats.trusted_endpoint_pairs, 1);
        assert_eq!(stats.existing_edge_pairs, 0);
        assert_eq!(stats.report_only_pairs, 1);
        assert_eq!(stats.boosted_edges, 0);
        assert_eq!(stats.candidates.len(), 1);
        assert_eq!(stats.candidates[0].from_node, "ACGT");
        assert_eq!(stats.candidates[0].to_node, "AAAA");
        assert_eq!(
            stats.candidates[0]
                .distance
                .as_ref()
                .map(|distance| distance.estimated_gap_bp),
            Some(2)
        );
        assert_eq!(
            stats.candidates[0]
                .distance
                .as_ref()
                .map(|distance| distance.confidence),
            Some(50)
        );
        assert_eq!(stats.candidates[0].score, 150);
        assert!(!stats.candidates[0].existing_dbg_edge);
        assert_eq!(graph.degree(b"ACGT"), 0);
    }

    #[test]
    fn mate_bridge_counts_competing_endpoint_conflicts() {
        let mut graph = DbgGraph::new(
            4,
            BTreeMap::from([
                (b"ACGT".to_vec(), 2),
                (b"AAAA".to_vec(), 2),
                (b"AAAC".to_vec(), 2),
            ]),
        );
        let reads = vec![
            Read {
                id: "r1/1".into(),
                sequence: b"ACGT".to_vec(),
            },
            Read {
                id: "r2/1".into(),
                sequence: b"ACGT".to_vec(),
            },
            Read {
                id: "r1/2".into(),
                sequence: b"TTTT".to_vec(),
            },
            Read {
                id: "r2/2".into(),
                sequence: b"AAAC".to_vec(),
            },
        ];

        let stats =
            boost_mate_pairs_on_existing_dbg_edges(&mut graph, &reads, 2, 4, Some(10), None);

        assert_eq!(stats.report_only_pairs, 2);
        assert_eq!(stats.candidates.len(), 2);
        assert!(stats.candidates.iter().all(|candidate| {
            candidate.from_node == "ACGT"
                && candidate.support_pairs == 1
                && candidate.conflict_pairs == 1
                && candidate.score == 50
        }));
    }
}
