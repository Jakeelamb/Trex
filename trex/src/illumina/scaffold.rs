//! Evidence-backed scaffold/path sidecar artifacts for Phase-2 Illumina.

use serde::{Deserialize, Serialize};

use crate::illumina::fragmentation::FragmentationReport;
use crate::illumina::mate::{
    MateBridgeCandidate, MateDistanceBin, MateDistanceEvidence, MateGraphContext,
    MatePairOrientation, MateSupportHistogram,
};
use crate::illumina::promotion::PromotionPolicySnapshot;
use crate::illumina::scaffold_path::ScaffoldPathBuilder;
pub use crate::illumina::scaffold_projection::{scaffold_fasta_records, scaffold_gfa_paths};
use crate::illumina::scaffold_promotion::{accepted_endpoint_join, ScaffoldPromotionEngine};

pub const SCAFFOLD_ARTIFACT_SCHEMA_VERSION: u64 = 6;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldStep {
    pub segment: String,
    pub orient: char,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldLink {
    pub constraint_id: String,
    pub from_segment: String,
    pub from_orient: char,
    pub to_segment: String,
    pub to_orient: char,
    pub from_context: MateGraphContext,
    pub to_context: MateGraphContext,
    pub orientation: MatePairOrientation,
    pub distance: Option<MateDistanceEvidence>,
    pub distance_bin: MateDistanceBin,
    pub overlap_cigar: String,
    pub support_pairs: usize,
    pub conflict_pairs: usize,
    pub support_histogram: MateSupportHistogram,
    pub blockers: Vec<String>,
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldPath {
    pub id: String,
    pub steps: Vec<ScaffoldStep>,
    pub links: Vec<ScaffoldLink>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EndpointJoinCandidate {
    pub id: String,
    pub constraint_id: String,
    pub from_contig: String,
    pub from_side: String,
    pub from_node: String,
    pub to_contig: String,
    pub to_side: String,
    pub to_node: String,
    pub from_context: MateGraphContext,
    pub to_context: MateGraphContext,
    pub orientation: MatePairOrientation,
    pub distance: Option<MateDistanceEvidence>,
    pub distance_bin: MateDistanceBin,
    pub support_pairs: usize,
    pub conflict_pairs: usize,
    pub support_histogram: MateSupportHistogram,
    pub conflict_cluster_size: usize,
    pub score: u64,
    pub existing_dbg_edge: bool,
    pub blockers: Vec<String>,
    pub promotion_stage: String,
    pub accepted: bool,
    pub rejection_reason: Option<String>,
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldArtifact {
    pub schema_version: u64,
    pub promotion_policy: PromotionPolicySnapshot,
    pub bridge_candidates: Vec<MateBridgeCandidate>,
    pub endpoint_join_candidates: Vec<EndpointJoinCandidate>,
    pub paths: Vec<ScaffoldPath>,
}

impl ScaffoldArtifact {
    pub fn empty() -> Self {
        Self {
            schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
            promotion_policy: PromotionPolicySnapshot::default(),
            bridge_candidates: Vec::new(),
            endpoint_join_candidates: Vec::new(),
            paths: Vec::new(),
        }
    }
}

pub fn build_scaffold_artifact(
    candidates: &[MateBridgeCandidate],
    unitig_paths: &[Vec<Vec<u8>>],
    k: usize,
    fragmentation: &FragmentationReport,
) -> ScaffoldArtifact {
    let mut paths = Vec::new();
    let path_builder = ScaffoldPathBuilder::new(unitig_paths, k);

    let promotion_policy = PromotionPolicySnapshot::default();
    let endpoint_join_candidates =
        ScaffoldPromotionEngine::new(candidates, fragmentation, &promotion_policy.endpoint_join)
            .ranked_endpoint_joins();
    for candidate in candidates {
        let accepted_join = accepted_endpoint_join(&endpoint_join_candidates, candidate);
        let id = format!("scf{:06}", paths.len() + 1);
        if let Some(path) = path_builder.path_for_candidate(candidate, accepted_join, id) {
            paths.push(path);
        }
    }

    ScaffoldArtifact {
        schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
        promotion_policy,
        bridge_candidates: candidates.to_vec(),
        endpoint_join_candidates,
        paths,
    }
}

#[cfg(test)]
mod tests {
    use super::{build_scaffold_artifact, scaffold_fasta_records, scaffold_gfa_paths};
    use crate::illumina::fragmentation::{
        FragmentEndpointReport, FragmentStopReason, FragmentationReport, FragmentationSummary,
        FRAGMENTATION_REPORT_SCHEMA_VERSION,
    };
    use crate::illumina::mate::{
        MateBridgeCandidate, MateBridgeCandidateParts, MateDistanceEvidence, MatePairOrientation,
    };

    fn fragmentation(endpoints: Vec<FragmentEndpointReport>) -> FragmentationReport {
        FragmentationReport {
            schema_version: FRAGMENTATION_REPORT_SCHEMA_VERSION,
            summary: FragmentationSummary {
                contigs: 2,
                contig_path_nodes: 4,
                endpoints: endpoints.len(),
                graph_dead_end_endpoints: endpoints
                    .iter()
                    .filter(|endpoint| endpoint.reason == FragmentStopReason::GraphDeadEnd)
                    .count(),
                ..FragmentationSummary::default()
            },
            endpoints,
        }
    }

    fn endpoint(contig: &str, side: &str, node: &str) -> FragmentEndpointReport {
        FragmentEndpointReport {
            contig: contig.to_string(),
            side: side.to_string(),
            node: Some(node.to_string()),
            degree: 1,
            multiplicity: 1,
            repeat_suspected: false,
            reason: FragmentStopReason::GraphDeadEnd,
        }
    }

    fn candidate(
        from_node: &str,
        to_node: &str,
        support_pairs: usize,
        existing_dbg_edge: bool,
    ) -> MateBridgeCandidate {
        candidate_with_conflicts(from_node, to_node, support_pairs, 0, existing_dbg_edge)
    }

    fn candidate_with_conflicts(
        from_node: &str,
        to_node: &str,
        support_pairs: usize,
        conflict_pairs: usize,
        existing_dbg_edge: bool,
    ) -> MateBridgeCandidate {
        let distance = Some(MateDistanceEvidence {
            insert_mean_bp: 10,
            insert_stddev_bp: Some(2),
            read_bases_bp: 6,
            estimated_gap_bp: 4,
            confidence: 80,
        });
        MateBridgeCandidate::from_constraint_parts(MateBridgeCandidateParts {
            constraint_id: "kbm000001".to_string(),
            from_node: from_node.to_string(),
            to_node: to_node.to_string(),
            orientation: MatePairOrientation::R1TailToR2Head,
            distance,
            support_pairs,
            same_from_pairs: support_pairs,
            same_to_pairs: support_pairs,
            same_pair_pairs: support_pairs,
            conflict_pairs,
            existing_dbg_edge,
        })
    }

    #[test]
    fn builds_path_for_explicit_unitig_tail_to_head_candidate() {
        let candidates = vec![candidate("AAC", "ACC", 3, true)];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];

        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
        ]);

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3, &frag);

        assert_eq!(artifact.bridge_candidates, candidates);
        assert_eq!(artifact.schema_version, 6);
        assert_eq!(artifact.endpoint_join_candidates.len(), 1);
        assert_eq!(artifact.paths.len(), 1);
        assert_eq!(artifact.paths[0].links[0].constraint_id, "kbm000001");
        assert_eq!(artifact.paths[0].links[0].from_context.side, "tail");
        assert_eq!(artifact.paths[0].links[0].to_context.side, "head");
        assert_eq!(artifact.paths[0].steps[0].segment, "utg000001");
        assert_eq!(artifact.paths[0].steps[1].segment, "utg000002");
        assert_eq!(artifact.paths[0].links[0].overlap_cigar, "2M");
        assert_eq!(artifact.paths[0].links[0].support_pairs, 3);
        assert_eq!(
            artifact.paths[0].links[0].orientation,
            MatePairOrientation::R1TailToR2Head
        );
        assert_eq!(
            artifact.paths[0].links[0]
                .distance
                .as_ref()
                .map(|distance| distance.estimated_gap_bp),
            Some(4)
        );
    }

    #[test]
    fn does_not_invent_path_when_candidate_is_not_unitig_boundary() {
        let candidates = vec![candidate("AAA", "CCC", 1, true)];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];

        let frag = fragmentation(vec![
            endpoint("ctg1", "left", "AAA"),
            endpoint("ctg2", "right", "CCC"),
        ]);

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3, &frag);

        assert_eq!(artifact.bridge_candidates.len(), 1);
        assert_eq!(artifact.endpoint_join_candidates.len(), 1);
        assert!(artifact.paths.is_empty());
    }

    #[test]
    fn ranks_report_only_dead_end_endpoint_joins_without_paths() {
        let candidates = vec![candidate("AAC", "ACC", 5, false)];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];
        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
        ]);

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3, &frag);

        assert_eq!(artifact.bridge_candidates.len(), 1);
        assert_eq!(artifact.endpoint_join_candidates.len(), 1);
        assert_eq!(artifact.endpoint_join_candidates[0].id, "ejc000001");
        assert_eq!(
            artifact.endpoint_join_candidates[0].constraint_id,
            "kbm000001"
        );
        assert_eq!(
            artifact.endpoint_join_candidates[0].from_context.node,
            "AAC"
        );
        assert_eq!(artifact.endpoint_join_candidates[0].to_context.node, "ACC");
        assert_eq!(
            artifact.endpoint_join_candidates[0]
                .support_histogram
                .support_pairs,
            5
        );
        assert!(artifact.endpoint_join_candidates[0]
            .blockers
            .contains(&"absent_dbg_edge_no_graph_edit".to_string()));
        assert!(artifact.endpoint_join_candidates[0].accepted);
        assert_eq!(
            artifact.endpoint_join_candidates[0].promotion_stage,
            "scaffold_artifact"
        );
        assert_eq!(artifact.promotion_policy.endpoint_join.min_support_pairs, 2);
        assert_eq!(
            artifact
                .promotion_policy
                .endpoint_join
                .min_distance_confidence,
            50
        );
        assert_eq!(
            artifact.endpoint_join_candidates[0].source,
            "mate_pair_endpoint_join_promoted_sidecar"
        );
        assert_eq!(
            artifact.endpoint_join_candidates[0].orientation,
            MatePairOrientation::R1TailToR2Head
        );
        assert_eq!(
            artifact.endpoint_join_candidates[0]
                .distance
                .as_ref()
                .map(|distance| distance.confidence),
            Some(80)
        );
        assert_eq!(artifact.paths.len(), 1);
        assert_eq!(
            artifact.paths[0].links[0].source,
            "mate_pair_endpoint_join_promoted_sidecar"
        );
    }

    #[test]
    fn rejects_endpoint_join_conflict_clusters_before_path_promotion() {
        let candidates = vec![
            candidate("AAC", "ACC", 5, false),
            candidate("AAC", "CCC", 4, false),
            candidate_with_conflicts("GGG", "TTT", 5, 1, false),
        ];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
            vec![b"GGG".to_vec()],
            vec![b"TTT".to_vec()],
        ];
        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
            endpoint("ctg3", "left", "CCC"),
            endpoint("ctg4", "right", "GGG"),
            endpoint("ctg5", "left", "TTT"),
        ]);

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3, &frag);

        assert_eq!(artifact.endpoint_join_candidates.len(), 3);
        assert!(artifact
            .endpoint_join_candidates
            .iter()
            .all(|candidate| !candidate.accepted));
        assert!(artifact
            .endpoint_join_candidates
            .iter()
            .all(|candidate| candidate.promotion_stage == "report_only_candidate"));
        assert!(artifact
            .endpoint_join_candidates
            .iter()
            .any(|candidate| candidate.rejection_reason.as_deref()
                == Some("endpoint participates in a competing join cluster")));
        assert!(artifact
            .endpoint_join_candidates
            .iter()
            .any(|candidate| candidate.rejection_reason.as_deref()
                == Some("candidate has conflicting mate-pair support")));
        assert!(artifact.paths.is_empty());
    }

    #[test]
    fn scaffold_fasta_trims_existing_dbg_overlap() {
        let candidates = vec![candidate("AAC", "ACC", 3, true)];
        let unitig_paths = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];
        let unitig_records = vec![
            ("utg000001".to_string(), b"AAAC".to_vec()),
            ("utg000002".to_string(), b"ACCC".to_vec()),
        ];
        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
        ]);
        let artifact = build_scaffold_artifact(&candidates, &unitig_paths, 3, &frag);

        let records = scaffold_fasta_records(&artifact, &unitig_records);

        assert_eq!(records, vec![("scf000001".to_string(), b"AAACCC".to_vec())]);
    }

    #[test]
    fn scaffold_fasta_uses_explicit_n_gap_for_promoted_endpoint_join() {
        let candidates = vec![candidate("AAC", "ACC", 5, false)];
        let unitig_paths = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];
        let unitig_records = vec![
            ("utg000001".to_string(), b"AAAC".to_vec()),
            ("utg000002".to_string(), b"ACCC".to_vec()),
        ];
        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
        ]);
        let artifact = build_scaffold_artifact(&candidates, &unitig_paths, 3, &frag);

        let records = scaffold_fasta_records(&artifact, &unitig_records);

        assert_eq!(
            records,
            vec![("scf000001".to_string(), b"AAACNNNNACCC".to_vec())]
        );
    }

    #[test]
    fn scaffold_gfa_paths_preserve_accepted_scaffold_steps() {
        let candidates = vec![candidate("AAC", "ACC", 5, false)];
        let unitig_paths = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];
        let frag = fragmentation(vec![
            endpoint("ctg1", "right", "AAC"),
            endpoint("ctg2", "left", "ACC"),
        ]);
        let artifact = build_scaffold_artifact(&candidates, &unitig_paths, 3, &frag);

        let paths = scaffold_gfa_paths(&artifact);

        assert_eq!(
            paths,
            vec![("scf000001".to_string(), vec![(1, '+'), (2, '+')])]
        );
    }
}
