//! Promotion policy for moving evidence into assembler outputs.
//!
//! This module is the seam between report-only evidence and behavior that emits
//! scaffold artifacts, GFA paths, FASTA scaffolds, graph edits, or future
//! polishing edits.

use serde::{Deserialize, Serialize};

use crate::illumina::mate::MateBridgeCandidate;

pub const ENDPOINT_JOIN_MIN_SUPPORT_PAIRS: usize = 2;
pub const ENDPOINT_JOIN_MIN_DISTANCE_CONFIDENCE: u64 = 50;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PromotionStage {
    ReportOnlyCandidate,
    ScaffoldArtifact,
    GfaPath,
    FastaScaffoldWithGaps,
    GraphEdit,
    PolishingEdit,
}

impl PromotionStage {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ReportOnlyCandidate => "report_only_candidate",
            Self::ScaffoldArtifact => "scaffold_artifact",
            Self::GfaPath => "gfa_path",
            Self::FastaScaffoldWithGaps => "fasta_scaffold_with_gaps",
            Self::GraphEdit => "graph_edit",
            Self::PolishingEdit => "polishing_edit",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PromotionDecision {
    pub target_stage: PromotionStage,
    pub accepted: bool,
    pub rejection_reason: Option<&'static str>,
}

impl PromotionDecision {
    fn accepted(target_stage: PromotionStage) -> Self {
        Self {
            target_stage,
            accepted: true,
            rejection_reason: None,
        }
    }

    fn rejected(reason: &'static str) -> Self {
        Self {
            target_stage: PromotionStage::ReportOnlyCandidate,
            accepted: false,
            rejection_reason: Some(reason),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EndpointJoinPromotionPolicy {
    pub min_support_pairs: usize,
    pub min_distance_confidence: u64,
    pub reject_conflicting_pairs: bool,
    pub reject_competing_endpoint_clusters: bool,
}

impl Default for EndpointJoinPromotionPolicy {
    fn default() -> Self {
        Self {
            min_support_pairs: ENDPOINT_JOIN_MIN_SUPPORT_PAIRS,
            min_distance_confidence: ENDPOINT_JOIN_MIN_DISTANCE_CONFIDENCE,
            reject_conflicting_pairs: true,
            reject_competing_endpoint_clusters: true,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct PromotionPolicySnapshot {
    pub endpoint_join: EndpointJoinPromotionPolicy,
}

impl EndpointJoinPromotionPolicy {
    pub fn evaluate(
        &self,
        candidate: &MateBridgeCandidate,
        conflict_cluster_size: usize,
    ) -> PromotionDecision {
        if candidate.existing_dbg_edge {
            return PromotionDecision::accepted(PromotionStage::GraphEdit);
        }
        if self.reject_conflicting_pairs && candidate.conflict_pairs > 0 {
            return PromotionDecision::rejected("candidate has conflicting mate-pair support");
        }
        if self.reject_competing_endpoint_clusters && conflict_cluster_size > 1 {
            return PromotionDecision::rejected(
                "endpoint participates in a competing join cluster",
            );
        }
        if candidate.support_pairs < self.min_support_pairs {
            return PromotionDecision::rejected("candidate has fewer than two supporting pairs");
        }
        let confidence = candidate
            .distance
            .as_ref()
            .map(|distance| distance.confidence)
            .unwrap_or(0);
        if confidence < self.min_distance_confidence {
            return PromotionDecision::rejected(
                "candidate distance confidence below promotion threshold",
            );
        }
        PromotionDecision::accepted(PromotionStage::ScaffoldArtifact)
    }
}

#[cfg(test)]
mod tests {
    use super::{EndpointJoinPromotionPolicy, PromotionStage};
    use crate::illumina::mate::{
        MateBridgeCandidate, MateBridgeCandidateParts, MateDistanceEvidence, MatePairOrientation,
    };

    fn candidate(
        support_pairs: usize,
        conflict_pairs: usize,
        confidence: u64,
        existing_dbg_edge: bool,
    ) -> MateBridgeCandidate {
        MateBridgeCandidate::from_constraint_parts(MateBridgeCandidateParts {
            constraint_id: "kbm000001".to_string(),
            from_node: "AAC".to_string(),
            to_node: "ACC".to_string(),
            orientation: MatePairOrientation::R1TailToR2Head,
            distance: Some(MateDistanceEvidence {
                insert_mean_bp: 10,
                insert_stddev_bp: None,
                read_bases_bp: 8,
                estimated_gap_bp: 2,
                confidence,
            }),
            support_pairs,
            same_from_pairs: support_pairs,
            same_to_pairs: support_pairs,
            same_pair_pairs: support_pairs,
            conflict_pairs,
            existing_dbg_edge,
        })
    }

    #[test]
    fn accepts_minimum_supported_endpoint_join() {
        let policy = EndpointJoinPromotionPolicy::default();

        let decision = policy.evaluate(&candidate(2, 0, 50, false), 1);

        assert!(decision.accepted);
        assert_eq!(decision.target_stage, PromotionStage::ScaffoldArtifact);
        assert_eq!(decision.rejection_reason, None);
    }

    #[test]
    fn rejects_in_precedence_order() {
        let policy = EndpointJoinPromotionPolicy::default();

        let conflicting = policy.evaluate(&candidate(1, 1, 49, false), 2);
        let cluster = policy.evaluate(&candidate(1, 0, 49, false), 2);
        let weak_support = policy.evaluate(&candidate(1, 0, 49, false), 1);
        let weak_confidence = policy.evaluate(&candidate(2, 0, 49, false), 1);

        assert_eq!(
            conflicting.rejection_reason,
            Some("candidate has conflicting mate-pair support")
        );
        assert_eq!(
            cluster.rejection_reason,
            Some("endpoint participates in a competing join cluster")
        );
        assert_eq!(
            weak_support.rejection_reason,
            Some("candidate has fewer than two supporting pairs")
        );
        assert_eq!(
            weak_confidence.rejection_reason,
            Some("candidate distance confidence below promotion threshold")
        );
    }

    #[test]
    fn treats_existing_dbg_edge_as_graph_edit_evidence() {
        let policy = EndpointJoinPromotionPolicy::default();

        let decision = policy.evaluate(&candidate(1, 0, 0, true), 0);

        assert!(decision.accepted);
        assert_eq!(decision.target_stage, PromotionStage::GraphEdit);
    }
}
