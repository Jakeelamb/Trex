//! Promotion engine for mate-pair endpoint joins in scaffold sidecars.

use std::collections::BTreeMap;

use crate::illumina::fragmentation::{
    FragmentEndpointReport, FragmentStopReason, FragmentationReport,
};
use crate::illumina::mate::MateBridgeCandidate;
use crate::illumina::promotion::EndpointJoinPromotionPolicy;
use crate::illumina::scaffold::EndpointJoinCandidate;

pub struct ScaffoldPromotionEngine<'a> {
    candidates: &'a [MateBridgeCandidate],
    fragmentation: &'a FragmentationReport,
    policy: &'a EndpointJoinPromotionPolicy,
}

impl<'a> ScaffoldPromotionEngine<'a> {
    pub fn new(
        candidates: &'a [MateBridgeCandidate],
        fragmentation: &'a FragmentationReport,
        policy: &'a EndpointJoinPromotionPolicy,
    ) -> Self {
        Self {
            candidates,
            fragmentation,
            policy,
        }
    }

    pub fn ranked_endpoint_joins(&self) -> Vec<EndpointJoinCandidate> {
        let endpoint_use = self.endpoint_use_counts();
        let mut joins: Vec<EndpointJoinCandidate> = self
            .candidates
            .iter()
            .filter_map(|candidate| self.endpoint_join_candidate(candidate, &endpoint_use))
            .collect();
        joins.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then_with(|| b.support_pairs.cmp(&a.support_pairs))
                .then_with(|| a.from_contig.cmp(&b.from_contig))
                .then_with(|| a.from_side.cmp(&b.from_side))
                .then_with(|| a.to_contig.cmp(&b.to_contig))
                .then_with(|| a.to_side.cmp(&b.to_side))
                .then_with(|| a.from_node.cmp(&b.from_node))
                .then_with(|| a.to_node.cmp(&b.to_node))
        });
        for (idx, join) in joins.iter_mut().enumerate() {
            join.id = format!("ejc{:06}", idx + 1);
        }
        joins
    }

    fn endpoint_use_counts(&self) -> BTreeMap<&'a str, usize> {
        let mut endpoint_use = BTreeMap::new();
        for candidate in self.candidates {
            if candidate.existing_dbg_edge {
                continue;
            }
            *endpoint_use
                .entry(candidate.from_node.as_str())
                .or_insert(0) += 1;
            *endpoint_use.entry(candidate.to_node.as_str()).or_insert(0) += 1;
        }
        endpoint_use
    }

    fn endpoint_join_candidate(
        &self,
        candidate: &MateBridgeCandidate,
        endpoint_use: &BTreeMap<&str, usize>,
    ) -> Option<EndpointJoinCandidate> {
        let from = self.dead_end_endpoint_for_node(&candidate.from_node)?;
        let to = self.dead_end_endpoint_for_node(&candidate.to_node)?;
        if from.contig == to.contig {
            return None;
        }
        let conflict_cluster_size = if candidate.existing_dbg_edge {
            0
        } else {
            endpoint_use
                .get(candidate.from_node.as_str())
                .copied()
                .unwrap_or(0)
                .max(
                    endpoint_use
                        .get(candidate.to_node.as_str())
                        .copied()
                        .unwrap_or(0),
                )
        };
        let decision = self.policy.evaluate(candidate, conflict_cluster_size);
        let accepted = decision.accepted;
        let promotion_stage = decision.target_stage.as_str().to_string();
        let rejection_reason = decision.rejection_reason.map(str::to_string);
        Some(EndpointJoinCandidate {
            id: String::new(),
            constraint_id: candidate.constraint_id.clone(),
            from_contig: from.contig.clone(),
            from_side: from.side.clone(),
            from_node: candidate.from_node.clone(),
            to_contig: to.contig.clone(),
            to_side: to.side.clone(),
            to_node: candidate.to_node.clone(),
            from_context: candidate.from_context.clone(),
            to_context: candidate.to_context.clone(),
            orientation: candidate.orientation,
            distance: candidate.distance.clone(),
            distance_bin: candidate.distance_bin.clone(),
            support_pairs: candidate.support_pairs,
            conflict_pairs: candidate.conflict_pairs,
            support_histogram: candidate.support_histogram.clone(),
            conflict_cluster_size,
            score: candidate.score,
            existing_dbg_edge: candidate.existing_dbg_edge,
            blockers: candidate.blockers.clone(),
            promotion_stage,
            accepted,
            rejection_reason,
            source: if candidate.existing_dbg_edge {
                "mate_bridge_existing_edge"
            } else if accepted {
                "mate_pair_endpoint_join_promoted_sidecar"
            } else {
                "mate_pair_endpoint_join_rejected"
            }
            .to_string(),
        })
    }

    fn dead_end_endpoint_for_node(&self, node: &str) -> Option<&'a FragmentEndpointReport> {
        self.fragmentation.endpoints.iter().find(|endpoint| {
            endpoint.node.as_deref() == Some(node)
                && endpoint.reason == FragmentStopReason::GraphDeadEnd
        })
    }
}

pub fn accepted_endpoint_join<'a>(
    joins: &'a [EndpointJoinCandidate],
    candidate: &MateBridgeCandidate,
) -> Option<&'a EndpointJoinCandidate> {
    joins.iter().find(|join| {
        join.accepted
            && !join.existing_dbg_edge
            && join.from_node == candidate.from_node
            && join.to_node == candidate.to_node
            && join.orientation == candidate.orientation
    })
}
