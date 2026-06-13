//! Scaffold path construction from promoted mate-pair endpoint joins.

use crate::illumina::mate::MateBridgeCandidate;
use crate::illumina::scaffold::{EndpointJoinCandidate, ScaffoldLink, ScaffoldPath, ScaffoldStep};

pub struct ScaffoldPathBuilder<'a> {
    unitig_paths: &'a [Vec<Vec<u8>>],
    overlap_cigar: String,
}

impl<'a> ScaffoldPathBuilder<'a> {
    pub fn new(unitig_paths: &'a [Vec<Vec<u8>>], k: usize) -> Self {
        Self {
            unitig_paths,
            overlap_cigar: format!("{}M", k.saturating_sub(1)),
        }
    }

    pub fn path_for_candidate(
        &self,
        candidate: &MateBridgeCandidate,
        accepted_join: Option<&EndpointJoinCandidate>,
        id: String,
    ) -> Option<ScaffoldPath> {
        if !candidate.existing_dbg_edge && accepted_join.is_none() {
            return None;
        }
        let (from_idx, from_orient, to_idx, to_orient) = if candidate.existing_dbg_edge {
            self.candidate_unitig_tail_head(candidate)?
        } else {
            self.candidate_unitig_endpoint_steps(accepted_join?)?
        };
        if from_idx == to_idx {
            return None;
        }
        let from_segment = format!("utg{:06}", from_idx + 1);
        let to_segment = format!("utg{:06}", to_idx + 1);
        Some(ScaffoldPath {
            id,
            steps: vec![
                ScaffoldStep {
                    segment: from_segment.clone(),
                    orient: from_orient,
                },
                ScaffoldStep {
                    segment: to_segment.clone(),
                    orient: to_orient,
                },
            ],
            links: vec![ScaffoldLink {
                constraint_id: candidate.constraint_id.clone(),
                from_segment,
                from_orient,
                to_segment,
                to_orient,
                from_context: candidate.from_context.clone(),
                to_context: candidate.to_context.clone(),
                orientation: candidate.orientation,
                distance: candidate.distance.clone(),
                distance_bin: candidate.distance_bin.clone(),
                overlap_cigar: self.overlap_cigar.clone(),
                support_pairs: candidate.support_pairs,
                conflict_pairs: candidate.conflict_pairs,
                support_histogram: candidate.support_histogram.clone(),
                blockers: candidate.blockers.clone(),
                source: if candidate.existing_dbg_edge {
                    "mate_bridge_existing_edge"
                } else {
                    "mate_pair_endpoint_join_promoted_sidecar"
                }
                .to_string(),
            }],
        })
    }

    fn candidate_unitig_tail_head(
        &self,
        candidate: &MateBridgeCandidate,
    ) -> Option<(usize, char, usize, char)> {
        let mut from_idx = None;
        let mut to_idx = None;
        for (idx, path) in self.unitig_paths.iter().enumerate() {
            if endpoint_node_matches(path.last(), &candidate.from_node) {
                from_idx = Some(idx);
            }
            if endpoint_node_matches(path.first(), &candidate.to_node) {
                to_idx = Some(idx);
            }
        }
        Some((from_idx?, '+', to_idx?, '+'))
    }

    fn candidate_unitig_endpoint_steps(
        &self,
        join: &EndpointJoinCandidate,
    ) -> Option<(usize, char, usize, char)> {
        let mut from_step = None;
        let mut to_step = None;
        for (idx, path) in self.unitig_paths.iter().enumerate() {
            let first = path.first();
            let last = path.last();
            if join.from_side == "left" && endpoint_node_matches(first, &join.from_node) {
                from_step = Some((idx, '-'));
            }
            if join.from_side == "right" && endpoint_node_matches(last, &join.from_node) {
                from_step = Some((idx, '+'));
            }
            if join.to_side == "left" && endpoint_node_matches(first, &join.to_node) {
                to_step = Some((idx, '+'));
            }
            if join.to_side == "right" && endpoint_node_matches(last, &join.to_node) {
                to_step = Some((idx, '-'));
            }
        }
        let (from_idx, from_orient) = from_step?;
        let (to_idx, to_orient) = to_step?;
        Some((from_idx, from_orient, to_idx, to_orient))
    }
}

fn endpoint_node_matches(node: Option<&Vec<u8>>, expected: &str) -> bool {
    node.map(|node| node.as_slice() == expected.as_bytes())
        .unwrap_or(false)
}
