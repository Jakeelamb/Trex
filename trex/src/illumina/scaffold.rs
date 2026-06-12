//! Evidence-backed scaffold/path sidecar artifacts for Phase-2 Illumina.

use serde::{Deserialize, Serialize};

use crate::illumina::fragmentation::{
    FragmentEndpointReport, FragmentStopReason, FragmentationReport,
};
use crate::illumina::mate::{MateBridgeCandidate, MateDistanceEvidence, MatePairOrientation};

pub const SCAFFOLD_ARTIFACT_SCHEMA_VERSION: u64 = 3;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldStep {
    pub segment: String,
    pub orient: char,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldLink {
    pub from_segment: String,
    pub from_orient: char,
    pub to_segment: String,
    pub to_orient: char,
    pub orientation: MatePairOrientation,
    pub distance: Option<MateDistanceEvidence>,
    pub overlap_cigar: String,
    pub support_pairs: usize,
    pub conflict_pairs: usize,
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
    pub from_contig: String,
    pub from_side: String,
    pub from_node: String,
    pub to_contig: String,
    pub to_side: String,
    pub to_node: String,
    pub orientation: MatePairOrientation,
    pub distance: Option<MateDistanceEvidence>,
    pub support_pairs: usize,
    pub conflict_pairs: usize,
    pub score: u64,
    pub existing_dbg_edge: bool,
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldArtifact {
    pub schema_version: u64,
    pub bridge_candidates: Vec<MateBridgeCandidate>,
    pub endpoint_join_candidates: Vec<EndpointJoinCandidate>,
    pub paths: Vec<ScaffoldPath>,
}

impl ScaffoldArtifact {
    pub fn empty() -> Self {
        Self {
            schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
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
    let overlap_cigar = format!("{}M", k.saturating_sub(1));

    for candidate in candidates {
        if !candidate.existing_dbg_edge {
            continue;
        }
        let Some((from_idx, to_idx)) = candidate_unitig_tail_head(candidate, unitig_paths) else {
            continue;
        };
        if from_idx == to_idx {
            continue;
        }
        let from_segment = format!("utg{:06}", from_idx + 1);
        let to_segment = format!("utg{:06}", to_idx + 1);
        let id = format!("scf{:06}", paths.len() + 1);
        paths.push(ScaffoldPath {
            id,
            steps: vec![
                ScaffoldStep {
                    segment: from_segment.clone(),
                    orient: '+',
                },
                ScaffoldStep {
                    segment: to_segment.clone(),
                    orient: '+',
                },
            ],
            links: vec![ScaffoldLink {
                from_segment,
                from_orient: '+',
                to_segment,
                to_orient: '+',
                orientation: candidate.orientation,
                distance: candidate.distance.clone(),
                overlap_cigar: overlap_cigar.clone(),
                support_pairs: candidate.support_pairs,
                conflict_pairs: candidate.conflict_pairs,
                source: "mate_bridge_existing_edge".to_string(),
            }],
        });
    }
    let endpoint_join_candidates = ranked_endpoint_joins(candidates, fragmentation);

    ScaffoldArtifact {
        schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
        bridge_candidates: candidates.to_vec(),
        endpoint_join_candidates,
        paths,
    }
}

fn ranked_endpoint_joins(
    candidates: &[MateBridgeCandidate],
    fragmentation: &FragmentationReport,
) -> Vec<EndpointJoinCandidate> {
    let mut joins: Vec<EndpointJoinCandidate> = candidates
        .iter()
        .filter_map(|candidate| {
            let from = dead_end_endpoint_for_node(fragmentation, &candidate.from_node)?;
            let to = dead_end_endpoint_for_node(fragmentation, &candidate.to_node)?;
            if from.contig == to.contig {
                return None;
            }
            Some(EndpointJoinCandidate {
                id: String::new(),
                from_contig: from.contig.clone(),
                from_side: from.side.clone(),
                from_node: candidate.from_node.clone(),
                to_contig: to.contig.clone(),
                to_side: to.side.clone(),
                to_node: candidate.to_node.clone(),
                orientation: candidate.orientation,
                distance: candidate.distance.clone(),
                support_pairs: candidate.support_pairs,
                conflict_pairs: candidate.conflict_pairs,
                score: candidate.score,
                existing_dbg_edge: candidate.existing_dbg_edge,
                source: if candidate.existing_dbg_edge {
                    "mate_bridge_existing_edge"
                } else {
                    "mate_pair_endpoint_join_report_only"
                }
                .to_string(),
            })
        })
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

fn dead_end_endpoint_for_node<'a>(
    fragmentation: &'a FragmentationReport,
    node: &str,
) -> Option<&'a FragmentEndpointReport> {
    fragmentation.endpoints.iter().find(|endpoint| {
        endpoint.node.as_deref() == Some(node)
            && endpoint.reason == FragmentStopReason::GraphDeadEnd
    })
}

fn candidate_unitig_tail_head(
    candidate: &MateBridgeCandidate,
    unitig_paths: &[Vec<Vec<u8>>],
) -> Option<(usize, usize)> {
    let mut from_idx = None;
    let mut to_idx = None;
    for (idx, path) in unitig_paths.iter().enumerate() {
        if path
            .last()
            .map(|node| String::from_utf8_lossy(node) == candidate.from_node)
            .unwrap_or(false)
        {
            from_idx = Some(idx);
        }
        if path
            .first()
            .map(|node| String::from_utf8_lossy(node) == candidate.to_node)
            .unwrap_or(false)
        {
            to_idx = Some(idx);
        }
    }
    Some((from_idx?, to_idx?))
}

#[cfg(test)]
mod tests {
    use super::build_scaffold_artifact;
    use crate::illumina::fragmentation::{
        FragmentEndpointReport, FragmentStopReason, FragmentationReport, FragmentationSummary,
        FRAGMENTATION_REPORT_SCHEMA_VERSION,
    };
    use crate::illumina::mate::{MateBridgeCandidate, MateDistanceEvidence, MatePairOrientation};

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
        let distance = Some(MateDistanceEvidence {
            insert_mean_bp: 10,
            insert_stddev_bp: Some(2),
            read_bases_bp: 6,
            estimated_gap_bp: 4,
            confidence: 80,
        });
        MateBridgeCandidate {
            from_node: from_node.to_string(),
            to_node: to_node.to_string(),
            orientation: MatePairOrientation::R1TailToR2Head,
            distance,
            support_pairs,
            conflict_pairs: 0,
            score: (support_pairs as u64) * 100 + 80,
            existing_dbg_edge,
        }
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
        assert_eq!(artifact.endpoint_join_candidates.len(), 1);
        assert_eq!(artifact.paths.len(), 1);
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
            artifact.endpoint_join_candidates[0].source,
            "mate_pair_endpoint_join_report_only"
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
        assert!(artifact.paths.is_empty());
    }
}
