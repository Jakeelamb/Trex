//! Evidence-backed scaffold/path sidecar artifacts for Phase-2 Illumina.

use serde::{Deserialize, Serialize};

use crate::illumina::fragmentation::{
    FragmentEndpointReport, FragmentStopReason, FragmentationReport,
};
use crate::illumina::mate::{
    MateBridgeCandidate, MateDistanceBin, MateDistanceEvidence, MateGraphContext,
    MatePairOrientation, MateSupportHistogram,
};
use crate::illumina::promotion::{EndpointJoinPromotionPolicy, PromotionPolicySnapshot};
use crate::kmer::reverse_complement;

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

pub fn scaffold_fasta_records(
    artifact: &ScaffoldArtifact,
    unitig_records: &[(String, Vec<u8>)],
) -> Vec<(String, Vec<u8>)> {
    artifact
        .paths
        .iter()
        .filter_map(|path| {
            scaffold_path_sequence(path, unitig_records).map(|seq| (path.id.clone(), seq))
        })
        .collect()
}

pub fn scaffold_gfa_paths(artifact: &ScaffoldArtifact) -> Vec<(String, Vec<(usize, char)>)> {
    artifact
        .paths
        .iter()
        .filter_map(|path| {
            let mut steps = Vec::with_capacity(path.steps.len());
            for step in &path.steps {
                let idx = parse_utg_segment(&step.segment)?;
                if !matches!(step.orient, '+' | '-') {
                    return None;
                }
                steps.push((idx + 1, step.orient));
            }
            (!steps.is_empty()).then_some((path.id.clone(), steps))
        })
        .collect()
}

fn scaffold_path_sequence(
    path: &ScaffoldPath,
    unitig_records: &[(String, Vec<u8>)],
) -> Option<Vec<u8>> {
    let first = path.steps.first()?;
    let mut seq = oriented_segment_sequence(first, unitig_records)?;
    for idx in 1..path.steps.len() {
        let step = &path.steps[idx];
        let link = path.links.get(idx - 1)?;
        let mut next = oriented_segment_sequence(step, unitig_records)?;
        if link.source == "mate_bridge_existing_edge" {
            let overlap = overlap_cigar_bases(&link.overlap_cigar);
            if overlap < next.len() {
                next = next[overlap..].to_vec();
            } else {
                next.clear();
            }
        } else if let Some(gap) = positive_gap_len(link) {
            seq.extend(std::iter::repeat(b'N').take(gap));
        }
        seq.extend(next);
    }
    Some(seq)
}

fn oriented_segment_sequence(
    step: &ScaffoldStep,
    unitig_records: &[(String, Vec<u8>)],
) -> Option<Vec<u8>> {
    let idx = parse_utg_segment(&step.segment)?;
    let seq = unitig_records.get(idx)?.1.clone();
    match step.orient {
        '+' => Some(seq),
        '-' => Some(reverse_complement(&seq)),
        _ => None,
    }
}

fn parse_utg_segment(segment: &str) -> Option<usize> {
    let suffix = segment.strip_prefix("utg")?;
    let one_based = suffix.parse::<usize>().ok()?;
    one_based.checked_sub(1)
}

fn positive_gap_len(link: &ScaffoldLink) -> Option<usize> {
    let gap = link.distance.as_ref()?.estimated_gap_bp;
    (gap > 0).then_some(gap as usize)
}

fn overlap_cigar_bases(cigar: &str) -> usize {
    cigar
        .strip_suffix('M')
        .and_then(|digits| digits.parse::<usize>().ok())
        .unwrap_or(0)
}

pub fn build_scaffold_artifact(
    candidates: &[MateBridgeCandidate],
    unitig_paths: &[Vec<Vec<u8>>],
    k: usize,
    fragmentation: &FragmentationReport,
) -> ScaffoldArtifact {
    let mut paths = Vec::new();
    let overlap_cigar = format!("{}M", k.saturating_sub(1));

    let promotion_policy = PromotionPolicySnapshot::default();
    let endpoint_join_candidates =
        ranked_endpoint_joins(candidates, fragmentation, &promotion_policy.endpoint_join);
    for candidate in candidates {
        let accepted_join = accepted_endpoint_join(&endpoint_join_candidates, candidate);
        if !candidate.existing_dbg_edge && accepted_join.is_none() {
            continue;
        }
        let Some((from_idx, from_orient, to_idx, to_orient)) = (if candidate.existing_dbg_edge {
            candidate_unitig_tail_head(candidate, unitig_paths)
        } else {
            accepted_join.and_then(|join| candidate_unitig_endpoint_steps(join, unitig_paths))
        }) else {
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
                overlap_cigar: overlap_cigar.clone(),
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
        });
    }

    ScaffoldArtifact {
        schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
        promotion_policy,
        bridge_candidates: candidates.to_vec(),
        endpoint_join_candidates,
        paths,
    }
}

fn ranked_endpoint_joins(
    candidates: &[MateBridgeCandidate],
    fragmentation: &FragmentationReport,
    policy: &EndpointJoinPromotionPolicy,
) -> Vec<EndpointJoinCandidate> {
    let mut endpoint_use: std::collections::BTreeMap<&str, usize> =
        std::collections::BTreeMap::new();
    for candidate in candidates {
        if candidate.existing_dbg_edge {
            continue;
        }
        *endpoint_use
            .entry(candidate.from_node.as_str())
            .or_insert(0) += 1;
        *endpoint_use.entry(candidate.to_node.as_str()).or_insert(0) += 1;
    }
    let mut joins: Vec<EndpointJoinCandidate> = candidates
        .iter()
        .filter_map(|candidate| {
            let from = dead_end_endpoint_for_node(fragmentation, &candidate.from_node)?;
            let to = dead_end_endpoint_for_node(fragmentation, &candidate.to_node)?;
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
            let decision = policy.evaluate(candidate, conflict_cluster_size);
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

fn accepted_endpoint_join<'a>(
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
) -> Option<(usize, char, usize, char)> {
    let mut from_idx = None;
    let mut to_idx = None;
    for (idx, path) in unitig_paths.iter().enumerate() {
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
    join: &EndpointJoinCandidate,
    unitig_paths: &[Vec<Vec<u8>>],
) -> Option<(usize, char, usize, char)> {
    let mut from_step = None;
    let mut to_step = None;
    for (idx, path) in unitig_paths.iter().enumerate() {
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

fn endpoint_node_matches(node: Option<&Vec<u8>>, expected: &str) -> bool {
    node.map(|node| node.as_slice() == expected.as_bytes())
        .unwrap_or(false)
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
