//! Evidence-backed scaffold/path sidecar artifacts for Phase-2 Illumina.

use serde::{Deserialize, Serialize};

use crate::illumina::mate::MateBridgeCandidate;

pub const SCAFFOLD_ARTIFACT_SCHEMA_VERSION: u64 = 1;

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
    pub overlap_cigar: String,
    pub support_pairs: usize,
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldPath {
    pub id: String,
    pub steps: Vec<ScaffoldStep>,
    pub links: Vec<ScaffoldLink>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaffoldArtifact {
    pub schema_version: u64,
    pub bridge_candidates: Vec<MateBridgeCandidate>,
    pub paths: Vec<ScaffoldPath>,
}

impl ScaffoldArtifact {
    pub fn empty() -> Self {
        Self {
            schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
            bridge_candidates: Vec::new(),
            paths: Vec::new(),
        }
    }
}

pub fn build_scaffold_artifact(
    candidates: &[MateBridgeCandidate],
    unitig_paths: &[Vec<Vec<u8>>],
    k: usize,
) -> ScaffoldArtifact {
    let mut paths = Vec::new();
    let overlap_cigar = format!("{}M", k.saturating_sub(1));

    for candidate in candidates {
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
                overlap_cigar: overlap_cigar.clone(),
                support_pairs: candidate.support_pairs,
                source: "mate_bridge_existing_edge".to_string(),
            }],
        });
    }

    ScaffoldArtifact {
        schema_version: SCAFFOLD_ARTIFACT_SCHEMA_VERSION,
        bridge_candidates: candidates.to_vec(),
        paths,
    }
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
    use crate::illumina::mate::MateBridgeCandidate;

    #[test]
    fn builds_path_for_explicit_unitig_tail_to_head_candidate() {
        let candidates = vec![MateBridgeCandidate {
            from_node: "AAC".to_string(),
            to_node: "ACC".to_string(),
            support_pairs: 3,
            score: 3,
            existing_dbg_edge: true,
        }];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3);

        assert_eq!(artifact.bridge_candidates, candidates);
        assert_eq!(artifact.paths.len(), 1);
        assert_eq!(artifact.paths[0].steps[0].segment, "utg000001");
        assert_eq!(artifact.paths[0].steps[1].segment, "utg000002");
        assert_eq!(artifact.paths[0].links[0].overlap_cigar, "2M");
        assert_eq!(artifact.paths[0].links[0].support_pairs, 3);
    }

    #[test]
    fn does_not_invent_path_when_candidate_is_not_unitig_boundary() {
        let candidates = vec![MateBridgeCandidate {
            from_node: "AAA".to_string(),
            to_node: "CCC".to_string(),
            support_pairs: 1,
            score: 1,
            existing_dbg_edge: true,
        }];
        let unitigs = vec![
            vec![b"AAA".to_vec(), b"AAC".to_vec()],
            vec![b"ACC".to_vec(), b"CCC".to_vec()],
        ];

        let artifact = build_scaffold_artifact(&candidates, &unitigs, 3);

        assert_eq!(artifact.bridge_candidates.len(), 1);
        assert!(artifact.paths.is_empty());
    }
}
