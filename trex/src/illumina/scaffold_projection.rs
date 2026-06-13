//! FASTA and GFA projections for scaffold artifacts.

use crate::illumina::scaffold::{ScaffoldArtifact, ScaffoldLink, ScaffoldPath, ScaffoldStep};
use crate::kmer::reverse_complement;

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
