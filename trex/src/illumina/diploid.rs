//! Diploid evidence summaries. These report parent-specific support without claiming full phasing.

use std::collections::{BTreeMap, HashSet};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{IngestError, TrexError};
use crate::illumina::fasta::parse_fasta;
use crate::illumina::io::read_maybe_gzip;
use crate::kmer::canonical_kmer;

pub const DIPLOID_EVIDENCE_SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone, Default)]
pub struct ParentReferenceParams {
    pub parent1: Option<std::path::PathBuf>,
    pub parent2: Option<std::path::PathBuf>,
}

impl ParentReferenceParams {
    pub fn has_pair(&self) -> bool {
        self.parent1.is_some() && self.parent2.is_some()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ParentAssignment {
    Parent1,
    Parent2,
    Shared,
    Mixed,
    Uninformative,
}

impl ParentAssignment {
    pub fn as_gfa_tag_value(self) -> &'static str {
        match self {
            ParentAssignment::Parent1 => "parent1",
            ParentAssignment::Parent2 => "parent2",
            ParentAssignment::Shared => "shared",
            ParentAssignment::Mixed => "mixed",
            ParentAssignment::Uninformative => "uninformative",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParentKmerSummary {
    pub parent1_only_kmers: usize,
    pub parent2_only_kmers: usize,
    pub shared_parent_kmers: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnitigParentEvidence {
    pub unitig_index: usize,
    pub node_count: usize,
    pub parent1_only_nodes: usize,
    pub parent2_only_nodes: usize,
    pub shared_parent_nodes: usize,
    pub uninformative_nodes: usize,
    pub assignment: ParentAssignment,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigParentEvidence {
    pub contig: String,
    pub parent1_only_kmers: usize,
    pub parent2_only_kmers: usize,
    pub shared_parent_kmers: usize,
    pub uninformative_kmers: usize,
    pub assignment: ParentAssignment,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiploidEvidenceSummary {
    pub diploid_enabled: bool,
    pub parent_references_supplied: bool,
    pub full_haplotype_fasta_claimed: bool,
    pub parent1_only_nodes: usize,
    pub parent2_only_nodes: usize,
    pub shared_parent_nodes: usize,
    pub mixed_unitigs: usize,
    pub parent_informative_unitigs: usize,
    pub parent_informative_contigs: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiploidEvidenceReport {
    pub schema_version: u64,
    pub summary: DiploidEvidenceSummary,
    pub parent_kmers: Option<ParentKmerSummary>,
    pub unitigs: Vec<UnitigParentEvidence>,
    pub contigs: Vec<ContigParentEvidence>,
}

#[derive(Debug)]
struct ParentKmers {
    parent1_only: HashSet<Vec<u8>>,
    parent2_only: HashSet<Vec<u8>>,
    shared: HashSet<Vec<u8>>,
}

pub fn build_diploid_evidence(
    enabled: bool,
    k: usize,
    parent_refs: &ParentReferenceParams,
    unitig_paths: &[Vec<Vec<u8>>],
    contig_records: &[(String, Vec<u8>)],
) -> Result<DiploidEvidenceReport, TrexError> {
    if !enabled {
        return Ok(empty_report(false, false));
    }
    let Some(parent1) = parent_refs.parent1.as_ref() else {
        return Ok(empty_report(true, false));
    };
    let Some(parent2) = parent_refs.parent2.as_ref() else {
        return Ok(empty_report(true, false));
    };

    let parent_kmers = load_parent_kmers(parent1, parent2, k)?;
    let unitigs = classify_unitigs(unitig_paths, &parent_kmers);
    let contigs = classify_contigs(contig_records, k, &parent_kmers);
    let summary = summarize(true, true, &unitigs, &contigs);
    Ok(DiploidEvidenceReport {
        schema_version: DIPLOID_EVIDENCE_SCHEMA_VERSION,
        summary,
        parent_kmers: Some(ParentKmerSummary {
            parent1_only_kmers: parent_kmers.parent1_only.len(),
            parent2_only_kmers: parent_kmers.parent2_only.len(),
            shared_parent_kmers: parent_kmers.shared.len(),
        }),
        unitigs,
        contigs,
    })
}

pub fn gfa_path_tags(report: &DiploidEvidenceReport) -> Vec<(String, String)> {
    report
        .contigs
        .iter()
        .map(|contig| {
            (
                contig.contig.clone(),
                format!("PS:Z:{}", contig.assignment.as_gfa_tag_value()),
            )
        })
        .collect()
}

pub fn gfa_unitig_tags(report: &DiploidEvidenceReport) -> Vec<(usize, String)> {
    report
        .unitigs
        .iter()
        .map(|unitig| {
            (
                unitig.unitig_index + 1,
                format!("PS:Z:{}", unitig.assignment.as_gfa_tag_value()),
            )
        })
        .collect()
}

fn empty_report(diploid_enabled: bool, parent_references_supplied: bool) -> DiploidEvidenceReport {
    DiploidEvidenceReport {
        schema_version: DIPLOID_EVIDENCE_SCHEMA_VERSION,
        summary: DiploidEvidenceSummary {
            diploid_enabled,
            parent_references_supplied,
            full_haplotype_fasta_claimed: false,
            parent1_only_nodes: 0,
            parent2_only_nodes: 0,
            shared_parent_nodes: 0,
            mixed_unitigs: 0,
            parent_informative_unitigs: 0,
            parent_informative_contigs: 0,
        },
        parent_kmers: None,
        unitigs: Vec::new(),
        contigs: Vec::new(),
    }
}

fn load_parent_kmers(parent1: &Path, parent2: &Path, k: usize) -> Result<ParentKmers, TrexError> {
    let p1 = fasta_kmers(parent1, k)?;
    let p2 = fasta_kmers(parent2, k)?;
    let parent1_only = p1.difference(&p2).cloned().collect();
    let parent2_only = p2.difference(&p1).cloned().collect();
    let shared = p1.intersection(&p2).cloned().collect();
    Ok(ParentKmers {
        parent1_only,
        parent2_only,
        shared,
    })
}

fn fasta_kmers(path: &Path, k: usize) -> Result<HashSet<Vec<u8>>, TrexError> {
    let bytes = read_maybe_gzip(path).map_err(|e| TrexError::Ingest(IngestError::Io(e)))?;
    let records = parse_fasta(&bytes)?;
    let mut out = HashSet::new();
    if k == 0 {
        return Ok(out);
    }
    for record in records {
        for window in record.sequence.windows(k) {
            if window
                .iter()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            {
                out.insert(canonical_kmer(window));
            }
        }
    }
    Ok(out)
}

fn classify_unitigs(
    unitig_paths: &[Vec<Vec<u8>>],
    parent_kmers: &ParentKmers,
) -> Vec<UnitigParentEvidence> {
    unitig_paths
        .iter()
        .enumerate()
        .map(|(idx, path)| {
            let mut counts = ParentCounts::default();
            for node in path {
                counts.observe(node, parent_kmers);
            }
            UnitigParentEvidence {
                unitig_index: idx,
                node_count: path.len(),
                parent1_only_nodes: counts.parent1_only,
                parent2_only_nodes: counts.parent2_only,
                shared_parent_nodes: counts.shared,
                uninformative_nodes: counts.uninformative,
                assignment: counts.assignment(),
            }
        })
        .collect()
}

fn classify_contigs(
    contig_records: &[(String, Vec<u8>)],
    k: usize,
    parent_kmers: &ParentKmers,
) -> Vec<ContigParentEvidence> {
    contig_records
        .iter()
        .enumerate()
        .map(|(idx, (name, seq))| {
            let contig = if name.is_empty() {
                format!("ctg{:06}", idx + 1)
            } else {
                name.clone()
            };
            let mut counts = ParentCounts::default();
            if k > 0 && seq.len() >= k {
                for window in seq.windows(k) {
                    counts.observe(&canonical_kmer(window), parent_kmers);
                }
            }
            ContigParentEvidence {
                contig,
                parent1_only_kmers: counts.parent1_only,
                parent2_only_kmers: counts.parent2_only,
                shared_parent_kmers: counts.shared,
                uninformative_kmers: counts.uninformative,
                assignment: counts.assignment(),
            }
        })
        .collect()
}

fn summarize(
    diploid_enabled: bool,
    parent_references_supplied: bool,
    unitigs: &[UnitigParentEvidence],
    contigs: &[ContigParentEvidence],
) -> DiploidEvidenceSummary {
    DiploidEvidenceSummary {
        diploid_enabled,
        parent_references_supplied,
        full_haplotype_fasta_claimed: false,
        parent1_only_nodes: unitigs.iter().map(|u| u.parent1_only_nodes).sum(),
        parent2_only_nodes: unitigs.iter().map(|u| u.parent2_only_nodes).sum(),
        shared_parent_nodes: unitigs.iter().map(|u| u.shared_parent_nodes).sum(),
        mixed_unitigs: unitigs
            .iter()
            .filter(|u| u.assignment == ParentAssignment::Mixed)
            .count(),
        parent_informative_unitigs: unitigs
            .iter()
            .filter(|u| {
                matches!(
                    u.assignment,
                    ParentAssignment::Parent1 | ParentAssignment::Parent2 | ParentAssignment::Mixed
                )
            })
            .count(),
        parent_informative_contigs: contigs
            .iter()
            .filter(|c| {
                matches!(
                    c.assignment,
                    ParentAssignment::Parent1 | ParentAssignment::Parent2 | ParentAssignment::Mixed
                )
            })
            .count(),
    }
}

#[derive(Default)]
struct ParentCounts {
    parent1_only: usize,
    parent2_only: usize,
    shared: usize,
    uninformative: usize,
}

impl ParentCounts {
    fn observe(&mut self, kmer: &[u8], parent_kmers: &ParentKmers) {
        if parent_kmers.parent1_only.contains(kmer) {
            self.parent1_only += 1;
        } else if parent_kmers.parent2_only.contains(kmer) {
            self.parent2_only += 1;
        } else if parent_kmers.shared.contains(kmer) {
            self.shared += 1;
        } else {
            self.uninformative += 1;
        }
    }

    fn assignment(&self) -> ParentAssignment {
        match (
            self.parent1_only > 0,
            self.parent2_only > 0,
            self.shared > 0,
        ) {
            (true, false, _) => ParentAssignment::Parent1,
            (false, true, _) => ParentAssignment::Parent2,
            (true, true, _) => ParentAssignment::Mixed,
            (false, false, true) => ParentAssignment::Shared,
            (false, false, false) => ParentAssignment::Uninformative,
        }
    }
}

pub fn parent_assignment_counts(report: &DiploidEvidenceReport) -> BTreeMap<String, usize> {
    let mut out = BTreeMap::new();
    for contig in &report.contigs {
        *out.entry(contig.assignment.as_gfa_tag_value().to_string())
            .or_insert(0) += 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::{
        build_diploid_evidence, gfa_path_tags, gfa_unitig_tags, ParentAssignment,
        ParentReferenceParams,
    };

    fn write_tmp(dir: &std::path::Path, name: &str, content: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, content).expect("write");
        path
    }

    #[test]
    fn classifies_parent_specific_unitigs_and_contigs_without_phasing_claim() {
        let dir = tempfile::tempdir().expect("tempdir");
        let p1 = write_tmp(dir.path(), "p1.fa", b">p1\nAAAAC\n");
        let p2 = write_tmp(dir.path(), "p2.fa", b">p2\nAAAAG\n");
        let refs = ParentReferenceParams {
            parent1: Some(p1),
            parent2: Some(p2),
        };
        let unitigs = vec![vec![b"AAA".to_vec(), b"AAC".to_vec()]];
        let contigs = vec![("".to_string(), b"AAAAC".to_vec())];

        let report = build_diploid_evidence(true, 3, &refs, &unitigs, &contigs).expect("report");

        assert!(report.summary.diploid_enabled);
        assert!(report.summary.parent_references_supplied);
        assert!(!report.summary.full_haplotype_fasta_claimed);
        assert_eq!(report.unitigs[0].assignment, ParentAssignment::Parent1);
        assert_eq!(report.contigs[0].assignment, ParentAssignment::Parent1);
        assert_eq!(
            gfa_path_tags(&report),
            vec![("ctg000001".to_string(), "PS:Z:parent1".to_string())]
        );
        assert_eq!(
            gfa_unitig_tags(&report),
            vec![(1, "PS:Z:parent1".to_string())]
        );
    }

    #[test]
    fn diploid_without_parent_refs_reports_metrics_but_no_parent_claim() {
        let refs = ParentReferenceParams::default();
        let report = build_diploid_evidence(true, 3, &refs, &[], &[]).expect("report");

        assert!(report.summary.diploid_enabled);
        assert!(!report.summary.parent_references_supplied);
        assert!(!report.summary.full_haplotype_fasta_claimed);
        assert!(report.parent_kmers.is_none());
        assert!(report.unitigs.is_empty());
        assert!(report.contigs.is_empty());
    }
}
