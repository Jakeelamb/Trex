//! Post-assembly audit reports. These flag suspicious output regions without repairing sequence.

use std::collections::{BTreeMap, HashMap};

use serde::{Deserialize, Serialize};

use crate::dbg::GraphAnnotations;
use crate::illumina::mate::MateBridgeStats;
use crate::kmer::canonical_kmer;

pub const ASSEMBLY_AUDIT_SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AuditFindingKind {
    LowReadSupport,
    AbnormalPairHint,
    CollapsedRepeatSuspicion,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuditFinding {
    pub kind: AuditFindingKind,
    pub severity: String,
    pub contig: Option<String>,
    pub start: Option<usize>,
    pub end: Option<usize>,
    pub message: String,
    pub counters: BTreeMap<String, u64>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AssemblyAuditSummary {
    pub k: usize,
    pub contigs: usize,
    pub contig_bases: usize,
    pub assembly_kmers: usize,
    pub trusted_supported_kmers: usize,
    pub low_support_kmers: usize,
    pub low_support_regions: usize,
    pub abnormal_pair_hints: usize,
    pub collapsed_repeat_suspicions: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AssemblyAuditReport {
    pub schema_version: u64,
    pub summary: AssemblyAuditSummary,
    pub findings: Vec<AuditFinding>,
}

pub fn audit_assembly(
    k: usize,
    contigs: &[(String, Vec<u8>)],
    trusted: &[(Vec<u8>, u64)],
    annotations: &GraphAnnotations,
    mate_stats: Option<&MateBridgeStats>,
) -> AssemblyAuditReport {
    let trusted_map: HashMap<&[u8], u64> = trusted
        .iter()
        .map(|(kmer, multiplicity)| (kmer.as_slice(), *multiplicity))
        .collect();
    let mut findings = Vec::new();
    let mut assembly_kmers = 0usize;
    let mut trusted_supported_kmers = 0usize;
    let mut low_support_kmers = 0usize;
    let mut low_support_regions = 0usize;
    let mut contig_bases = 0usize;

    for (idx, (name, seq)) in contigs.iter().enumerate() {
        contig_bases += seq.len();
        let contig_name = if name.is_empty() {
            format!("ctg{:06}", idx + 1)
        } else {
            name.clone()
        };
        let mut region_start: Option<usize> = None;
        let mut region_kmers = 0usize;

        if k > 0 && seq.len() >= k {
            for (pos, window) in seq.windows(k).enumerate() {
                assembly_kmers += 1;
                let canonical = canonical_kmer(window);
                if trusted_map.contains_key(canonical.as_slice()) {
                    trusted_supported_kmers += 1;
                    if let Some(start) = region_start.take() {
                        low_support_regions += 1;
                        findings.push(low_support_region_finding(
                            &contig_name,
                            start,
                            pos + k - 1,
                            region_kmers,
                        ));
                        region_kmers = 0;
                    }
                } else {
                    low_support_kmers += 1;
                    if region_start.is_none() {
                        region_start = Some(pos);
                    }
                    region_kmers += 1;
                }
            }
        }

        if let Some(start) = region_start {
            low_support_regions += 1;
            findings.push(low_support_region_finding(
                &contig_name,
                start,
                seq.len(),
                region_kmers,
            ));
        }
    }

    let abnormal_pair_hints = mate_stats
        .map(|stats| {
            stats
                .trusted_endpoint_pairs
                .saturating_sub(stats.existing_edge_pairs)
        })
        .unwrap_or(0);
    if abnormal_pair_hints > 0 {
        findings.push(abnormal_pair_finding(abnormal_pair_hints));
    }

    let collapsed_repeat_suspicions = usize::from(
        annotations.summary.repeat_suspected_nodes > 0
            || annotations.summary.repeat_suspected_unitigs > 0,
    );
    if collapsed_repeat_suspicions > 0 {
        findings.push(repeat_suspicion_finding(annotations));
    }

    AssemblyAuditReport {
        schema_version: ASSEMBLY_AUDIT_SCHEMA_VERSION,
        summary: AssemblyAuditSummary {
            k,
            contigs: contigs.len(),
            contig_bases,
            assembly_kmers,
            trusted_supported_kmers,
            low_support_kmers,
            low_support_regions,
            abnormal_pair_hints,
            collapsed_repeat_suspicions,
        },
        findings,
    }
}

pub fn audit_tsv(report: &AssemblyAuditReport) -> String {
    let mut out = String::from("kind\tseverity\tcontig\tstart\tend\tmessage\n");
    for finding in &report.findings {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            finding_kind_str(&finding.kind),
            finding.severity,
            finding.contig.as_deref().unwrap_or("."),
            finding
                .start
                .map(|v| v.to_string())
                .unwrap_or_else(|| ".".to_string()),
            finding
                .end
                .map(|v| v.to_string())
                .unwrap_or_else(|| ".".to_string()),
            finding.message.replace('\t', " ")
        ));
    }
    out
}

fn low_support_region_finding(
    contig_name: &str,
    start: usize,
    end: usize,
    kmer_count: usize,
) -> AuditFinding {
    let mut counters = BTreeMap::new();
    counters.insert("low_support_kmers".to_string(), kmer_count as u64);
    AuditFinding {
        kind: AuditFindingKind::LowReadSupport,
        severity: "warning".to_string(),
        contig: Some(contig_name.to_string()),
        start: Some(start),
        end: Some(end),
        message: "assembly k-mers in this interval are absent from trusted read k-mers".to_string(),
        counters,
    }
}

fn abnormal_pair_finding(abnormal_pair_hints: usize) -> AuditFinding {
    let mut counters = BTreeMap::new();
    counters.insert(
        "abnormal_pair_hints".to_string(),
        abnormal_pair_hints as u64,
    );
    AuditFinding {
        kind: AuditFindingKind::AbnormalPairHint,
        severity: "info".to_string(),
        contig: None,
        start: None,
        end: None,
        message: "trusted mate endpoint pairs were observed without existing DBG edge support"
            .to_string(),
        counters,
    }
}

fn repeat_suspicion_finding(annotations: &GraphAnnotations) -> AuditFinding {
    let mut counters = BTreeMap::new();
    counters.insert(
        "repeat_suspected_nodes".to_string(),
        annotations.summary.repeat_suspected_nodes as u64,
    );
    counters.insert(
        "repeat_suspected_unitigs".to_string(),
        annotations.summary.repeat_suspected_unitigs as u64,
    );
    AuditFinding {
        kind: AuditFindingKind::CollapsedRepeatSuspicion,
        severity: "info".to_string(),
        contig: None,
        start: None,
        end: None,
        message: "repeat-like graph annotations remain in the emitted assembly graph".to_string(),
        counters,
    }
}

fn finding_kind_str(kind: &AuditFindingKind) -> &'static str {
    match kind {
        AuditFindingKind::LowReadSupport => "low_read_support",
        AuditFindingKind::AbnormalPairHint => "abnormal_pair_hint",
        AuditFindingKind::CollapsedRepeatSuspicion => "collapsed_repeat_suspicion",
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{audit_assembly, audit_tsv, AuditFindingKind};
    use crate::dbg::{
        GraphAnnotationSummary, GraphAnnotations, NodeAnnotation, NodeDepthClass, UnitigAnnotation,
    };

    fn empty_annotations() -> GraphAnnotations {
        GraphAnnotations {
            nodes: BTreeMap::new(),
            summary: GraphAnnotationSummary {
                baseline_multiplicity: 1,
                node_count: 0,
                low_copy_nodes: 0,
                single_copy_nodes: 0,
                high_copy_nodes: 0,
                repeat_suspected_nodes: 0,
                unitig_count: 0,
                repeat_suspected_unitigs: 0,
                unitigs: Vec::new(),
            },
        }
    }

    #[test]
    fn flags_contig_kmers_absent_from_trusted_reads_without_repairing_sequence() {
        let contigs = vec![("ctg1".to_string(), b"AAAAC".to_vec())];
        let trusted = vec![(b"AAA".to_vec(), 2)];
        let report = audit_assembly(3, &contigs, &trusted, &empty_annotations(), None);

        assert_eq!(report.summary.assembly_kmers, 3);
        assert_eq!(report.summary.trusted_supported_kmers, 2);
        assert_eq!(report.summary.low_support_kmers, 1);
        assert_eq!(report.summary.low_support_regions, 1);
        assert_eq!(report.findings[0].kind, AuditFindingKind::LowReadSupport);
        assert_eq!(report.findings[0].contig.as_deref(), Some("ctg1"));
        assert_eq!(report.findings[0].start, Some(2));
        assert_eq!(report.findings[0].end, Some(5));
    }

    #[test]
    fn reports_repeat_suspicion_from_graph_annotations() {
        let mut annotations = empty_annotations();
        annotations.nodes.insert(
            "AAA".to_string(),
            NodeAnnotation {
                multiplicity: 10,
                degree: 3,
                depth_class: NodeDepthClass::HighCopy,
                repeat_suspected: true,
            },
        );
        annotations.summary.repeat_suspected_nodes = 1;
        annotations.summary.repeat_suspected_unitigs = 1;
        annotations.summary.unitigs.push(UnitigAnnotation {
            unitig_index: 0,
            node_count: 1,
            min_multiplicity: 10,
            max_multiplicity: 10,
            mean_multiplicity: 10.0,
            repeat_suspected: true,
        });

        let report = audit_assembly(3, &[], &[], &annotations, None);
        assert_eq!(report.summary.collapsed_repeat_suspicions, 1);
        assert!(report
            .findings
            .iter()
            .any(|finding| finding.kind == AuditFindingKind::CollapsedRepeatSuspicion));
        assert!(audit_tsv(&report).contains("collapsed_repeat_suspicion"));
    }
}
