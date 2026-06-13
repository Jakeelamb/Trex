//! Output and sidecar writing for the Illumina assembler.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use serde::Serialize;

use crate::dbg::{
    unitig_adjacency_links, write_contigs_fasta, write_gfa1, write_unitigs_fasta, DbgGraph,
    GfaWriteOptions,
};
use crate::dbg::{GraphAnnotations, SimplifyDecisionLog};
use crate::error::{GraphError, TrexError};
use crate::evidence::EvidenceLedger;
use crate::illumina::audit::{audit_tsv, AssemblyAuditReport};
use crate::illumina::diploid::{gfa_path_tags, gfa_unitig_tags, DiploidEvidenceReport};
use crate::illumina::fragmentation::FragmentationReport;
use crate::illumina::multik::MultiKSelectionReport;
use crate::illumina::pipeline::AssembleOutputs;
use crate::illumina::scaffold::{scaffold_fasta_records, scaffold_gfa_paths, ScaffoldArtifact};
use crate::illumina::trust::TrustDiagnosticsReport;

pub struct AssemblyOutputBundle<'a> {
    pub graph: &'a DbgGraph,
    pub unitig_records: &'a [(String, Vec<u8>)],
    pub contig_records: &'a [(String, Vec<u8>)],
    pub unitig_paths: &'a [Vec<Vec<u8>>],
    pub primary_contig_paths_gfa: &'a [(String, Vec<(usize, char)>)],
    pub evidence: &'a EvidenceLedger,
    pub trust_report: &'a TrustDiagnosticsReport,
    pub graph_annotations: &'a GraphAnnotations,
    pub simplify_decisions: &'a SimplifyDecisionLog,
    pub scaffold_artifact: &'a ScaffoldArtifact,
    pub multi_k_selection: &'a MultiKSelectionReport,
    pub fragmentation_report: &'a FragmentationReport,
    pub audit_report: &'a AssemblyAuditReport,
    pub diploid_evidence: &'a DiploidEvidenceReport,
    pub diploid_enabled: bool,
}

pub struct AssemblyOutputWriter<'a> {
    outputs: &'a AssembleOutputs,
}

impl<'a> AssemblyOutputWriter<'a> {
    pub fn new(outputs: &'a AssembleOutputs) -> Self {
        Self { outputs }
    }

    pub fn write_all(&self, bundle: AssemblyOutputBundle<'_>) -> Result<(), TrexError> {
        self.prepare_directories()?;

        write_unitigs_fasta(&self.outputs.unitigs_path(), bundle.unitig_records)?;
        write_contigs_fasta(&self.outputs.contigs_path(), bundle.contig_records)?;

        let diploid_gfa_links = if bundle.diploid_enabled {
            Some(unitig_adjacency_links(bundle.graph, bundle.unitig_paths))
        } else {
            None
        };
        let parent_unitig_tags = gfa_unitig_tags(bundle.diploid_evidence);
        let parent_path_tags = gfa_path_tags(bundle.diploid_evidence);
        let scaffold_paths_gfa = scaffold_gfa_paths(bundle.scaffold_artifact);
        write_gfa1(
            &self.outputs.gfa_path_resolved(),
            bundle.unitig_records,
            GfaWriteOptions {
                phase2_illumina_diploid: bundle.diploid_enabled,
                diploid_unitig_links: diploid_gfa_links.as_deref(),
                primary_contig_paths: bundle.primary_contig_paths_gfa,
                scaffold_paths: &scaffold_paths_gfa,
                phase2_unphased_hap_paths: bundle.diploid_enabled,
                parent_unitig_tags: Some(parent_unitig_tags.as_slice()),
                parent_path_tags: Some(parent_path_tags.as_slice()),
            },
        )?;

        write_json_pretty(&self.outputs.evidence_path(), bundle.evidence)?;
        write_json_pretty(&self.outputs.trust_path(), bundle.trust_report)?;
        write_json_pretty(&self.outputs.annotations_path(), bundle.graph_annotations)?;
        write_json_pretty(
            &self.outputs.simplification_path(),
            bundle.simplify_decisions,
        )?;
        write_json_pretty(&self.outputs.scaffolds_path(), bundle.scaffold_artifact)?;
        write_scaffolds_fasta(
            &self.outputs.scaffolds_fasta_path(),
            &scaffold_fasta_records(bundle.scaffold_artifact, bundle.unitig_records),
        )?;
        write_multi_k_json(&self.outputs.multi_k_path(), bundle.multi_k_selection)?;
        write_json_pretty(
            &self.outputs.fragmentation_path(),
            bundle.fragmentation_report,
        )?;
        write_json_pretty(&self.outputs.audit_json_path(), bundle.audit_report)?;
        write_audit_tsv(&self.outputs.audit_tsv_path(), bundle.audit_report)?;
        write_json_pretty(&self.outputs.diploid_path(), bundle.diploid_evidence)?;
        Ok(())
    }

    fn prepare_directories(&self) -> Result<(), TrexError> {
        if self.outputs.out_dir.as_os_str() != "-" {
            std::fs::create_dir_all(&self.outputs.out_dir)
                .map_err(|e| TrexError::from(GraphError::from(e)))?;
        }

        for path in [
            self.outputs.unitigs_path(),
            self.outputs.contigs_path(),
            self.outputs.gfa_path_resolved(),
        ] {
            if path.as_os_str() == "-" {
                continue;
            }
            if let Some(parent) = path.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent)
                        .map_err(|e| TrexError::from(GraphError::from(e)))?;
                }
            }
        }
        Ok(())
    }
}

fn write_scaffolds_fasta(path: &Path, records: &[(String, Vec<u8>)]) -> Result<(), TrexError> {
    if records.is_empty() || path.as_os_str() == "-" {
        return Ok(());
    }
    write_contigs_fasta(path, records).map_err(TrexError::Graph)
}

fn write_multi_k_json(path: &Path, selection: &MultiKSelectionReport) -> Result<(), TrexError> {
    if !selection.enabled || path.as_os_str() == "-" {
        return Ok(());
    }
    write_json_pretty(path, selection)
}

fn write_audit_tsv(path: &Path, report: &AssemblyAuditReport) -> Result<(), TrexError> {
    if path.as_os_str() == "-" {
        return Ok(());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| TrexError::Graph(GraphError::Io(e)))?;
        }
    }
    std::fs::write(path, audit_tsv(report)).map_err(|e| TrexError::Graph(GraphError::Io(e)))?;
    Ok(())
}

fn write_json_pretty<T: Serialize>(path: &Path, value: &T) -> Result<(), TrexError> {
    if path.as_os_str() == "-" {
        return Ok(());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| TrexError::Graph(GraphError::Io(e)))?;
        }
    }
    let file = File::create(path).map_err(|e| TrexError::Graph(GraphError::Io(e)))?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, value).map_err(|e| TrexError::Graph(e.into()))?;
    writer
        .flush()
        .map_err(|e| TrexError::Graph(GraphError::Io(e)))?;
    Ok(())
}
