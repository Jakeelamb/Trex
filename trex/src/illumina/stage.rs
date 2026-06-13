//! Named stage runner for the Illumina assembler pipeline.

use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyStage {
    LoadReads,
    SelectK,
    WritePreprocessCheckpoint,
    CountKmers,
    TrustDiagnostics,
    WriteCountsCheckpoint,
    BuildSimplifiedGraph,
    StitchUnitigsAndContigs,
    PrimaryContigCollapse,
    PrimaryContigPathsForGfa,
    AnnotateGraph,
    DiagnoseFragmentation,
    BuildScaffoldArtifact,
    BuildDiploidEvidence,
    AuditAssembly,
    WriteExportCheckpoint,
    WriteOutputs,
}

impl AssemblyStage {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::LoadReads => "load_reads",
            Self::SelectK => "select_k",
            Self::WritePreprocessCheckpoint => "write_preprocess_checkpoint",
            Self::CountKmers => "count_kmers",
            Self::TrustDiagnostics => "build_trust_diagnostics",
            Self::WriteCountsCheckpoint => "write_counts_checkpoint",
            Self::BuildSimplifiedGraph => "build_simplified_graph",
            Self::StitchUnitigsAndContigs => "stitch_unitigs_and_contigs",
            Self::PrimaryContigCollapse => "phase2_primary_contig_collapse",
            Self::PrimaryContigPathsForGfa => "primary_contig_paths_for_gfa",
            Self::AnnotateGraph => "annotate_graph",
            Self::DiagnoseFragmentation => "diagnose_fragmentation",
            Self::BuildScaffoldArtifact => "build_scaffold_artifact",
            Self::BuildDiploidEvidence => "build_diploid_evidence",
            Self::AuditAssembly => "audit_assembly",
            Self::WriteExportCheckpoint => "write_export_checkpoint",
            Self::WriteOutputs => "write_outputs",
        }
    }
}

#[derive(Debug, Default)]
pub struct AssemblyStageRunner;

impl AssemblyStageRunner {
    pub fn new() -> Self {
        Self
    }

    pub fn run<T, E>(
        &self,
        stage: AssemblyStage,
        f: impl FnOnce() -> Result<T, E>,
    ) -> Result<T, E> {
        let start = Instant::now();
        let result = f();
        if result.is_ok() {
            tracing::info!(
                stage = stage.as_str(),
                elapsed_ms = start.elapsed().as_millis(),
                "illumina pipeline stage complete"
            );
        }
        result
    }

    pub fn observe(&self, stage: AssemblyStage, f: impl FnOnce()) {
        let start = Instant::now();
        f();
        tracing::info!(
            stage = stage.as_str(),
            elapsed_ms = start.elapsed().as_millis(),
            "illumina pipeline stage complete"
        );
    }

    pub fn log_stage(&self, stage: AssemblyStage, elapsed_ms: u128, message: &'static str) {
        tracing::info!(stage = stage.as_str(), elapsed_ms, message);
    }
}
