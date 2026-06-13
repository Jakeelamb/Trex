//! Named stage runner for the Illumina assembler pipeline.

use std::cell::RefCell;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
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
pub struct AssemblyStageRunner {
    reports: RefCell<Vec<AssemblyStageReport>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AssemblyStageOutcome {
    Complete,
    Failed,
}

impl AssemblyStageOutcome {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Failed => "failed",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct AssemblyStageReport {
    pub stage: AssemblyStage,
    pub outcome: AssemblyStageOutcome,
    pub elapsed_ms: u128,
}

impl AssemblyStageRunner {
    pub fn new() -> Self {
        Self {
            reports: RefCell::new(Vec::new()),
        }
    }

    pub fn run<T, E>(
        &self,
        stage: AssemblyStage,
        f: impl FnOnce() -> Result<T, E>,
    ) -> Result<T, E> {
        let start = Instant::now();
        let result = f();
        let elapsed_ms = start.elapsed().as_millis();
        let outcome = if result.is_ok() {
            AssemblyStageOutcome::Complete
        } else {
            AssemblyStageOutcome::Failed
        };
        self.record(stage, outcome, elapsed_ms);
        tracing::info!(
            stage = stage.as_str(),
            outcome = outcome.as_str(),
            elapsed_ms,
            "illumina pipeline stage finished"
        );
        result
    }

    pub fn observe(&self, stage: AssemblyStage, f: impl FnOnce()) {
        let start = Instant::now();
        f();
        let elapsed_ms = start.elapsed().as_millis();
        self.record(stage, AssemblyStageOutcome::Complete, elapsed_ms);
        tracing::info!(
            stage = stage.as_str(),
            outcome = AssemblyStageOutcome::Complete.as_str(),
            elapsed_ms,
            "illumina pipeline stage finished"
        );
    }

    pub fn log_stage(&self, stage: AssemblyStage, elapsed_ms: u128, message: &'static str) {
        tracing::info!(stage = stage.as_str(), elapsed_ms, message);
    }

    pub fn reports(&self) -> Vec<AssemblyStageReport> {
        self.reports.borrow().clone()
    }

    fn record(&self, stage: AssemblyStage, outcome: AssemblyStageOutcome, elapsed_ms: u128) {
        self.reports.borrow_mut().push(AssemblyStageReport {
            stage,
            outcome,
            elapsed_ms,
        });
    }
}
