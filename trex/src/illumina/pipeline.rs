//! End-to-end **Illumina** preprocess → *k*-mer counts → trusted DBG → simplify → unitigs → contigs → exports.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use tracing::instrument;

use crate::dbg::{
    annotate_graph, assert_no_self_loops, build_dbg, extract_unitigs, forward_representatives,
    primary_contig_paths_for_gfa, reference_contig_paths, run_simplification_schedule,
    stitch_sequence, ContigWalkTieBreak, DbgGraph, DiploidSimplifyMode, GraphAnnotations,
    SimplifyDecisionLog, SimplifyParams, SimplifyStats,
};
use crate::error::{GraphError, IngestError, TrexError};
use crate::evidence::EvidenceLedger;
use crate::illumina::audit::{audit_assembly, AssemblyAuditReport};
use crate::illumina::checkpoint::{CheckpointStore, GraphCheckpointIdentity, GraphStageArtifacts};
use crate::illumina::counts::enumerate_sorted_counts;
use crate::illumina::diploid::{
    build_diploid_evidence, DiploidEvidenceReport, ParentReferenceParams,
};
use crate::illumina::fragmentation::{diagnose_fragmentation, FragmentationReport};
use crate::illumina::io::{read_fasta_records_maybe_gzip, read_fastq_records_maybe_gzip};
use crate::illumina::mate;
use crate::illumina::multik::{select_k, MultiKParams, MultiKSelectionReport};
use crate::illumina::output::{AssemblyOutputBundle, AssemblyOutputWriter};
use crate::illumina::paired::validate_pair_parity;
use crate::illumina::phase2_primary;
use crate::illumina::read::Read;
use crate::illumina::scaffold::{build_scaffold_artifact, ScaffoldArtifact};
use crate::illumina::stage::{AssemblyStage, AssemblyStageReport, AssemblyStageRunner};
use crate::illumina::trust::{build_trust_diagnostics, TrustDiagnosticsReport};
use crate::kmer::apply_trusted_threshold;

type SequenceRecords = Vec<(String, Vec<u8>)>;
type ContigPaths = Vec<Vec<Vec<u8>>>;
type StitchedAssembly = (SequenceRecords, SequenceRecords, ContigPaths, ContigPaths);

/// Default output layout (**Phase-1 export layout**): separate files under `out_dir`.
#[derive(Debug, Clone)]
pub struct AssembleOutputs {
    pub out_dir: PathBuf,
    pub unitigs_fasta: PathBuf,
    pub contigs_fasta: PathBuf,
    pub gfa_path: PathBuf,
}

impl Default for AssembleOutputs {
    fn default() -> Self {
        Self {
            out_dir: PathBuf::from("."),
            unitigs_fasta: PathBuf::from("unitigs.fa"),
            contigs_fasta: PathBuf::from("contigs.fa"),
            gfa_path: PathBuf::from("graph.gfa"),
        }
    }
}

impl AssembleOutputs {
    fn resolve(&self, rel: &Path) -> PathBuf {
        if rel.as_os_str() == "-" {
            return PathBuf::from("-");
        }
        if rel.is_absolute() {
            rel.to_path_buf()
        } else {
            self.out_dir.join(rel)
        }
    }

    pub fn unitigs_path(&self) -> PathBuf {
        self.resolve(&self.unitigs_fasta)
    }

    pub fn contigs_path(&self) -> PathBuf {
        self.resolve(&self.contigs_fasta)
    }

    pub fn gfa_path_resolved(&self) -> PathBuf {
        self.resolve(&self.gfa_path)
    }

    pub fn evidence_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("evidence.json"))
        }
    }

    pub fn trust_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("trust.json"))
        }
    }

    pub fn annotations_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("annotations.json"))
        }
    }

    pub fn simplification_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("simplification.json"))
        }
    }

    pub fn scaffolds_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("scaffolds.json"))
        }
    }

    pub fn scaffolds_fasta_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("scaffolds.fa"))
        }
    }

    pub fn multi_k_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("multi_k.json"))
        }
    }

    pub fn fragmentation_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("fragmentation.json"))
        }
    }

    pub fn audit_json_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("audit.json"))
        }
    }

    pub fn audit_tsv_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("audit.tsv"))
        }
    }

    pub fn diploid_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("diploid.json"))
        }
    }

    pub fn stages_path(&self) -> PathBuf {
        if self.out_dir.as_os_str() == "-" {
            PathBuf::from("-")
        } else {
            self.resolve(Path::new("stages.json"))
        }
    }
}

/// Optional overrides for **Phase-1 graph simplification** (tips + bounded diamond bubbles +
/// low-copy component pruning).
#[derive(Debug, Clone, Default)]
pub struct SimplifyOverrides {
    pub max_tip_bases: Option<usize>,
    pub tip_max_multiplicity: Option<u64>,
    pub max_bubble_vertices: Option<usize>,
    pub max_bubble_internal_bases: Option<usize>,
    pub max_low_coverage_component_bases: Option<usize>,
    pub low_coverage_component_max_multiplicity: Option<u64>,
}

fn graph_checkpoint_identity(params: &AssembleParams) -> GraphCheckpointIdentity {
    GraphCheckpointIdentity {
        diploid_enabled: params.diploid.enabled,
        diploid_paired_end: params.r2_path.is_some(),
        diploid_insert_mean_bp: params.diploid.insert_mean_bp,
        diploid_insert_stddev_bp: params.diploid.insert_stddev_bp,
        phase2_mate_bridge_v1: params.diploid.enabled
            && params.r2_path.is_some()
            && params.diploid.insert_mean_bp.is_some(),
    }
}

fn merge_simplify_params(k: usize, overrides: &SimplifyOverrides) -> SimplifyParams {
    let mut p = SimplifyParams::for_k(k);
    if let Some(v) = overrides.max_tip_bases {
        p.max_tip_bases = v;
    }
    if let Some(v) = overrides.tip_max_multiplicity {
        p.tip_max_multiplicity = v;
    }
    if let Some(v) = overrides.max_bubble_vertices {
        p.max_bubble_vertices = v;
    }
    if let Some(v) = overrides.max_bubble_internal_bases {
        p.max_bubble_internal_bases = v;
    }
    if let Some(v) = overrides.max_low_coverage_component_bases {
        p.max_low_coverage_component_bases = v;
    }
    if let Some(v) = overrides.low_coverage_component_max_multiplicity {
        p.low_coverage_component_max_multiplicity = v;
    }
    p
}

fn apply_phase2_mate_bridge(
    graph: &mut DbgGraph,
    reads: &[Read],
    paired_r1_len: Option<usize>,
    k: usize,
    insert_mean_bp: Option<u64>,
    insert_stddev_bp: Option<u64>,
    enabled: bool,
) -> Option<mate::MateBridgeStats> {
    if !enabled {
        return None;
    }
    let Some(n) = paired_r1_len else {
        tracing::warn!(
            "Phase-2 mate bridge skipped: missing preprocess/pair_layout.json; re-run without --resume to record paired layout"
        );
        return None;
    };

    let stats = mate::boost_mate_pairs_on_existing_dbg_edges(
        graph,
        reads,
        n,
        k,
        insert_mean_bp,
        insert_stddev_bp,
    );
    tracing::info!(
        pairs_seen = stats.pairs_seen,
        pairs_with_endpoint_kmers = stats.pairs_with_endpoint_kmers,
        trusted_endpoint_pairs = stats.trusted_endpoint_pairs,
        existing_edge_pairs = stats.existing_edge_pairs,
        report_only_pairs = stats.report_only_pairs,
        boosted_edges = stats.boosted_edges,
        "Phase-2 Illumina mate-pair bridge evidence"
    );
    Some(stats)
}

fn stitch_unitigs_and_contigs(
    graph: &crate::dbg::DbgGraph,
    forward: &HashMap<Vec<u8>, Vec<u8>>,
    k: usize,
    tie_break: ContigWalkTieBreak,
) -> Result<StitchedAssembly, TrexError> {
    let raw_unitig_paths = extract_unitigs(graph);
    let mut unitig_paths: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut unitig_records: Vec<(String, Vec<u8>)> = Vec::new();
    for p in raw_unitig_paths {
        match stitch_sequence(&p, forward, k) {
            Ok(seq) => {
                unitig_paths.push(p);
                unitig_records.push((String::new(), seq));
            }
            Err(GraphError::OrientationConflict) => {}
            Err(e) => return Err(TrexError::Graph(e)),
        }
    }
    let contig_paths = reference_contig_paths(graph, forward, k, tie_break)?;
    let mut contig_records: Vec<(String, Vec<u8>)> = Vec::new();
    for p in &contig_paths {
        let seq = stitch_sequence(p, forward, k)?;
        contig_records.push((String::new(), seq));
    }
    Ok((unitig_records, contig_records, unitig_paths, contig_paths))
}

struct BuiltGraph {
    graph: DbgGraph,
    simplify_stats: SimplifyStats,
    simplify_decisions: SimplifyDecisionLog,
    mate_bridge_stats: Option<mate::MateBridgeStats>,
}

impl BuiltGraph {
    fn stage_artifacts(&self) -> GraphStageArtifacts {
        GraphStageArtifacts {
            simplify_stats: self.simplify_stats,
            simplify_decisions: self.simplify_decisions.clone(),
            mate_bridge_stats: self.mate_bridge_stats.clone(),
        }
    }
}

struct BuildGraphInputs<'a> {
    reads: &'a [Read],
    trusted: &'a [(Vec<u8>, u64)],
    paired_r1_len: Option<usize>,
    selected_k: usize,
    simplify_params: &'a SimplifyParams,
    diploid_simplify: Option<DiploidSimplifyMode>,
    graph_identity: GraphCheckpointIdentity,
    diploid: &'a DiploidParams,
}

fn build_simplified_graph(inputs: BuildGraphInputs<'_>) -> Result<BuiltGraph, TrexError> {
    let mut graph = build_dbg(inputs.reads, inputs.selected_k, inputs.trusted)?;
    let mate_bridge_stats = apply_phase2_mate_bridge(
        &mut graph,
        inputs.reads,
        inputs.paired_r1_len,
        inputs.selected_k,
        inputs.diploid.insert_mean_bp,
        inputs.diploid.insert_stddev_bp,
        inputs.graph_identity.phase2_mate_bridge_v1,
    );
    let (simplify_stats, simplify_decisions) =
        run_simplification_schedule(&mut graph, inputs.simplify_params, inputs.diploid_simplify);
    Ok(BuiltGraph {
        graph,
        simplify_stats,
        simplify_decisions,
        mate_bridge_stats,
    })
}

/// **Phase-2 Illumina diploid** options (experimental). When `enabled`, diamond bubble simplification
/// retains near-balanced branches and GFA exports carry a `trex-phase2-illumina` header tag.
#[derive(Debug, Clone, Default)]
pub struct DiploidParams {
    pub enabled: bool,
    /// Reserved for mate-distance–aware simplification; stored in graph checkpoints for resume identity.
    pub insert_mean_bp: Option<u64>,
    pub insert_stddev_bp: Option<u64>,
    /// Optional parent references used only for evidence/reporting; they do not change graph shape.
    pub parent_references: ParentReferenceParams,
}

/// Parameters for the current **Phase-1** `assemble` slice (ingest → counts → DBG → exports).
#[derive(Debug, Clone)]
pub struct AssembleParams {
    pub r1_path: PathBuf,
    pub r2_path: Option<PathBuf>,
    pub k: usize,
    /// Global trusted multiplicity floor *T* (**Phase-1 trusted k-mer rule**).
    pub trusted_threshold: u64,
    pub checkpoint_root: Option<PathBuf>,
    pub resume: bool,
    pub strict_checkpoints: bool,
    pub simplify: SimplifyOverrides,
    pub diploid: DiploidParams,
    pub multi_k: MultiKParams,
    pub outputs: AssembleOutputs,
}

/// Outputs after a full **assemble** run.
#[derive(Debug, Clone)]
pub struct AssembleResult {
    pub reads: Vec<Read>,
    pub trusted_kmers: Vec<(Vec<u8>, u64)>,
    pub total_unique_kmers: usize,
    pub trust_report: TrustDiagnosticsReport,
    pub multi_k_selection: MultiKSelectionReport,
    pub simplify_stats: SimplifyStats,
    pub simplify_decisions: SimplifyDecisionLog,
    pub evidence: EvidenceLedger,
    pub graph_annotations: GraphAnnotations,
    pub scaffold_artifact: ScaffoldArtifact,
    pub fragmentation_report: FragmentationReport,
    pub audit_report: AssemblyAuditReport,
    pub diploid_evidence: DiploidEvidenceReport,
    pub mate_bridge_stats: Option<mate::MateBridgeStats>,
    pub stage_reports: Vec<AssemblyStageReport>,
    pub unitig_count: usize,
    pub contig_count: usize,
    pub outputs: AssembleOutputs,
}

fn input_is_fasta(path: &Path) -> bool {
    let s = path.to_string_lossy().to_ascii_lowercase();
    s.ends_with(".fa")
        || s.ends_with(".fasta")
        || s.ends_with(".fna")
        || s.ends_with(".fa.gz")
        || s.ends_with(".fasta.gz")
        || s.ends_with(".fna.gz")
}

#[instrument(
    skip(params),
    fields(
        k = params.k,
        T = params.trusted_threshold,
        diploid = params.diploid.enabled,
        paired = params.r2_path.is_some(),
    )
)]
pub fn assemble_illumina(params: &AssembleParams) -> Result<AssembleResult, TrexError> {
    let stages = AssemblyStageRunner::new();
    let ck = params
        .checkpoint_root
        .as_ref()
        .map(|p| CheckpointStore::new(p.clone(), params.strict_checkpoints));

    let (reads, paired_r1_len) = stages.run(AssemblyStage::LoadReads, || {
        if params.resume {
            match &ck {
                Some(c) => match c.load_preprocess()? {
                    Some(r) => {
                        let pl = c.load_pair_layout()?;
                        Ok::<_, TrexError>((r, pl))
                    }
                    None => Ok(load_reads(params)?),
                },
                None => Err(TrexError::from(IngestError::ResumeRequiresCheckpointRoot)),
            }
        } else {
            Ok(load_reads(params)?)
        }
    })?;

    if let Some(s) = Read::shortest_length(&reads) {
        if !params.multi_k.enabled() && s < params.k {
            return Err(IngestError::KTooLarge {
                k: params.k,
                shortest: s,
            }
            .into());
        }
    } else {
        return Err(IngestError::FastqFormat("no reads ingested".into()).into());
    }

    let multi_k_selection = stages.run(AssemblyStage::SelectK, || {
        select_k(&reads, params.k, params.trusted_threshold, &params.multi_k)
    })?;
    let selected_k = multi_k_selection.selected_k;
    if multi_k_selection.enabled {
        tracing::info!(
            requested_k = multi_k_selection.requested_k,
            selected_k = selected_k,
            candidates = multi_k_selection.candidates.len(),
            "multi-k graph selection"
        );
    }
    let selected_ck = ck.as_ref().map(|c| {
        if multi_k_selection.enabled {
            c.selected_k(selected_k)
        } else {
            c.clone()
        }
    });

    if let Some(ref c) = ck {
        stages.run(AssemblyStage::WritePreprocessCheckpoint, || {
            c.write_preprocess(&reads, paired_r1_len)
        })?;
    }

    let merged = stages.run(AssemblyStage::CountKmers, || {
        if params.resume {
            match &selected_ck {
                Some(c) => match c.load_counts(selected_k)? {
                    Some(rows) => Ok::<_, TrexError>(rows),
                    None => Ok(enumerate_sorted_counts(&reads, selected_k)?),
                },
                None => Ok(enumerate_sorted_counts(&reads, selected_k)?),
            }
        } else {
            Ok(enumerate_sorted_counts(&reads, selected_k)?)
        }
    })?;

    let total_unique = merged.len();
    let trust_report = stages.run(AssemblyStage::TrustDiagnostics, || {
        Ok::<_, TrexError>(build_trust_diagnostics(
            selected_k,
            params.trusted_threshold,
            &merged,
        ))
    })?;

    if let Some(ref c) = selected_ck {
        stages.run(AssemblyStage::WriteCountsCheckpoint, || {
            c.write_counts(selected_k, &merged)
        })?;
    }

    let trusted = apply_trusted_threshold(merged, params.trusted_threshold);

    let forward = forward_representatives(&reads, selected_k)?;
    let simplify_params = merge_simplify_params(selected_k, &params.simplify);
    let diploid_simplify = params.diploid.enabled.then_some(DiploidSimplifyMode);
    let graph_ck_id = graph_checkpoint_identity(params);
    let mut simplify_stats = SimplifyStats::default();
    let mut simplify_decisions = SimplifyDecisionLog::default();
    let mut mate_bridge_stats = None;
    let mut evidence = EvidenceLedger::new();

    let graph = stages.run(AssemblyStage::BuildSimplifiedGraph, || {
        if let Some(ref c) = selected_ck {
            if params.resume {
                match c.load_graph(selected_k, &graph_ck_id)? {
                    Some(g) => match c.load_graph_stage_artifacts(selected_k)? {
                        Some(artifacts) => {
                            simplify_stats = artifacts.simplify_stats;
                            simplify_decisions = artifacts.simplify_decisions;
                            mate_bridge_stats = artifacts.mate_bridge_stats;
                            Ok::<_, TrexError>(g)
                        }
                        None => {
                            tracing::warn!(
                                "graph checkpoint missing stage artifacts; rebuilding graph to avoid stale sidecars"
                            );
                            let built = build_simplified_graph(BuildGraphInputs {
                                reads: &reads,
                                trusted: &trusted,
                                paired_r1_len,
                                selected_k,
                                simplify_params: &simplify_params,
                                diploid_simplify,
                                graph_identity: graph_ck_id,
                                diploid: &params.diploid,
                            })?;
                            simplify_stats = built.simplify_stats;
                            simplify_decisions = built.simplify_decisions.clone();
                            mate_bridge_stats = built.mate_bridge_stats.clone();
                            c.write_graph(&built.graph, &graph_ck_id)?;
                            c.write_graph_stage_artifacts(selected_k, &built.stage_artifacts())?;
                            Ok(built.graph)
                        }
                    },
                    None => {
                        let built = build_simplified_graph(BuildGraphInputs {
                            reads: &reads,
                            trusted: &trusted,
                            paired_r1_len,
                            selected_k,
                            simplify_params: &simplify_params,
                            diploid_simplify,
                            graph_identity: graph_ck_id,
                            diploid: &params.diploid,
                        })?;
                        simplify_stats = built.simplify_stats;
                        simplify_decisions = built.simplify_decisions.clone();
                        mate_bridge_stats = built.mate_bridge_stats.clone();
                        c.write_graph(&built.graph, &graph_ck_id)?;
                        c.write_graph_stage_artifacts(selected_k, &built.stage_artifacts())?;
                        Ok(built.graph)
                    }
                }
            } else {
                let built = build_simplified_graph(BuildGraphInputs {
                    reads: &reads,
                    trusted: &trusted,
                    paired_r1_len,
                    selected_k,
                    simplify_params: &simplify_params,
                    diploid_simplify,
                    graph_identity: graph_ck_id,
                    diploid: &params.diploid,
                })?;
                simplify_stats = built.simplify_stats;
                simplify_decisions = built.simplify_decisions.clone();
                mate_bridge_stats = built.mate_bridge_stats.clone();
                c.write_graph(&built.graph, &graph_ck_id)?;
                c.write_graph_stage_artifacts(selected_k, &built.stage_artifacts())?;
                Ok(built.graph)
            }
        } else {
            let built = build_simplified_graph(BuildGraphInputs {
                reads: &reads,
                trusted: &trusted,
                paired_r1_len,
                selected_k,
                simplify_params: &simplify_params,
                diploid_simplify,
                graph_identity: graph_ck_id,
                diploid: &params.diploid,
            })?;
            simplify_stats = built.simplify_stats;
            simplify_decisions = built.simplify_decisions;
            mate_bridge_stats = built.mate_bridge_stats;
            Ok(built.graph)
        }
    })?;

    if let Some(stats) = mate_bridge_stats.as_ref() {
        evidence.push(stats.evidence_record());
    }

    assert_no_self_loops(&graph)?;
    tracing::info!(
        tips_removed = simplify_stats.tips_removed,
        diamond_bubbles_resolved = simplify_stats.diamond_bubbles_resolved,
        low_coverage_components_removed = simplify_stats.low_coverage_components_removed,
        diploid_diamonds_retained = simplify_stats.diploid_diamonds_retained,
        repeat_guarded_diamonds_retained = simplify_stats.repeat_guarded_diamonds_retained,
        ambiguous_k22_diamonds_skipped = simplify_stats.ambiguous_k22_diamonds_skipped,
        "graph simplification decisions"
    );

    if params.diploid.enabled {
        tracing::info!(
            diploid_enabled = true,
            paired_end = params.r2_path.is_some(),
            insert_mean_bp = ?params.diploid.insert_mean_bp,
            insert_stddev_bp = ?params.diploid.insert_stddev_bp,
            "Phase-2 Illumina diploid profile (experimental)"
        );
    }

    let contig_tie_break = if params.diploid.enabled {
        ContigWalkTieBreak::Phase2DiploidNodeMul
    } else {
        ContigWalkTieBreak::Phase1Lex
    };

    let (unitig_records, mut contig_records, unitig_paths, contig_paths) =
        stages.run(AssemblyStage::StitchUnitigsAndContigs, || {
            match (&selected_ck, params.resume) {
                (Some(c), true) => match c.load_export(selected_k)? {
                    Some((ur, cr)) => {
                        let up = extract_unitigs(&graph);
                        let cp =
                            reference_contig_paths(&graph, &forward, selected_k, contig_tie_break)?;
                        Ok::<_, TrexError>((ur, cr, up, cp))
                    }
                    None => Ok(stitch_unitigs_and_contigs(
                        &graph,
                        &forward,
                        selected_k,
                        contig_tie_break,
                    )?),
                },
                _ => Ok(stitch_unitigs_and_contigs(
                    &graph,
                    &forward,
                    selected_k,
                    contig_tie_break,
                )?),
            }
        })?;

    if params.diploid.enabled {
        stages.observe(AssemblyStage::PrimaryContigCollapse, || {
            let trusted_map: HashMap<Vec<u8>, u64> = trusted.iter().cloned().collect();
            for (_hdr, seq) in contig_records.iter_mut() {
                phase2_primary::collapse_primary_contig_by_trusted_kmers(
                    seq,
                    selected_k,
                    &trusted_map,
                );
            }
        });
    }

    let primary_contig_paths_gfa = stages.run(AssemblyStage::PrimaryContigPathsForGfa, || {
        Ok::<_, TrexError>(primary_contig_paths_for_gfa(&contig_paths, &unitig_paths))
    })?;
    let graph_annotations = stages.run(AssemblyStage::AnnotateGraph, || {
        Ok::<_, TrexError>(annotate_graph(&graph, &unitig_paths))
    })?;
    let fragmentation_report = stages.run(AssemblyStage::DiagnoseFragmentation, || {
        Ok::<_, TrexError>(diagnose_fragmentation(
            &graph,
            &contig_paths,
            &graph_annotations,
        ))
    })?;
    let scaffold_artifact = stages.run(AssemblyStage::BuildScaffoldArtifact, || {
        Ok::<_, TrexError>(build_scaffold_artifact(
            mate_bridge_stats
                .as_ref()
                .map(|stats| stats.candidates.as_slice())
                .unwrap_or(&[]),
            &unitig_paths,
            selected_k,
            &fragmentation_report,
        ))
    })?;
    let diploid_evidence = stages.run(AssemblyStage::BuildDiploidEvidence, || {
        build_diploid_evidence(
            params.diploid.enabled,
            selected_k,
            &params.diploid.parent_references,
            &unitig_paths,
            &contig_records,
        )
    })?;
    let audit_report = stages.run(AssemblyStage::AuditAssembly, || {
        Ok::<_, TrexError>(audit_assembly(
            selected_k,
            &contig_records,
            &trusted,
            &graph_annotations,
            mate_bridge_stats.as_ref(),
        ))
    })?;
    tracing::info!(
        baseline_multiplicity = graph_annotations.summary.baseline_multiplicity,
        high_copy_nodes = graph_annotations.summary.high_copy_nodes,
        repeat_suspected_nodes = graph_annotations.summary.repeat_suspected_nodes,
        repeat_suspected_unitigs = graph_annotations.summary.repeat_suspected_unitigs,
        "graph copy-number and repeat annotations"
    );
    tracing::info!(
        graph_dead_end_endpoints = fragmentation_report.summary.graph_dead_end_endpoints,
        branch_tangle_endpoints = fragmentation_report.summary.branch_tangle_endpoints,
        repeat_suspected_endpoints = fragmentation_report.summary.repeat_suspected_endpoints,
        empty_contig_paths = fragmentation_report.summary.empty_contig_paths,
        "contig fragmentation endpoint diagnosis"
    );
    tracing::info!(
        audit_findings = audit_report.findings.len(),
        low_support_kmers = audit_report.summary.low_support_kmers,
        low_support_regions = audit_report.summary.low_support_regions,
        abnormal_pair_hints = audit_report.summary.abnormal_pair_hints,
        collapsed_repeat_suspicions = audit_report.summary.collapsed_repeat_suspicions,
        "post-assembly audit report"
    );
    if params.diploid.enabled {
        tracing::info!(
            parent_refs = diploid_evidence.summary.parent_references_supplied,
            parent1_only_nodes = diploid_evidence.summary.parent1_only_nodes,
            parent2_only_nodes = diploid_evidence.summary.parent2_only_nodes,
            mixed_unitigs = diploid_evidence.summary.mixed_unitigs,
            full_haplotype_fasta_claimed = diploid_evidence.summary.full_haplotype_fasta_claimed,
            "diploid parent-specific k-mer evidence"
        );
    }

    if let Some(ref c) = selected_ck {
        stages.run(AssemblyStage::WriteExportCheckpoint, || {
            c.write_export(selected_k, &unitig_records, &contig_records)
        })?;
    }

    stages.run(AssemblyStage::WriteOutputs, || {
        let stage_reports = stages.reports();
        AssemblyOutputWriter::new(&params.outputs).write_all(AssemblyOutputBundle {
            graph: &graph,
            unitig_records: &unitig_records,
            contig_records: &contig_records,
            unitig_paths: &unitig_paths,
            primary_contig_paths_gfa: &primary_contig_paths_gfa,
            evidence: &evidence,
            trust_report: &trust_report,
            graph_annotations: &graph_annotations,
            simplify_decisions: &simplify_decisions,
            scaffold_artifact: &scaffold_artifact,
            multi_k_selection: &multi_k_selection,
            fragmentation_report: &fragmentation_report,
            audit_report: &audit_report,
            diploid_evidence: &diploid_evidence,
            diploid_enabled: params.diploid.enabled,
            stage_reports: &stage_reports,
        })
    })?;

    Ok(AssembleResult {
        reads,
        trusted_kmers: trusted,
        total_unique_kmers: total_unique,
        trust_report,
        multi_k_selection,
        simplify_stats,
        simplify_decisions,
        evidence,
        graph_annotations,
        scaffold_artifact,
        fragmentation_report,
        audit_report,
        diploid_evidence,
        mate_bridge_stats,
        stage_reports: stages.reports(),
        unitig_count: unitig_records.len(),
        contig_count: contig_records.len(),
        outputs: params.outputs.clone(),
    })
}

fn load_reads(params: &AssembleParams) -> Result<(Vec<Read>, Option<usize>), TrexError> {
    let r1_raw = if input_is_fasta(&params.r1_path) {
        read_fasta_records_maybe_gzip(&params.r1_path)?
    } else {
        read_fastq_records_maybe_gzip(&params.r1_path)?
    };
    let r1_reads = crate::illumina::preprocess::preprocess_records(r1_raw)?;

    match &params.r2_path {
        None => Ok((r1_reads, None)),
        Some(p) => {
            let r2_raw = if input_is_fasta(p) {
                read_fasta_records_maybe_gzip(p)?
            } else {
                read_fastq_records_maybe_gzip(p)?
            };
            let r2_reads = crate::illumina::preprocess::preprocess_records(r2_raw)?;
            validate_pair_parity(&r1_reads, &r2_reads)?;
            let r1_count = r1_reads.len();
            let mut combined = r1_reads;
            combined.extend(r2_reads);
            Ok((combined, Some(r1_count)))
        }
    }
}
