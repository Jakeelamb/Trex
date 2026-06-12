//! End-to-end **Illumina** preprocess → *k*-mer counts → trusted DBG → simplify → unitigs → contigs → exports.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use serde::Serialize;
use tracing::instrument;

use crate::dbg::{
    annotate_graph, assert_no_self_loops, build_dbg, extract_unitigs, forward_representatives,
    primary_contig_paths_for_gfa, reference_contig_paths, run_simplification_schedule,
    stitch_sequence, unitig_adjacency_links, write_contigs_fasta, write_gfa1, write_unitigs_fasta,
    ContigWalkTieBreak, DbgGraph, DiploidSimplifyMode, GfaWriteOptions, GraphAnnotations,
    SimplifyDecisionLog, SimplifyParams, SimplifyStats,
};
use crate::error::{GraphError, IngestError, TrexError};
use crate::evidence::EvidenceLedger;
use crate::illumina::audit::{audit_assembly, audit_tsv, AssemblyAuditReport};
use crate::illumina::checkpoint::{self, CheckpointRoot, GraphCheckpointIdentity};
use crate::illumina::counts::enumerate_sorted_counts;
use crate::illumina::diploid::{
    build_diploid_evidence, gfa_path_tags, gfa_unitig_tags, DiploidEvidenceReport,
    ParentReferenceParams,
};
use crate::illumina::fasta::parse_fasta;
use crate::illumina::fastq::parse_fastq;
use crate::illumina::fragmentation::{diagnose_fragmentation, FragmentationReport};
use crate::illumina::io::read_maybe_gzip;
use crate::illumina::mate;
use crate::illumina::multik::{select_k, MultiKParams, MultiKSelectionReport};
use crate::illumina::paired::validate_pair_parity;
use crate::illumina::phase2_primary;
use crate::illumina::read::Read;
use crate::illumina::scaffold::{
    build_scaffold_artifact, scaffold_fasta_records, scaffold_gfa_paths, ScaffoldArtifact,
};
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
}

/// Optional overrides for **Phase-1 graph simplification** (tips + bounded diamond bubbles).
#[derive(Debug, Clone, Default)]
pub struct SimplifyOverrides {
    pub max_tip_bases: Option<usize>,
    pub tip_max_multiplicity: Option<u64>,
    pub max_bubble_vertices: Option<usize>,
    pub max_bubble_internal_bases: Option<usize>,
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
    let ck = params
        .checkpoint_root
        .as_ref()
        .map(|p| CheckpointRoot::new(p.clone()));

    let (reads, paired_r1_len) = if params.resume {
        match &ck {
            Some(c) => {
                match checkpoint::load_preprocess_checkpoint(c, params.strict_checkpoints)? {
                    Some(r) => {
                        let pl = checkpoint::load_pair_layout_checkpoint(c)?;
                        (r, pl)
                    }
                    None => load_reads(params)?,
                }
            }
            None => return Err(IngestError::ResumeRequiresCheckpointRoot.into()),
        }
    } else {
        load_reads(params)?
    };

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

    let multi_k_selection = select_k(&reads, params.k, params.trusted_threshold, &params.multi_k)?;
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
            c.selected_k_root(selected_k)
        } else {
            c.clone()
        }
    });

    if let Some(ref c) = ck {
        checkpoint::write_preprocess_checkpoint(
            c,
            &reads,
            params.strict_checkpoints,
            paired_r1_len,
        )?;
    }

    let merged = if params.resume {
        match &selected_ck {
            Some(c) => {
                match checkpoint::load_counts_checkpoint(c, params.strict_checkpoints, selected_k)?
                {
                    Some(rows) => rows,
                    None => enumerate_sorted_counts(&reads, selected_k)?,
                }
            }
            None => enumerate_sorted_counts(&reads, selected_k)?,
        }
    } else {
        enumerate_sorted_counts(&reads, selected_k)?
    };

    let total_unique = merged.len();
    let trust_report = build_trust_diagnostics(selected_k, params.trusted_threshold, &merged);

    if let Some(ref c) = selected_ck {
        checkpoint::write_counts_checkpoint(c, selected_k, &merged, params.strict_checkpoints)?;
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

    let graph = if let Some(ref c) = selected_ck {
        if params.resume {
            match checkpoint::load_graph_checkpoint(
                c,
                params.strict_checkpoints,
                selected_k,
                &graph_ck_id,
            )? {
                Some(g) => g,
                None => {
                    let mut g = build_dbg(&reads, selected_k, &trusted)?;
                    mate_bridge_stats = apply_phase2_mate_bridge(
                        &mut g,
                        &reads,
                        paired_r1_len,
                        selected_k,
                        params.diploid.insert_mean_bp,
                        params.diploid.insert_stddev_bp,
                        graph_ck_id.phase2_mate_bridge_v1,
                    );
                    (simplify_stats, simplify_decisions) =
                        run_simplification_schedule(&mut g, &simplify_params, diploid_simplify);
                    checkpoint::write_graph_checkpoint(
                        c,
                        &g,
                        params.strict_checkpoints,
                        &graph_ck_id,
                    )?;
                    g
                }
            }
        } else {
            let mut g = build_dbg(&reads, selected_k, &trusted)?;
            mate_bridge_stats = apply_phase2_mate_bridge(
                &mut g,
                &reads,
                paired_r1_len,
                selected_k,
                params.diploid.insert_mean_bp,
                params.diploid.insert_stddev_bp,
                graph_ck_id.phase2_mate_bridge_v1,
            );
            (simplify_stats, simplify_decisions) =
                run_simplification_schedule(&mut g, &simplify_params, diploid_simplify);
            checkpoint::write_graph_checkpoint(c, &g, params.strict_checkpoints, &graph_ck_id)?;
            g
        }
    } else {
        let mut g = build_dbg(&reads, selected_k, &trusted)?;
        mate_bridge_stats = apply_phase2_mate_bridge(
            &mut g,
            &reads,
            paired_r1_len,
            selected_k,
            params.diploid.insert_mean_bp,
            params.diploid.insert_stddev_bp,
            graph_ck_id.phase2_mate_bridge_v1,
        );
        (simplify_stats, simplify_decisions) =
            run_simplification_schedule(&mut g, &simplify_params, diploid_simplify);
        g
    };

    if let Some(stats) = mate_bridge_stats.as_ref() {
        evidence.push(stats.evidence_record());
    }

    assert_no_self_loops(&graph)?;
    tracing::info!(
        tips_removed = simplify_stats.tips_removed,
        diamond_bubbles_resolved = simplify_stats.diamond_bubbles_resolved,
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

    let stage_start = Instant::now();
    let (unitig_records, mut contig_records, unitig_paths, contig_paths) =
        match (&selected_ck, params.resume) {
            (Some(c), true) => {
                match checkpoint::load_export_checkpoint(c, params.strict_checkpoints, selected_k)?
                {
                    Some((ur, cr)) => {
                        let up = extract_unitigs(&graph);
                        let cp =
                            reference_contig_paths(&graph, &forward, selected_k, contig_tie_break)?;
                        (ur, cr, up, cp)
                    }
                    None => {
                        stitch_unitigs_and_contigs(&graph, &forward, selected_k, contig_tie_break)?
                    }
                }
            }
            _ => stitch_unitigs_and_contigs(&graph, &forward, selected_k, contig_tie_break)?,
        };
    tracing::info!(
        stage = "stitch_unitigs_and_contigs",
        elapsed_ms = stage_start.elapsed().as_millis(),
        unitigs = unitig_records.len(),
        contigs = contig_records.len(),
        "illumina pipeline stage complete"
    );

    if params.diploid.enabled {
        let stage_start = Instant::now();
        let trusted_map: HashMap<Vec<u8>, u64> = trusted.iter().cloned().collect();
        for (_hdr, seq) in contig_records.iter_mut() {
            phase2_primary::collapse_primary_contig_by_trusted_kmers(seq, selected_k, &trusted_map);
        }
        tracing::info!(
            elapsed_ms = stage_start.elapsed().as_millis(),
            "Phase-2 Illumina: primary contig FASTA column collapse applied (trusted k-mer multiplicity)"
        );
    }

    let stage_start = Instant::now();
    let primary_contig_paths_gfa = primary_contig_paths_for_gfa(&contig_paths, &unitig_paths);
    tracing::info!(
        stage = "primary_contig_paths_for_gfa",
        elapsed_ms = stage_start.elapsed().as_millis(),
        gfa_paths = primary_contig_paths_gfa.len(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let graph_annotations = annotate_graph(&graph, &unitig_paths);
    tracing::info!(
        stage = "annotate_graph",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let fragmentation_report = diagnose_fragmentation(&graph, &contig_paths, &graph_annotations);
    tracing::info!(
        stage = "diagnose_fragmentation",
        elapsed_ms = stage_start.elapsed().as_millis(),
        graph_dead_end_endpoints = fragmentation_report.summary.graph_dead_end_endpoints,
        branch_tangle_endpoints = fragmentation_report.summary.branch_tangle_endpoints,
        repeat_suspected_endpoints = fragmentation_report.summary.repeat_suspected_endpoints,
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let scaffold_artifact = build_scaffold_artifact(
        mate_bridge_stats
            .as_ref()
            .map(|stats| stats.candidates.as_slice())
            .unwrap_or(&[]),
        &unitig_paths,
        selected_k,
        &fragmentation_report,
    );
    tracing::info!(
        stage = "build_scaffold_artifact",
        elapsed_ms = stage_start.elapsed().as_millis(),
        scaffold_paths = scaffold_artifact.paths.len(),
        bridge_candidates = scaffold_artifact.bridge_candidates.len(),
        endpoint_join_candidates = scaffold_artifact.endpoint_join_candidates.len(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let diploid_evidence = build_diploid_evidence(
        params.diploid.enabled,
        selected_k,
        &params.diploid.parent_references,
        &unitig_paths,
        &contig_records,
    )?;
    tracing::info!(
        stage = "build_diploid_evidence",
        elapsed_ms = stage_start.elapsed().as_millis(),
        parent_informative_unitigs = diploid_evidence.summary.parent_informative_unitigs,
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let audit_report = audit_assembly(
        selected_k,
        &contig_records,
        &trusted,
        &graph_annotations,
        mate_bridge_stats.as_ref(),
    );
    tracing::info!(
        stage = "audit_assembly",
        elapsed_ms = stage_start.elapsed().as_millis(),
        audit_findings = audit_report.findings.len(),
        "illumina pipeline stage complete"
    );
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
        let stage_start = Instant::now();
        checkpoint::write_export_checkpoint(
            c,
            selected_k,
            &unitig_records,
            &contig_records,
            params.strict_checkpoints,
        )?;
        tracing::info!(
            stage = "write_export_checkpoint",
            elapsed_ms = stage_start.elapsed().as_millis(),
            "illumina pipeline stage complete"
        );
    }

    let stage_start = Instant::now();
    if params.outputs.out_dir.as_os_str() != "-" {
        std::fs::create_dir_all(&params.outputs.out_dir)
            .map_err(|e| TrexError::from(GraphError::from(e)))?;
    }

    let uf = params.outputs.unitigs_path();
    let cf = params.outputs.contigs_path();
    let gf = params.outputs.gfa_path_resolved();

    for path in [&uf, &cf, &gf] {
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
    tracing::info!(
        stage = "prepare_output_directories",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );

    let stage_start = Instant::now();
    write_unitigs_fasta(&uf, &unitig_records)?;
    tracing::info!(
        stage = "write_unitigs_fasta",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    write_contigs_fasta(&cf, &contig_records)?;
    tracing::info!(
        stage = "write_contigs_fasta",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    let diploid_gfa_links = if params.diploid.enabled {
        Some(unitig_adjacency_links(&graph, &unitig_paths))
    } else {
        None
    };
    let parent_unitig_tags = gfa_unitig_tags(&diploid_evidence);
    let parent_path_tags = gfa_path_tags(&diploid_evidence);
    let scaffold_paths_gfa = scaffold_gfa_paths(&scaffold_artifact);
    tracing::info!(
        stage = "prepare_gfa_metadata",
        elapsed_ms = stage_start.elapsed().as_millis(),
        diploid_links = diploid_gfa_links
            .as_ref()
            .map(|links| links.len())
            .unwrap_or(0),
        parent_unitig_tags = parent_unitig_tags.len(),
        parent_path_tags = parent_path_tags.len(),
        scaffold_paths = scaffold_paths_gfa.len(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    write_gfa1(
        &gf,
        &unitig_records,
        GfaWriteOptions {
            phase2_illumina_diploid: params.diploid.enabled,
            diploid_unitig_links: diploid_gfa_links.as_deref(),
            primary_contig_paths: &primary_contig_paths_gfa,
            scaffold_paths: &scaffold_paths_gfa,
            phase2_unphased_hap_paths: params.diploid.enabled,
            parent_unitig_tags: Some(parent_unitig_tags.as_slice()),
            parent_path_tags: Some(parent_path_tags.as_slice()),
        },
    )?;
    tracing::info!(
        stage = "write_gfa1",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );
    let stage_start = Instant::now();
    write_evidence_json(&params.outputs.evidence_path(), &evidence)?;
    write_trust_json(&params.outputs.trust_path(), &trust_report)?;
    write_annotations_json(&params.outputs.annotations_path(), &graph_annotations)?;
    write_simplification_json(&params.outputs.simplification_path(), &simplify_decisions)?;
    write_scaffolds_json(&params.outputs.scaffolds_path(), &scaffold_artifact)?;
    write_scaffolds_fasta(
        &params.outputs.scaffolds_fasta_path(),
        &scaffold_fasta_records(&scaffold_artifact, &unitig_records),
    )?;
    write_multi_k_json(&params.outputs.multi_k_path(), &multi_k_selection)?;
    write_fragmentation_json(&params.outputs.fragmentation_path(), &fragmentation_report)?;
    write_audit_json(&params.outputs.audit_json_path(), &audit_report)?;
    write_audit_tsv(&params.outputs.audit_tsv_path(), &audit_report)?;
    write_diploid_json(&params.outputs.diploid_path(), &diploid_evidence)?;
    tracing::info!(
        stage = "write_json_sidecars",
        elapsed_ms = stage_start.elapsed().as_millis(),
        "illumina pipeline stage complete"
    );

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
        unitig_count: unitig_records.len(),
        contig_count: contig_records.len(),
        outputs: params.outputs.clone(),
    })
}

fn write_evidence_json(path: &Path, evidence: &EvidenceLedger) -> Result<(), TrexError> {
    write_json_pretty(path, evidence)
}

fn write_trust_json(path: &Path, report: &TrustDiagnosticsReport) -> Result<(), TrexError> {
    write_json_pretty(path, report)
}

fn write_annotations_json(path: &Path, annotations: &GraphAnnotations) -> Result<(), TrexError> {
    write_json_pretty(path, annotations)
}

fn write_simplification_json(
    path: &Path,
    decisions: &SimplifyDecisionLog,
) -> Result<(), TrexError> {
    write_json_pretty(path, decisions)
}

fn write_scaffolds_json(path: &Path, artifact: &ScaffoldArtifact) -> Result<(), TrexError> {
    write_json_pretty(path, artifact)
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

fn write_fragmentation_json(path: &Path, report: &FragmentationReport) -> Result<(), TrexError> {
    write_json_pretty(path, report)
}

fn write_audit_json(path: &Path, report: &AssemblyAuditReport) -> Result<(), TrexError> {
    write_json_pretty(path, report)
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

fn write_diploid_json(path: &Path, report: &DiploidEvidenceReport) -> Result<(), TrexError> {
    write_json_pretty(path, report)
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

fn load_reads(params: &AssembleParams) -> Result<(Vec<Read>, Option<usize>), TrexError> {
    let r1_bytes =
        read_maybe_gzip(&params.r1_path).map_err(|e| TrexError::Ingest(IngestError::Io(e)))?;
    let r1_raw = if input_is_fasta(&params.r1_path) {
        parse_fasta(&r1_bytes)?
    } else {
        parse_fastq(&r1_bytes)?
    };
    let r1_reads = crate::illumina::preprocess::preprocess_records(r1_raw)?;

    match &params.r2_path {
        None => Ok((r1_reads, None)),
        Some(p) => {
            let r2_bytes = read_maybe_gzip(p).map_err(|e| TrexError::Ingest(IngestError::Io(e)))?;
            let r2_raw = if input_is_fasta(p) {
                parse_fasta(&r2_bytes)?
            } else {
                parse_fastq(&r2_bytes)?
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
