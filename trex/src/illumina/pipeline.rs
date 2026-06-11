//! End-to-end **Illumina** preprocess → *k*-mer counts → trusted DBG → simplify → unitigs → contigs → exports.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use tracing::instrument;

use crate::dbg::{
    assert_no_self_loops, build_dbg, extract_unitigs, forward_representatives,
    primary_contig_paths_for_gfa, reference_contig_paths, remove_diamond_bubbles_ext, remove_tips,
    stitch_sequence, unitig_adjacency_links, write_contigs_fasta, write_gfa1, write_unitigs_fasta,
    ContigWalkTieBreak, DiploidSimplifyMode, SimplifyParams,
};
use crate::error::{GraphError, IngestError, TrexError};
use crate::illumina::checkpoint::{self, CheckpointRoot, GraphCheckpointIdentity};
use crate::illumina::counts::enumerate_sorted_counts;
use crate::illumina::fasta::parse_fasta;
use crate::illumina::fastq::parse_fastq;
use crate::illumina::io::read_maybe_gzip;
use crate::illumina::mate;
use crate::illumina::paired::validate_pair_parity;
use crate::illumina::phase2_primary;
use crate::illumina::read::Read;
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
    pub outputs: AssembleOutputs,
}

/// Outputs after a full **assemble** run.
#[derive(Debug, Clone)]
pub struct AssembleResult {
    pub reads: Vec<Read>,
    pub trusted_kmers: Vec<(Vec<u8>, u64)>,
    pub total_unique_kmers: usize,
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
        if s < params.k {
            return Err(IngestError::KTooLarge {
                k: params.k,
                shortest: s,
            }
            .into());
        }
    } else {
        return Err(IngestError::FastqFormat("no reads ingested".into()).into());
    }

    if let Some(ref c) = ck {
        checkpoint::write_preprocess_checkpoint(
            c,
            &reads,
            params.strict_checkpoints,
            paired_r1_len,
        )?;
    }

    let merged = if params.resume {
        match &ck {
            Some(c) => {
                match checkpoint::load_counts_checkpoint(c, params.strict_checkpoints, params.k)? {
                    Some(rows) => rows,
                    None => enumerate_sorted_counts(&reads, params.k)?,
                }
            }
            None => enumerate_sorted_counts(&reads, params.k)?,
        }
    } else {
        enumerate_sorted_counts(&reads, params.k)?
    };

    let total_unique = merged.len();

    if let Some(ref c) = ck {
        checkpoint::write_counts_checkpoint(c, params.k, &merged, params.strict_checkpoints)?;
    }

    let trusted = apply_trusted_threshold(merged, params.trusted_threshold);

    let forward = forward_representatives(&reads, params.k)?;
    let simplify_params = merge_simplify_params(params.k, &params.simplify);
    let diploid_simplify = params.diploid.enabled.then_some(DiploidSimplifyMode);
    let graph_ck_id = graph_checkpoint_identity(params);

    let graph = if let Some(ref c) = ck {
        if params.resume {
            match checkpoint::load_graph_checkpoint(
                c,
                params.strict_checkpoints,
                params.k,
                &graph_ck_id,
            )? {
                Some(g) => g,
                None => {
                    let mut g = build_dbg(&reads, params.k, &trusted)?;
                    if graph_ck_id.phase2_mate_bridge_v1 {
                        if let Some(n) = paired_r1_len {
                            let nb = mate::boost_mate_pairs_on_existing_dbg_edges(
                                &mut g, &reads, n, params.k,
                            );
                            if nb > 0 {
                                tracing::info!(
                                    mate_bridge_boosts = nb,
                                    "Phase-2 Illumina mate-pair edge boosts (existing DBG edges only)"
                                );
                            }
                        } else {
                            tracing::warn!(
                                "Phase-2 mate bridge skipped: missing preprocess/pair_layout.json; re-run without --resume to record paired layout"
                            );
                        }
                    }
                    remove_tips(&mut g, &simplify_params);
                    remove_diamond_bubbles_ext(&mut g, &simplify_params, diploid_simplify);
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
            let mut g = build_dbg(&reads, params.k, &trusted)?;
            if graph_ck_id.phase2_mate_bridge_v1 {
                if let Some(n) = paired_r1_len {
                    let nb =
                        mate::boost_mate_pairs_on_existing_dbg_edges(&mut g, &reads, n, params.k);
                    if nb > 0 {
                        tracing::info!(
                            mate_bridge_boosts = nb,
                            "Phase-2 Illumina mate-pair edge boosts (existing DBG edges only)"
                        );
                    }
                } else {
                    tracing::warn!(
                        "Phase-2 mate bridge skipped: missing preprocess/pair_layout.json; re-run without --resume to record paired layout"
                    );
                }
            }
            remove_tips(&mut g, &simplify_params);
            remove_diamond_bubbles_ext(&mut g, &simplify_params, diploid_simplify);
            checkpoint::write_graph_checkpoint(c, &g, params.strict_checkpoints, &graph_ck_id)?;
            g
        }
    } else {
        let mut g = build_dbg(&reads, params.k, &trusted)?;
        if graph_ck_id.phase2_mate_bridge_v1 {
            if let Some(n) = paired_r1_len {
                let nb = mate::boost_mate_pairs_on_existing_dbg_edges(&mut g, &reads, n, params.k);
                if nb > 0 {
                    tracing::info!(
                        mate_bridge_boosts = nb,
                        "Phase-2 Illumina mate-pair edge boosts (existing DBG edges only)"
                    );
                }
            } else {
                tracing::warn!(
                    "Phase-2 mate bridge skipped: missing preprocess/pair_layout.json; re-run without --resume to record paired layout"
                );
            }
        }
        remove_tips(&mut g, &simplify_params);
        remove_diamond_bubbles_ext(&mut g, &simplify_params, diploid_simplify);
        g
    };

    assert_no_self_loops(&graph)?;

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
        match (&ck, params.resume) {
            (Some(c), true) => {
                match checkpoint::load_export_checkpoint(c, params.strict_checkpoints, params.k)? {
                    Some((ur, cr)) => {
                        let up = extract_unitigs(&graph);
                        let cp =
                            reference_contig_paths(&graph, &forward, params.k, contig_tie_break)?;
                        (ur, cr, up, cp)
                    }
                    None => {
                        stitch_unitigs_and_contigs(&graph, &forward, params.k, contig_tie_break)?
                    }
                }
            }
            _ => stitch_unitigs_and_contigs(&graph, &forward, params.k, contig_tie_break)?,
        };

    if params.diploid.enabled {
        let trusted_map: HashMap<Vec<u8>, u64> = trusted.iter().cloned().collect();
        for (_hdr, seq) in contig_records.iter_mut() {
            phase2_primary::collapse_primary_contig_by_trusted_kmers(seq, params.k, &trusted_map);
        }
        tracing::info!(
            "Phase-2 Illumina: primary contig FASTA column collapse applied (trusted k-mer multiplicity)"
        );
    }

    let primary_contig_paths_gfa = primary_contig_paths_for_gfa(&contig_paths, &unitig_paths);

    if let Some(ref c) = ck {
        checkpoint::write_export_checkpoint(
            c,
            params.k,
            &unitig_records,
            &contig_records,
            params.strict_checkpoints,
        )?;
    }

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

    write_unitigs_fasta(&uf, &unitig_records)?;
    write_contigs_fasta(&cf, &contig_records)?;
    let diploid_gfa_links = if params.diploid.enabled {
        Some(unitig_adjacency_links(&graph, &unitig_paths))
    } else {
        None
    };
    write_gfa1(
        &gf,
        &unitig_records,
        params.diploid.enabled,
        diploid_gfa_links.as_deref(),
        &primary_contig_paths_gfa,
        params.diploid.enabled,
    )?;

    Ok(AssembleResult {
        reads,
        trusted_kmers: trusted,
        total_unique_kmers: total_unique,
        unitig_count: unitig_records.len(),
        contig_count: contig_records.len(),
        outputs: params.outputs.clone(),
    })
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
