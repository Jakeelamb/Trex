//! `trex` CLI: **Tokio** async I/O boundary; CPU work runs in **`spawn_blocking`** on the sync **`trex`** core.

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use serde::Deserialize;
use tokio::task::spawn_blocking;
use tracing_subscriber::EnvFilter;

use trex::illumina::pipeline::{
    assemble_illumina, AssembleOutputs, AssembleParams, DiploidParams, SimplifyOverrides,
};
use trex::{IngestError, TrexError};

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct SimplifyFilePartial {
    max_tip_bases: Option<usize>,
    tip_max_multiplicity: Option<u64>,
    max_bubble_vertices: Option<usize>,
    max_bubble_internal_bases: Option<usize>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct DiploidFilePartial {
    enabled: Option<bool>,
    insert_mean_bp: Option<u64>,
    insert_stddev_bp: Option<u64>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct AssembleFileConfig {
    r2: Option<PathBuf>,
    k: Option<usize>,
    trusted_threshold: Option<u64>,
    checkpoint_root: Option<PathBuf>,
    resume: Option<bool>,
    strict_checkpoints: Option<bool>,
    out_dir: Option<PathBuf>,
    unitigs_fasta: Option<PathBuf>,
    contigs_fasta: Option<PathBuf>,
    gfa: Option<PathBuf>,
    simplify: Option<SimplifyFilePartial>,
    diploid: Option<DiploidFilePartial>,
}

fn parse_assemble_file_config(text: &str) -> Result<AssembleFileConfig, toml::de::Error> {
    let value: toml::Value = toml::from_str(text)?;
    match value.get("assemble") {
        Some(inner) => AssembleFileConfig::deserialize(inner.clone()),
        None => AssembleFileConfig::deserialize(value),
    }
}

#[derive(Parser)]
#[command(
    name = "trex",
    version,
    about = "Trex genome assembler (Phase-1 Illumina; optional Phase-2 diploid experimental)"
)]
struct Cli {
    #[command(subcommand)]
    command: TopCommand,
}

#[derive(Subcommand)]
enum TopCommand {
    /// Illumina short-read modes.
    Illumina {
        #[command(subcommand)]
        cmd: IlluminaCmd,
    },
}

#[derive(Subcommand)]
enum IlluminaCmd {
    /// Ingest → preprocess → *k*-mer counts → DBG → unitigs → contigs → FASTA/GFA.
    Assemble {
        #[arg(long)]
        r1: PathBuf,
        #[arg(long)]
        r2: Option<PathBuf>,
        /// *k*-mer size (required unless set in `--config`).
        #[arg(short = 'k', long = "kmer-size")]
        k: Option<usize>,
        #[arg(short = 'T', long = "trusted-threshold")]
        trusted_threshold: Option<u64>,
        #[arg(long)]
        checkpoint_root: Option<PathBuf>,
        #[arg(long, default_value_t = false)]
        resume: bool,
        /// Force **off** resume even if the config file sets `resume = true`.
        #[arg(long, default_value_t = false)]
        no_resume: bool,
        #[arg(long, default_value_t = false)]
        strict_checkpoints: bool,
        /// Force **off** strict checkpoint verification even if the config sets it.
        #[arg(long, default_value_t = false)]
        no_strict_checkpoints: bool,
        #[arg(long)]
        out_dir: Option<PathBuf>,
        #[arg(long)]
        unitigs_fasta: Option<PathBuf>,
        #[arg(long)]
        contigs_fasta: Option<PathBuf>,
        #[arg(long)]
        gfa: Option<PathBuf>,
        /// Maximum tip chain length (bases) for clipping (**Phase-1 tip clipping**).
        #[arg(long)]
        max_tip_bases: Option<usize>,
        /// Remove a tip leaf when its trusted multiplicity is **≤** this floor.
        #[arg(long)]
        tip_max_multiplicity: Option<u64>,
        /// Maximum vertices in an automatic diamond bubble motif.
        #[arg(long)]
        max_bubble_vertices: Option<usize>,
        /// Internal sequence-span budget (bases) for automatic bubble resolution.
        #[arg(long)]
        max_bubble_internal_bases: Option<usize>,
        /// **Phase-2 Illumina diploid** (experimental): retain near-balanced diamond bubbles; tag GFA header.
        #[arg(long, default_value_t = false)]
        diploid: bool,
        /// Optional insert-size prior mean (bp); stored for checkpoint identity and future mate-aware logic.
        #[arg(long)]
        insert_mean_bp: Option<u64>,
        #[arg(long)]
        insert_stddev_bp: Option<u64>,
        #[arg(long)]
        config: Option<PathBuf>,
    },
}

#[tokio::main]
async fn main() -> std::process::ExitCode {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .init();

    let cli = Cli::parse();
    let res = match cli.command {
        TopCommand::Illumina { cmd } => match cmd {
            IlluminaCmd::Assemble {
                r1,
                r2,
                k,
                trusted_threshold,
                checkpoint_root,
                resume,
                no_resume,
                strict_checkpoints,
                no_strict_checkpoints,
                out_dir,
                unitigs_fasta,
                contigs_fasta,
                gfa,
                max_tip_bases,
                tip_max_multiplicity,
                max_bubble_vertices,
                max_bubble_internal_bases,
                diploid,
                insert_mean_bp,
                insert_stddev_bp,
                config,
            } => {
                let file_cfg = if let Some(path) = &config {
                    let bytes = match tokio::fs::read(path).await {
                        Ok(b) => b,
                        Err(e) => {
                            tracing::error!(path = %path.display(), error = %e, "failed to read config");
                            return std::process::ExitCode::from(1);
                        }
                    };
                    let text = match std::str::from_utf8(&bytes) {
                        Ok(t) => t,
                        Err(e) => {
                            tracing::error!(path = %path.display(), error = %e, "config is not valid UTF-8");
                            return std::process::ExitCode::from(1);
                        }
                    };
                    match parse_assemble_file_config(text) {
                        Ok(c) => {
                            tracing::info!(path = %path.display(), "loaded assemble defaults from config");
                            c
                        }
                        Err(e) => {
                            tracing::error!(path = %path.display(), error = %e, "invalid TOML config");
                            return std::process::ExitCode::from(1);
                        }
                    }
                } else {
                    AssembleFileConfig::default()
                };

                let k = k.or(file_cfg.k);
                let Some(k) = k else {
                    tracing::error!(
                        "k-mer size missing: pass `--kmer-size` / `-k` or set `k` in config (optionally under `[assemble]`)"
                    );
                    return std::process::ExitCode::from(1);
                };

                let trusted_threshold = trusted_threshold
                    .or(file_cfg.trusted_threshold)
                    .unwrap_or(2);
                let r2 = r2.or(file_cfg.r2);
                let checkpoint_root = checkpoint_root.or(file_cfg.checkpoint_root);
                let resume = if no_resume {
                    false
                } else {
                    resume || file_cfg.resume.unwrap_or(false)
                };
                let strict_checkpoints = if no_strict_checkpoints {
                    false
                } else {
                    strict_checkpoints || file_cfg.strict_checkpoints.unwrap_or(false)
                };
                let out_dir = out_dir
                    .or(file_cfg.out_dir)
                    .unwrap_or_else(|| PathBuf::from("."));
                let unitigs_fasta = unitigs_fasta
                    .or(file_cfg.unitigs_fasta)
                    .unwrap_or_else(|| PathBuf::from("unitigs.fa"));
                let contigs_fasta = contigs_fasta
                    .or(file_cfg.contigs_fasta)
                    .unwrap_or_else(|| PathBuf::from("contigs.fa"));
                let gfa_path = gfa
                    .or(file_cfg.gfa)
                    .unwrap_or_else(|| PathBuf::from("graph.gfa"));

                let simplify = SimplifyOverrides {
                    max_tip_bases: max_tip_bases.or(file_cfg.simplify.as_ref().and_then(|s| s.max_tip_bases)),
                    tip_max_multiplicity: tip_max_multiplicity.or(
                        file_cfg
                            .simplify
                            .as_ref()
                            .and_then(|s| s.tip_max_multiplicity),
                    ),
                    max_bubble_vertices: max_bubble_vertices.or(
                        file_cfg
                            .simplify
                            .as_ref()
                            .and_then(|s| s.max_bubble_vertices),
                    ),
                    max_bubble_internal_bases: max_bubble_internal_bases.or(
                        file_cfg
                            .simplify
                            .as_ref()
                            .and_then(|s| s.max_bubble_internal_bases),
                    ),
                };

                let diploid_params = DiploidParams {
                    enabled: diploid
                        || file_cfg
                            .diploid
                            .as_ref()
                            .and_then(|d| d.enabled)
                            .unwrap_or(false),
                    insert_mean_bp: insert_mean_bp.or(
                        file_cfg
                            .diploid
                            .as_ref()
                            .and_then(|d| d.insert_mean_bp),
                    ),
                    insert_stddev_bp: insert_stddev_bp.or(
                        file_cfg
                            .diploid
                            .as_ref()
                            .and_then(|d| d.insert_stddev_bp),
                    ),
                };

                let params = AssembleParams {
                    r1_path: r1,
                    r2_path: r2,
                    k,
                    trusted_threshold,
                    checkpoint_root,
                    resume,
                    strict_checkpoints,
                    simplify,
                    diploid: diploid_params,
                    outputs: AssembleOutputs {
                        out_dir,
                        unitigs_fasta,
                        contigs_fasta,
                        gfa_path,
                    },
                };
                run_assemble(params).await
            }
        },
    };

    match res {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            tracing::error!(error = %e, "trex failed");
            std::process::ExitCode::from(1)
        }
    }
}

async fn run_assemble(params: AssembleParams) -> Result<(), TrexError> {
    let log_diploid = params.diploid.enabled;
    let out = spawn_blocking(move || assemble_illumina(&params))
        .await
        .map_err(|join_err| {
            TrexError::Ingest(IngestError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                join_err,
            )))
        })??;
    tracing::info!(
        reads = out.reads.len(),
        unique_kmers = out.total_unique_kmers,
        trusted_kmers = out.trusted_kmers.len(),
        unitigs = out.unitig_count,
        contigs = out.contig_count,
        diploid = log_diploid,
        unitigs_fasta = %out.outputs.unitigs_path().display(),
        contigs_fasta = %out.outputs.contigs_path().display(),
        gfa = %out.outputs.gfa_path_resolved().display(),
        "assemble complete"
    );
    Ok(())
}
