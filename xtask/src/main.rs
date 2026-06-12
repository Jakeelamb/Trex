use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use sha2::{Digest, Sha256};

type DynResult<T> = Result<T, Box<dyn Error>>;

#[derive(Parser)]
#[command(name = "xtask", version, about = "Trex repository automation")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run all repository contract validators.
    Validate,
    /// Validate tools/benchmark_matrix.toml.
    ValidateMatrix,
    /// Validate docs/CAPABILITIES.md against CLI flags and scripts.
    ValidateCapabilities,
    /// Validate tools/benchmark_data.toml.
    ValidateData,
    /// Validate tools/assembler_framework.toml.
    ValidateFramework,
    /// Validate docs/DEVELOPMENT.md and docs/ROADMAP.md.
    ValidateDevelopment,
    /// Run benchmark matrix scripts for one tier and write a JSON artifact.
    Bench {
        #[arg(long, value_enum, default_value_t = Tier::Pr)]
        tier: Tier,
        #[arg(long)]
        row: Option<String>,
        #[arg(long, default_value = "target/benchmarks/matrix.json")]
        out: PathBuf,
    },
    /// Generate deterministic FASTQ reads from a FASTA reference.
    GenerateReads {
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value_t = 100)]
        read_len: usize,
        #[arg(long, default_value_t = 25)]
        step: usize,
        #[arg(long, default_value_t = false)]
        circular: bool,
    },
    /// Run a Rust-owned gate for a development tier.
    Gate {
        #[arg(long, value_enum, default_value_t = Tier::Pr)]
        tier: Tier,
    },
    /// Fetch and prepare ignored external benchmark data subsets.
    FetchData {
        #[arg(long)]
        dataset: Option<String>,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
enum Tier {
    Pr,
    Main,
    Nightly,
    Manual,
}

#[derive(Debug, Deserialize)]
struct Matrix {
    schema_version: u64,
    #[serde(default)]
    rows: Vec<Row>,
}

#[derive(Debug, Deserialize)]
struct Row {
    id: Option<String>,
    technology: Option<String>,
    organism: Option<String>,
    license: Option<String>,
    provenance: Option<String>,
    depth_class: Option<String>,
    ci_tier: Option<Tier>,
    product_claim: Option<String>,
    claim_level: Option<String>,
    reference_availability: Option<String>,
    artifact_policy: Option<String>,
    fixtures: Option<Vec<String>>,
    external_data: Option<String>,
    digest_manifest: Option<String>,
    manifest_table: Option<String>,
    pr_scripts: Option<Vec<String>>,
    main_scripts: Option<Vec<String>>,
    nightly_scripts: Option<Vec<String>>,
    manual_scripts: Option<Vec<String>>,
    required_tools: Option<Vec<String>>,
    optional_tools: Option<Vec<String>>,
    artifacts: Option<Vec<String>>,
    pr_artifacts: Option<Vec<String>>,
    main_artifacts: Option<Vec<String>>,
    nightly_artifacts: Option<Vec<String>>,
    manual_artifacts: Option<Vec<String>>,
    trex: Option<TrexBench>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct TrexBench {
    tiers: Vec<Tier>,
    args: Vec<String>,
    out_dir: String,
    reference: Option<String>,
    #[serde(default)]
    parental_references: Vec<ParentReference>,
}

#[derive(Debug, Deserialize)]
struct ParentReference {
    label: String,
    path: String,
}

#[derive(Debug, Deserialize)]
struct DataCatalog {
    schema_version: u64,
    #[serde(default)]
    datasets: Vec<DataSet>,
}

#[derive(Debug, Deserialize)]
struct DataSet {
    id: String,
    scientific_name: String,
    source: String,
    study_accession: String,
    run_accession: String,
    sample_accession: String,
    experiment_accession: String,
    library_strategy: String,
    library_layout: String,
    instrument_model: String,
    read_count: u64,
    base_count: u64,
    ploidy: String,
    license: String,
    provenance_url: String,
    ploidy_provenance_url: Option<String>,
    notes: String,
    #[serde(default)]
    files: Vec<DataFile>,
    #[serde(default)]
    alignments: Vec<DataAlignmentFile>,
    #[serde(default)]
    references: Vec<DataReference>,
    #[serde(default)]
    reference_slices: Vec<DataReferenceSlice>,
    #[serde(default)]
    truth: Vec<DataTruthFile>,
    #[serde(default)]
    interval_pairs: Vec<IntervalPairDataFile>,
    #[serde(default)]
    local_files: Vec<LocalDataFile>,
    #[serde(default)]
    prepared: Vec<PreparedDataFile>,
}

#[derive(Debug, Deserialize)]
struct DataFile {
    role: String,
    url: String,
    md5: Option<String>,
    bytes: u64,
}

#[derive(Debug, Deserialize)]
struct LocalDataFile {
    role: String,
    source_path: String,
    path: String,
    records: Option<usize>,
    bases: Option<u64>,
    sha256: Option<String>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct PreparedDataFile {
    role: String,
    source_role: String,
    path: String,
    records: usize,
    sha256: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DataAlignmentFile {
    role: String,
    url: String,
    index_url: Option<String>,
    bytes: Option<u64>,
    index_bytes: Option<u64>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct IntervalPairDataFile {
    role: String,
    source_role: String,
    region: String,
    records: usize,
    r1_path: String,
    r2_path: String,
    r1_sha256: Option<String>,
    r2_sha256: Option<String>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DataReference {
    role: String,
    url: String,
    path: String,
    sha256: Option<String>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DataReferenceSlice {
    role: String,
    source_role: String,
    region: String,
    path: String,
    sha256: Option<String>,
    notes: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DataTruthFile {
    role: String,
    url: String,
    path: String,
    md5: Option<String>,
    bytes: Option<u64>,
    notes: Option<String>,
}

#[derive(Debug, Eq, PartialEq)]
struct GenomicRegion {
    contig: String,
    start: usize,
    end: usize,
}

#[derive(Debug, Deserialize)]
struct AssemblerFramework {
    schema_version: u64,
    #[serde(default)]
    modules: Vec<FrameworkModule>,
}

#[derive(Debug, Deserialize)]
struct FrameworkModule {
    id: String,
    phase: String,
    status: String,
    summary: String,
    papers: Vec<String>,
    current_files: Vec<String>,
    promotion_stages: Vec<String>,
    gates: Vec<String>,
    next_work: Vec<String>,
}

const PROMOTION_STAGES: &[&str] = &[
    "report_only_candidate",
    "scaffold_artifact",
    "gfa_path",
    "fasta_scaffold_with_gaps",
    "graph_edit",
    "polishing_edit",
];

#[derive(Serialize)]
struct BenchReport {
    schema_version: u64,
    tier: Tier,
    generated_unix_ms: u128,
    rows: Vec<BenchRowReport>,
}

#[derive(Serialize)]
struct BenchRowReport {
    id: String,
    ci_tier: Tier,
    product_claim: String,
    claim_level: String,
    reference_availability: String,
    artifact_policy: String,
    tools: ToolAvailability,
    layer_status: RowLayerStatus,
    scripts: Vec<ScriptReport>,
    trex_runs: Vec<TrexRunReport>,
    artifacts: Vec<ArtifactReport>,
}

#[derive(Serialize)]
struct ToolAvailability {
    required_available: Vec<String>,
    required_unavailable: Vec<String>,
    optional_available: Vec<String>,
    optional_unavailable: Vec<String>,
}

#[derive(Serialize)]
struct RowLayerStatus {
    status: String,
    failed_layer: Option<String>,
}

#[derive(Serialize)]
struct ScriptReport {
    path: String,
    exit_code: Option<i32>,
    success: bool,
    wall_ms: u128,
    time_seconds: Option<f64>,
    max_rss_kib: Option<u64>,
}

#[derive(Serialize)]
struct TrexRunReport {
    command: Vec<String>,
    exit_code: Option<i32>,
    success: bool,
    wall_ms: u128,
    time_seconds: Option<f64>,
    max_rss_kib: Option<u64>,
    observed: Option<TrexObserved>,
    metrics: Option<AssemblyMetrics>,
    quast: Option<ScriptReport>,
}

#[derive(Serialize)]
struct TrexObserved {
    reads: Option<usize>,
    unique_kmers: Option<usize>,
    trusted_kmers: Option<usize>,
    unitigs: Option<usize>,
    contigs: Option<usize>,
}

#[derive(Serialize)]
struct AssemblyMetrics {
    kmer_size: Option<usize>,
    candidate_kmers: Option<usize>,
    reads: Option<FastqStats>,
    r2_reads: Option<FastqStats>,
    reference: Option<FastaStats>,
    contigs: Option<FastaStats>,
    unitigs: Option<FastaStats>,
    gfa: Option<GfaStats>,
    evidence: Option<Value>,
    annotations: Option<Value>,
    simplification: Option<Value>,
    scaffolds: Option<Value>,
    multi_k: Option<Value>,
    fragmentation: Option<Value>,
    audit: Option<Value>,
    diploid: Option<Value>,
    reference_quality: Option<ReferenceQuality>,
    parental_reference_quality: Vec<NamedReferenceQuality>,
    read_assembly_quality: Option<ReadAssemblyKmerQuality>,
}

#[derive(Serialize)]
struct FastqStats {
    records: usize,
    total_bases: usize,
    min_len: usize,
    max_len: usize,
    candidate_kmers: Option<usize>,
}

#[derive(Clone, Serialize)]
struct FastaStats {
    records: usize,
    total_bases: usize,
    min_len: usize,
    max_len: usize,
    n50: usize,
}

#[derive(Serialize)]
struct GfaStats {
    s_lines: usize,
    l_lines: usize,
    p_lines: usize,
    bytes: u64,
}

#[derive(Serialize)]
struct ReferenceQuality {
    kmer_size: usize,
    contig_kmers: usize,
    contig_kmers_in_reference: usize,
    contig_kmer_reference_fraction: f64,
    contigs_with_reference_hit: usize,
    contig_records: usize,
    contig_bases: usize,
    contig_bases_with_reference_hit: usize,
}

#[derive(Serialize)]
struct NamedReferenceQuality {
    label: String,
    path: String,
    reference: FastaStats,
    quality: ReferenceQuality,
}

#[derive(Serialize)]
struct ReadAssemblyKmerQuality {
    kmer_size: usize,
    reliable_threshold: usize,
    read_distinct_kmers: usize,
    reliable_read_kmers: usize,
    assembly_distinct_kmers: usize,
    assembly_total_kmers: usize,
    reliable_read_kmers_in_assembly: usize,
    reliable_read_containment_fraction: f64,
    assembly_only_kmers: usize,
    assembly_only_total_kmers: usize,
    assembly_only_fraction: f64,
    assembly_only_total_fraction: f64,
    qv_error_rate: f64,
    approximate_qv: Option<f64>,
}

#[derive(Serialize)]
struct ArtifactReport {
    path: String,
    exists: bool,
    bytes: Option<u64>,
}

fn main() -> DynResult<()> {
    let cli = Cli::parse();
    let root = repo_root()?;
    match cli.cmd {
        Cmd::Validate => {
            validate_matrix(&root)?;
            validate_capabilities(&root)?;
            validate_data_catalog(&root)?;
            validate_assembler_framework(&root)?;
            validate_development_docs(&root)?;
        }
        Cmd::ValidateMatrix => validate_matrix(&root)?,
        Cmd::ValidateCapabilities => validate_capabilities(&root)?,
        Cmd::ValidateData => validate_data_catalog(&root)?,
        Cmd::ValidateFramework => validate_assembler_framework(&root)?,
        Cmd::ValidateDevelopment => validate_development_docs(&root)?,
        Cmd::Bench { tier, row, out } => run_bench(&root, tier, row.as_deref(), &out)?,
        Cmd::GenerateReads {
            reference,
            out,
            read_len,
            step,
            circular,
        } => generate_reads(&root, &reference, &out, read_len, step, circular)?,
        Cmd::Gate { tier } => run_gate(&root, tier)?,
        Cmd::FetchData { dataset } => fetch_data(&root, dataset.as_deref())?,
    }
    Ok(())
}

fn repo_root() -> DynResult<PathBuf> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let root = manifest_dir
        .parent()
        .ok_or("xtask manifest directory has no parent")?;
    Ok(root.to_path_buf())
}

fn load_matrix(root: &Path) -> DynResult<Matrix> {
    let path = root.join("tools/benchmark_matrix.toml");
    let text = fs::read_to_string(path)?;
    Ok(toml::from_str(&text)?)
}

fn load_data_catalog(root: &Path) -> DynResult<DataCatalog> {
    let path = root.join("tools/benchmark_data.toml");
    let text = fs::read_to_string(path)?;
    Ok(toml::from_str(&text)?)
}

fn load_assembler_framework(root: &Path) -> DynResult<AssemblerFramework> {
    let path = root.join("tools/assembler_framework.toml");
    let text = fs::read_to_string(path)?;
    Ok(toml::from_str(&text)?)
}

fn validate_matrix(root: &Path) -> DynResult<()> {
    let matrix = load_matrix(root)?;
    let data_catalog = load_data_catalog(root).ok();
    let data_ids: BTreeSet<String> = data_catalog
        .as_ref()
        .map(|catalog| {
            catalog
                .datasets
                .iter()
                .map(|dataset| dataset.id.clone())
                .collect()
        })
        .unwrap_or_default();
    if matrix.schema_version != 1 {
        return Err("benchmark_matrix.toml schema_version must be 1".into());
    }
    if matrix.rows.is_empty() {
        return Err("benchmark_matrix.toml must contain at least one [[rows]] entry".into());
    }

    let mut seen = BTreeSet::new();
    for (idx, row) in matrix.rows.iter().enumerate() {
        let row_id = required_str(row.id.as_deref(), idx, "id")?;
        if !seen.insert(row_id.to_string()) {
            return Err(format!("duplicate benchmark row id {row_id:?}").into());
        }

        required_str(row.technology.as_deref(), idx, "technology")?;
        required_str(row.organism.as_deref(), idx, "organism")?;
        required_str(row.license.as_deref(), idx, "license")?;
        required_str(row.provenance.as_deref(), idx, "provenance")?;
        required_str(row.depth_class.as_deref(), idx, "depth_class")?;
        required_str(row.product_claim.as_deref(), idx, "product_claim")?;
        let claim_level = required_str(row.claim_level.as_deref(), idx, "claim_level")?;
        validate_claim_level(row_id, claim_level)?;
        let reference_availability = required_str(
            row.reference_availability.as_deref(),
            idx,
            "reference_availability",
        )?;
        validate_reference_availability(row_id, reference_availability)?;
        let artifact_policy = required_str(row.artifact_policy.as_deref(), idx, "artifact_policy")?;
        validate_artifact_policy(row_id, artifact_policy)?;
        let ci_tier = row
            .ci_tier
            .ok_or_else(|| format!("{row_id}: missing ci_tier"))?;

        let external_data = row.external_data.as_deref();
        if let Some(data_id) = external_data {
            if !data_ids.contains(data_id) {
                return Err(format!("{row_id}: unknown external_data id {data_id}").into());
            }
        }
        let is_external = external_data.is_some();

        let fixtures = required_list(row.fixtures.as_ref(), row_id, "fixtures")?;
        for fixture in fixtures {
            require_rel_path(root, row_id, "fixture", fixture, !is_external)?;
            if is_external && !fixture.starts_with("data/") {
                return Err(format!(
                    "{row_id}: external fixture path must live under data/: {fixture}"
                )
                .into());
            }
        }

        let script_paths = scripts_for_any_tier(row);
        if script_paths.is_empty() && row.trex.is_none() {
            return Err(format!(
                "{row_id}: at least one *_scripts list or [rows.trex] is required"
            )
            .into());
        }
        if ci_tier == Tier::Pr && row.pr_scripts.as_deref().unwrap_or(&[]).is_empty() {
            let has_pr_trex = row
                .trex
                .as_ref()
                .map(|trex| trex.tiers.contains(&Tier::Pr))
                .unwrap_or(false);
            if !has_pr_trex {
                return Err(
                    format!("{row_id}: ci_tier=pr requires pr_scripts or pr trex tier").into(),
                );
            }
        }
        for script in script_paths {
            require_rel_path(root, row_id, "script", script, true)?;
            if !script.starts_with("scripts/") {
                return Err(
                    format!("{row_id}: script path must live under scripts/: {script}").into(),
                );
            }
        }

        if let Some(manifest) = row.digest_manifest.as_deref() {
            require_rel_path(root, row_id, "digest_manifest", manifest, true)?;
            required_str(row.manifest_table.as_deref(), idx, "manifest_table")?;
            verify_manifest_digests(
                root,
                row_id,
                manifest,
                row.manifest_table.as_deref().unwrap(),
                fixtures,
                is_external,
            )?;
        }

        let artifacts = required_list(row.artifacts.as_ref(), row_id, "artifacts")?;
        for artifact in artifacts {
            if artifact.trim().is_empty() {
                return Err(format!("{row_id}: artifacts contains an empty entry").into());
            }
            require_rel_path(root, row_id, "artifact", artifact, false)?;
        }
        for artifact in artifacts_for_any_tier(row) {
            if artifact.trim().is_empty() {
                return Err(format!("{row_id}: tier artifacts contain an empty entry").into());
            }
            require_rel_path(root, row_id, "tier_artifact", artifact, false)?;
        }

        let required_tools = required_list(row.required_tools.as_ref(), row_id, "required_tools")?;
        validate_tool_list(row_id, "required_tools", Some(required_tools))?;
        validate_tool_list(row_id, "optional_tools", row.optional_tools.as_deref())?;
        if let Some(notes) = row.notes.as_deref() {
            if notes.trim().is_empty() {
                return Err(format!("{row_id}: notes must not be empty when present").into());
            }
        }
        if let Some(trex) = row.trex.as_ref() {
            validate_trex_bench(root, row_id, trex, is_external)?;
        }
    }

    println!("xtask validate-matrix: OK ({} rows)", matrix.rows.len());
    Ok(())
}

fn verify_manifest_digests(
    root: &Path,
    row_id: &str,
    manifest: &str,
    manifest_table: &str,
    fixtures: &[String],
    allow_missing_external: bool,
) -> DynResult<()> {
    let text = fs::read_to_string(root.join(manifest))?;
    let value: toml::Value = toml::from_str(&text)?;
    let mut table = &value;
    for part in manifest_table.split('.') {
        table = table
            .get(part)
            .ok_or_else(|| format!("{row_id}: missing manifest table {manifest_table}"))?;
    }
    for fixture in fixtures {
        let key = digest_key_for_fixture(fixture)?;
        let expected = table
            .get(&key)
            .and_then(toml::Value::as_str)
            .ok_or_else(|| format!("{row_id}: missing digest key {manifest_table}.{key}"))?;
        let path = root.join(fixture);
        if allow_missing_external && !path.exists() {
            continue;
        }
        let got = sha256_file(&path)?;
        if got != expected {
            return Err(format!(
                "{row_id}: digest mismatch for {fixture}: got {got}, expected {expected}"
            )
            .into());
        }
    }
    Ok(())
}

fn digest_key_for_fixture(path: &str) -> DynResult<String> {
    let path = Path::new(path);
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or("fixture path has no file name")?;
    Ok(format!("{}_sha256", name.replace('.', "_")))
}

fn sha256_file(path: &Path) -> DynResult<String> {
    let bytes = fs::read(path)?;
    let digest = Sha256::digest(bytes);
    Ok(format!("{digest:x}"))
}

fn validate_capabilities(root: &Path) -> DynResult<()> {
    let doc_path = root.join("docs/CAPABILITIES.md");
    let text = fs::read_to_string(doc_path)?;
    let flags = extract_assemble_flags(&root.join("trex-cli/src/main.rs"))?;
    let missing_flags: Vec<_> = flags
        .iter()
        .filter(|flag| !text.contains(flag.as_str()))
        .collect();
    if !missing_flags.is_empty() {
        return Err(format!("docs/CAPABILITIES.md missing CLI flags: {missing_flags:?}").into());
    }

    let mut scripts = Vec::new();
    for entry in fs::read_dir(root.join("scripts"))? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("sh") {
            if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                scripts.push(name.to_string());
            }
        }
    }
    scripts.sort();
    let missing_scripts: Vec<_> = scripts
        .iter()
        .filter(|name| !text.contains(name.as_str()))
        .collect();
    if !missing_scripts.is_empty() {
        return Err(format!("docs/CAPABILITIES.md missing scripts: {missing_scripts:?}").into());
    }

    for phrase in [
        "Phase-1 default",
        "Phase-2 Illumina --diploid",
        "Future / deferred",
        "tools/benchmark_matrix.toml",
        "cargo run -p xtask -- bench",
    ] {
        if !text.contains(phrase) {
            return Err(format!("docs/CAPABILITIES.md missing required phrase: {phrase}").into());
        }
    }

    println!(
        "xtask validate-capabilities: OK ({} scripts, {} flags)",
        scripts.len(),
        flags.len()
    );
    Ok(())
}

fn validate_data_catalog(root: &Path) -> DynResult<()> {
    let catalog = load_data_catalog(root)?;
    if catalog.schema_version != 1 {
        return Err("benchmark_data.toml schema_version must be 1".into());
    }
    if catalog.datasets.is_empty() {
        return Err("benchmark_data.toml must contain at least one [[datasets]] entry".into());
    }

    let mut ids = BTreeSet::new();
    for dataset in &catalog.datasets {
        if dataset.id.trim().is_empty() {
            return Err("benchmark_data.toml dataset id must not be empty".into());
        }
        if !ids.insert(dataset.id.clone()) {
            return Err(format!("duplicate benchmark dataset id {:?}", dataset.id).into());
        }
        for field in [
            &dataset.scientific_name,
            &dataset.source,
            &dataset.study_accession,
            &dataset.run_accession,
            &dataset.sample_accession,
            &dataset.experiment_accession,
            &dataset.library_strategy,
            &dataset.library_layout,
            &dataset.instrument_model,
            &dataset.ploidy,
            &dataset.license,
            &dataset.provenance_url,
            &dataset.notes,
        ] {
            if field.trim().is_empty() {
                return Err(format!("{}: catalog fields must not be empty", dataset.id).into());
            }
        }
        if dataset.read_count == 0 || dataset.base_count == 0 {
            return Err(format!("{}: read_count/base_count must be non-zero", dataset.id).into());
        }
        if !dataset.provenance_url.starts_with("https://") {
            return Err(format!("{}: provenance_url must be https", dataset.id).into());
        }
        if let Some(url) = dataset.ploidy_provenance_url.as_deref() {
            if !url.starts_with("https://") {
                return Err(format!("{}: ploidy_provenance_url must be https", dataset.id).into());
            }
        }

        let mut roles = BTreeSet::new();
        for file in &dataset.files {
            if !roles.insert(file.role.clone()) {
                return Err(format!("{}: duplicate file role {}", dataset.id, file.role).into());
            }
            if !file.url.starts_with("https://") {
                return Err(format!("{}: file URL must be https: {}", dataset.id, file.url).into());
            }
            if let Some(md5) = file.md5.as_deref() {
                if md5.len() != 32 || !md5.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(format!("{}: invalid md5 for {}", dataset.id, file.role).into());
                }
            }
            if file.bytes == 0 {
                return Err(
                    format!("{}: bytes must be non-zero for {}", dataset.id, file.role).into(),
                );
            }
        }
        let mut alignment_roles = BTreeSet::new();
        for alignment in &dataset.alignments {
            if !alignment_roles.insert(alignment.role.clone()) {
                return Err(format!(
                    "{}: duplicate alignment role {}",
                    dataset.id, alignment.role
                )
                .into());
            }
            if !alignment.url.starts_with("https://") {
                return Err(format!(
                    "{}: alignment URL must be https: {}",
                    dataset.id, alignment.url
                )
                .into());
            }
            if let Some(index_url) = alignment.index_url.as_deref() {
                if !index_url.starts_with("https://") {
                    return Err(format!(
                        "{}: alignment index URL must be https: {}",
                        dataset.id, index_url
                    )
                    .into());
                }
            }
            if alignment.bytes == Some(0) {
                return Err(format!(
                    "{}: alignment bytes must be non-zero for {}",
                    dataset.id, alignment.role
                )
                .into());
            }
            if alignment.index_bytes == Some(0) {
                return Err(format!(
                    "{}: alignment index bytes must be non-zero for {}",
                    dataset.id, alignment.role
                )
                .into());
            }
            if let Some(notes) = alignment.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: alignment {} notes must not be empty when present",
                        dataset.id, alignment.role
                    )
                    .into());
                }
            }
        }

        let mut reference_roles = BTreeSet::new();
        for reference in &dataset.references {
            if !reference_roles.insert(reference.role.clone()) {
                return Err(format!(
                    "{}: duplicate reference role {}",
                    dataset.id, reference.role
                )
                .into());
            }
            if !reference.url.starts_with("https://") {
                return Err(format!(
                    "{}: reference URL must be https: {}",
                    dataset.id, reference.url
                )
                .into());
            }
            require_rel_path(root, &dataset.id, "reference.path", &reference.path, false)?;
            if !reference.path.starts_with("data/") {
                return Err(format!(
                    "{}: reference data must live under data/: {}",
                    dataset.id, reference.path
                )
                .into());
            }
            if let Some(expected) = reference.sha256.as_deref() {
                if expected.len() != 64 || !expected.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(format!(
                        "{}: invalid sha256 for reference {}",
                        dataset.id, reference.role
                    )
                    .into());
                }
                let path = root.join(&reference.path);
                if path.exists() {
                    let got = sha256_file(&path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: reference {} digest mismatch: got {}, expected {}",
                            dataset.id, reference.role, got, expected
                        )
                        .into());
                    }
                }
            }
            if let Some(notes) = reference.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: reference {} notes must not be empty when present",
                        dataset.id, reference.role
                    )
                    .into());
                }
            }
        }

        let mut reference_slice_roles = BTreeSet::new();
        for slice in &dataset.reference_slices {
            if !reference_slice_roles.insert(slice.role.clone()) {
                return Err(format!(
                    "{}: duplicate reference slice role {}",
                    dataset.id, slice.role
                )
                .into());
            }
            if !reference_roles.contains(&slice.source_role) {
                return Err(format!(
                    "{}: reference slice {} references unknown reference role {}",
                    dataset.id, slice.role, slice.source_role
                )
                .into());
            }
            parse_region(&slice.region)?;
            require_rel_path(
                root,
                &dataset.id,
                "reference_slice.path",
                &slice.path,
                false,
            )?;
            if !slice.path.starts_with("data/") {
                return Err(format!(
                    "{}: reference slice data must live under data/: {}",
                    dataset.id, slice.path
                )
                .into());
            }
            if let Some(expected) = slice.sha256.as_deref() {
                if expected.len() != 64 || !expected.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(format!(
                        "{}: invalid sha256 for reference slice {}",
                        dataset.id, slice.role
                    )
                    .into());
                }
                let path = root.join(&slice.path);
                if path.exists() {
                    let got = sha256_file(&path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: reference slice {} digest mismatch: got {}, expected {}",
                            dataset.id, slice.role, got, expected
                        )
                        .into());
                    }
                }
            }
            if let Some(notes) = slice.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: reference slice {} notes must not be empty when present",
                        dataset.id, slice.role
                    )
                    .into());
                }
            }
        }

        let mut truth_roles = BTreeSet::new();
        for truth in &dataset.truth {
            if !truth_roles.insert(truth.role.clone()) {
                return Err(format!("{}: duplicate truth role {}", dataset.id, truth.role).into());
            }
            if !truth.url.starts_with("https://") {
                return Err(
                    format!("{}: truth URL must be https: {}", dataset.id, truth.url).into(),
                );
            }
            require_rel_path(root, &dataset.id, "truth.path", &truth.path, false)?;
            if !truth.path.starts_with("data/") {
                return Err(format!(
                    "{}: truth data must live under data/: {}",
                    dataset.id, truth.path
                )
                .into());
            }
            if let Some(md5) = truth.md5.as_deref() {
                if md5.len() != 32 || !md5.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(
                        format!("{}: invalid md5 for truth {}", dataset.id, truth.role).into(),
                    );
                }
            }
            if truth.bytes == Some(0) {
                return Err(format!(
                    "{}: truth bytes must be non-zero for {}",
                    dataset.id, truth.role
                )
                .into());
            }
            if let Some(notes) = truth.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: truth {} notes must not be empty when present",
                        dataset.id, truth.role
                    )
                    .into());
                }
            }
        }

        let mut interval_pair_roles = BTreeSet::new();
        for pair in &dataset.interval_pairs {
            if !interval_pair_roles.insert(pair.role.clone()) {
                return Err(
                    format!("{}: duplicate interval pair role {}", dataset.id, pair.role).into(),
                );
            }
            if !alignment_roles.contains(&pair.source_role) {
                return Err(format!(
                    "{}: interval pair {} references unknown alignment role {}",
                    dataset.id, pair.role, pair.source_role
                )
                .into());
            }
            if pair.region.trim().is_empty()
                || pair.region.contains(char::is_whitespace)
                || !pair.region.contains(':')
                || !pair.region.contains('-')
            {
                return Err(format!(
                    "{}: interval pair {} has invalid region {}",
                    dataset.id, pair.role, pair.region
                )
                .into());
            }
            if pair.records == 0 {
                return Err(format!(
                    "{}: interval pair {} records must be non-zero",
                    dataset.id, pair.role
                )
                .into());
            }
            for (field, path) in [("r1_path", &pair.r1_path), ("r2_path", &pair.r2_path)] {
                require_rel_path(root, &dataset.id, field, path, false)?;
                if !path.starts_with("data/") {
                    return Err(format!(
                        "{}: interval pair data must live under data/: {}",
                        dataset.id, path
                    )
                    .into());
                }
            }
            for (label, sha) in [
                ("r1_sha256", pair.r1_sha256.as_deref()),
                ("r2_sha256", pair.r2_sha256.as_deref()),
            ] {
                if let Some(sha) = sha {
                    if sha.len() != 64 || !sha.chars().all(|c| c.is_ascii_hexdigit()) {
                        return Err(format!(
                            "{}: invalid {} for interval pair {}",
                            dataset.id, label, pair.role
                        )
                        .into());
                    }
                }
            }
            if let Some(notes) = pair.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: interval pair {} notes must not be empty when present",
                        dataset.id, pair.role
                    )
                    .into());
                }
            }
            let r1_path = root.join(&pair.r1_path);
            let r2_path = root.join(&pair.r2_path);
            if r1_path.exists() || r2_path.exists() {
                validate_paired_fastq_record_count(&r1_path, &r2_path, pair.records)?;
                if let Some(expected) = pair.r1_sha256.as_deref() {
                    let got = sha256_file(&r1_path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: interval pair {} r1 digest mismatch: got {}, expected {}",
                            dataset.id, pair.role, got, expected
                        )
                        .into());
                    }
                }
                if let Some(expected) = pair.r2_sha256.as_deref() {
                    let got = sha256_file(&r2_path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: interval pair {} r2 digest mismatch: got {}, expected {}",
                            dataset.id, pair.role, got, expected
                        )
                        .into());
                    }
                }
            }
        }

        let mut local_roles = BTreeSet::new();
        for local in &dataset.local_files {
            if !local_roles.insert(local.role.clone()) {
                return Err(
                    format!("{}: duplicate local file role {}", dataset.id, local.role).into(),
                );
            }
            if local.source_path.trim().is_empty() {
                return Err(format!(
                    "{}: local file {} source_path must not be empty",
                    dataset.id, local.role
                )
                .into());
            }
            require_rel_path(root, &dataset.id, "local_file.path", &local.path, false)?;
            if !local.path.starts_with("data/") {
                return Err(format!(
                    "{}: local file data must live under data/: {}",
                    dataset.id, local.path
                )
                .into());
            }
            if local.records == Some(0) {
                return Err(format!(
                    "{}: local file {} records must be non-zero when present",
                    dataset.id, local.role
                )
                .into());
            }
            if local.bases == Some(0) {
                return Err(format!(
                    "{}: local file {} bases must be non-zero when present",
                    dataset.id, local.role
                )
                .into());
            }
            if let Some(expected) = local.sha256.as_deref() {
                if expected.len() != 64 || !expected.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(format!(
                        "{}: invalid sha256 for local file {}",
                        dataset.id, local.role
                    )
                    .into());
                }
                let path = root.join(&local.path);
                if path.exists() {
                    let got = sha256_file(&path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: local file {} digest mismatch: got {}, expected {}",
                            dataset.id, local.role, got, expected
                        )
                        .into());
                    }
                }
            }
            if let Some(notes) = local.notes.as_deref() {
                if notes.trim().is_empty() {
                    return Err(format!(
                        "{}: local file {} notes must not be empty when present",
                        dataset.id, local.role
                    )
                    .into());
                }
            }
        }

        if roles.is_empty()
            && alignment_roles.is_empty()
            && reference_roles.is_empty()
            && truth_roles.is_empty()
            && local_roles.is_empty()
        {
            return Err(format!("{}: at least one source file is required", dataset.id).into());
        }

        for prepared in &dataset.prepared {
            if !roles.contains(&prepared.source_role) {
                return Err(format!(
                    "{}: prepared {} references unknown source role {}",
                    dataset.id, prepared.role, prepared.source_role
                )
                .into());
            }
            if prepared.records == 0 {
                return Err(format!(
                    "{}: prepared {} records must be non-zero",
                    dataset.id, prepared.role
                )
                .into());
            }
            require_rel_path(root, &dataset.id, "prepared.path", &prepared.path, false)?;
            if !prepared.path.starts_with("data/") {
                return Err(format!(
                    "{}: prepared data must live under data/: {}",
                    dataset.id, prepared.path
                )
                .into());
            }
            if let Some(expected) = prepared.sha256.as_deref() {
                if expected.len() != 64 || !expected.chars().all(|c| c.is_ascii_hexdigit()) {
                    return Err(
                        format!("{}: invalid sha256 for {}", dataset.id, prepared.role).into(),
                    );
                }
                let path = root.join(&prepared.path);
                if path.exists() {
                    let got = sha256_file(&path)?;
                    if got != expected {
                        return Err(format!(
                            "{}: prepared {} digest mismatch: got {}, expected {}",
                            dataset.id, prepared.role, got, expected
                        )
                        .into());
                    }
                }
            }
        }
    }

    println!(
        "xtask validate-data: OK ({} external datasets)",
        catalog.datasets.len()
    );
    Ok(())
}

fn validate_assembler_framework(root: &Path) -> DynResult<()> {
    let framework = load_assembler_framework(root)?;
    if framework.schema_version != 1 {
        return Err("assembler_framework.toml schema_version must be 1".into());
    }
    if framework.modules.is_empty() {
        return Err("assembler_framework.toml must contain at least one [[modules]] entry".into());
    }

    let doc = fs::read_to_string(root.join("docs/ASSEMBLER_FRAMEWORK.md"))?;
    let blueprint = fs::read_to_string(root.join("docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md"))?;
    for stage in PROMOTION_STAGES {
        if !blueprint.contains(stage) {
            return Err(format!(
                "docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md missing promotion stage {stage}"
            )
            .into());
        }
    }

    let mut seen = BTreeSet::new();
    for module in &framework.modules {
        if module.id.trim().is_empty() {
            return Err("assembler framework module id must not be empty".into());
        }
        if !seen.insert(module.id.clone()) {
            return Err(format!("duplicate assembler framework module id {:?}", module.id).into());
        }
        for (field, value) in [
            ("phase", &module.phase),
            ("status", &module.status),
            ("summary", &module.summary),
        ] {
            if value.trim().is_empty() {
                return Err(format!("{}: {field} must not be empty", module.id).into());
            }
        }
        if !matches!(
            module.phase.as_str(),
            "phase1" | "phase2_illumina" | "phase2_deferral"
        ) {
            return Err(format!("{}: unsupported phase {}", module.id, module.phase).into());
        }
        if !matches!(module.status.as_str(), "active" | "planned" | "deferred") {
            return Err(format!("{}: unsupported status {}", module.id, module.status).into());
        }
        validate_framework_list(root, &module.id, "papers", &module.papers, true)?;
        validate_framework_list(
            root,
            &module.id,
            "current_files",
            &module.current_files,
            true,
        )?;
        validate_framework_list(
            root,
            &module.id,
            "promotion_stages",
            &module.promotion_stages,
            false,
        )?;
        for stage in &module.promotion_stages {
            if !PROMOTION_STAGES.contains(&stage.as_str()) {
                return Err(format!("{}: unsupported promotion stage {stage}", module.id).into());
            }
        }
        validate_framework_list(root, &module.id, "gates", &module.gates, false)?;
        validate_framework_list(root, &module.id, "next_work", &module.next_work, false)?;

        let module_words = module.id.replace('_', " ");
        if !doc.contains(&module_words)
            && !doc.contains(&module.id)
            && !blueprint.contains(&module_words)
            && !blueprint.contains(&module.id)
        {
            return Err(format!(
                "framework docs do not mention framework module {}",
                module.id
            )
            .into());
        }
    }

    println!(
        "xtask validate-framework: OK ({} modules)",
        framework.modules.len()
    );
    Ok(())
}

fn validate_development_docs(root: &Path) -> DynResult<()> {
    let development_path = root.join("docs/DEVELOPMENT.md");
    let roadmap_path = root.join("docs/ROADMAP.md");
    let development = fs::read_to_string(&development_path)?;
    let roadmap = fs::read_to_string(&roadmap_path)?;
    let blueprint = fs::read_to_string(root.join("docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md"))?;
    let readme = fs::read_to_string(root.join("README.md"))?;
    let capabilities = fs::read_to_string(root.join("docs/CAPABILITIES.md"))?;

    for phrase in [
        "Orchestrator",
        "Worker Rules",
        "Worker Packet Template",
        "Claim Levels",
        "Acceptance Checklist",
        "Commit Discipline",
        "No `git`.",
        "No `cargo`.",
        "cargo fmt --all --check",
        "cargo clippy --workspace --all-features -- -D warnings",
        "cargo test --workspace --all-features",
        "cargo run -p xtask -- validate",
    ] {
        if !development.contains(phrase) {
            return Err(
                format!("docs/DEVELOPMENT.md missing required protocol phrase: {phrase}").into(),
            );
        }
    }

    for header in [
        "Objective:",
        "Files Read:",
        "Paper/Technique Basis:",
        "Findings:",
        "Proposed Change:",
        "Patch Sketch:",
        "Tests Orchestrator Should Run:",
        "Risks / Unknowns:",
    ] {
        if !development.contains(header) {
            return Err(
                format!("docs/DEVELOPMENT.md missing worker packet header: {header}").into(),
            );
        }
    }

    for level in [
        "`observed`",
        "`tested`",
        "`benchmarked`",
        "`production-gated`",
    ] {
        if !development.contains(level) {
            return Err(format!("docs/DEVELOPMENT.md missing claim level: {level}").into());
        }
    }

    for phrase in [
        "read correction / trusted k-mer model",
        "multi-k graph ladder",
        "repeat-aware graph annotation",
        "decision-first simplification",
        "mate-pair distance/orientation evidence",
        "scaffold/path promotion",
        "polishing/audit loop",
        "diploid ambiguity handling",
        "quality gates and benchmark claims",
    ] {
        if !blueprint.contains(phrase) {
            return Err(format!(
                "docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md missing architecture phrase: {phrase}"
            )
            .into());
        }
    }

    for wave in ["Wave A", "Wave B", "Wave C", "Wave D", "Wave E", "Wave F"] {
        if !roadmap.contains(wave) || !blueprint.contains(wave) {
            return Err(
                format!("framework wave {wave} must appear in roadmap and blueprint").into(),
            );
        }
    }

    for lane in [
        "Quality gates",
        "Evidence ledger",
        "Graph IR",
        "Multi-k selection",
        "Simplification policy",
        "Assembly audit",
        "Diploid semantics",
        "Path/scaffold builder",
        "Benchmark matrix",
        "Literature-derived future adapters",
    ] {
        if !roadmap.contains(lane) {
            return Err(format!("docs/ROADMAP.md missing roadmap lane: {lane}").into());
        }
    }

    for doc in [
        "docs/DEVELOPMENT.md",
        "docs/ROADMAP.md",
        "docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md",
    ] {
        if !readme.contains(doc) {
            return Err(format!("README.md missing link to {doc}").into());
        }
    }
    if !capabilities.contains("validate-development") {
        return Err("docs/CAPABILITIES.md missing validate-development command".into());
    }
    if !capabilities.contains("docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md") {
        return Err("docs/CAPABILITIES.md missing blueprint link".into());
    }

    println!("xtask validate-development: OK");
    Ok(())
}

fn validate_framework_list(
    root: &Path,
    module_id: &str,
    field: &str,
    values: &[String],
    paths_must_exist: bool,
) -> DynResult<()> {
    if values.is_empty() || values.iter().any(|value| value.trim().is_empty()) {
        return Err(format!("{module_id}: {field} must be a non-empty list").into());
    }
    if paths_must_exist {
        for value in values {
            require_rel_path(root, module_id, field, value, true)?;
        }
    }
    Ok(())
}

fn run_bench(root: &Path, tier: Tier, row_filter: Option<&str>, out: &Path) -> DynResult<()> {
    validate_matrix(root)?;
    let matrix = load_matrix(root)?;
    let mut rows = Vec::new();
    let mut failed = false;
    let mut matched_filter = false;

    for row in &matrix.rows {
        let row_id = row.id.clone().ok_or("validated row missing id")?;
        if row_filter.map(|filter| filter != row_id).unwrap_or(false) {
            continue;
        }
        matched_filter = true;
        let ci_tier = row.ci_tier.ok_or("validated row missing ci_tier")?;
        let product_claim = row
            .product_claim
            .clone()
            .ok_or("validated row missing product_claim")?;
        let claim_level = row
            .claim_level
            .clone()
            .ok_or("validated row missing claim_level")?;
        let reference_availability = row
            .reference_availability
            .clone()
            .ok_or("validated row missing reference_availability")?;
        let artifact_policy = row
            .artifact_policy
            .clone()
            .ok_or("validated row missing artifact_policy")?;
        let tools = tool_availability(
            row.required_tools.as_deref().unwrap_or(&[]),
            row.optional_tools.as_deref().unwrap_or(&[]),
        );
        let missing_required_tools = !tools.required_unavailable.is_empty();
        if missing_required_tools {
            failed = true;
        }
        let scripts = scripts_for_tier(row, tier);
        let should_run_trex = row
            .trex
            .as_ref()
            .map(|trex| trex.tiers.contains(&tier))
            .unwrap_or(false);
        if scripts.is_empty() && !should_run_trex {
            continue;
        }

        let mut script_reports = Vec::new();
        if !missing_required_tools {
            for script in scripts {
                let report = run_script(root, script)?;
                if !report.success {
                    failed = true;
                }
                script_reports.push(report);
            }
        }

        let mut trex_reports = Vec::new();
        if !missing_required_tools {
            if let Some(trex) = row.trex.as_ref() {
                if should_run_trex {
                    let report = run_trex_bench(root, trex)?;
                    let quast_failed = report
                        .quast
                        .as_ref()
                        .map(|quast| !quast.success)
                        .unwrap_or(false);
                    if !report.success || quast_failed {
                        failed = true;
                    }
                    trex_reports.push(report);
                }
            }
        }
        let layer_status = row_layer_status(&tools, &script_reports, &trex_reports);

        let artifacts = artifacts_for_tier(row, tier)
            .into_iter()
            .map(|path| {
                let meta = fs::metadata(root.join(path)).ok();
                ArtifactReport {
                    path: path.to_string(),
                    exists: meta.is_some(),
                    bytes: meta.map(|m| m.len()),
                }
            })
            .collect();

        rows.push(BenchRowReport {
            id: row_id,
            ci_tier,
            product_claim,
            claim_level,
            reference_availability,
            artifact_policy,
            tools,
            layer_status,
            scripts: script_reports,
            trex_runs: trex_reports,
            artifacts,
        });
    }
    if let Some(filter) = row_filter {
        if !matched_filter {
            return Err(format!("xtask bench: row {filter:?} not found").into());
        }
    }

    let report = BenchReport {
        schema_version: 1,
        tier,
        generated_unix_ms: SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
        rows,
    };
    let out_path = if out.is_absolute() {
        out.to_path_buf()
    } else {
        root.join(out)
    };
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&out_path, serde_json::to_string_pretty(&report)? + "\n")?;
    println!("xtask bench: wrote {}", out_path.display());

    if failed {
        return Err("xtask bench: one or more rows failed".into());
    }
    Ok(())
}

fn tool_availability(required_tools: &[String], optional_tools: &[String]) -> ToolAvailability {
    let (required_available, required_unavailable) = partition_tools(required_tools);
    let (optional_available, optional_unavailable) = partition_tools(optional_tools);
    ToolAvailability {
        required_available,
        required_unavailable,
        optional_available,
        optional_unavailable,
    }
}

fn partition_tools(tools: &[String]) -> (Vec<String>, Vec<String>) {
    let mut available = Vec::new();
    let mut unavailable = Vec::new();
    for tool in tools {
        if command_available(tool) {
            available.push(tool.clone());
        } else {
            unavailable.push(tool.clone());
        }
    }
    (available, unavailable)
}

fn command_available(name: &str) -> bool {
    let path = Path::new(name);
    if path.components().count() > 1 {
        return path.exists();
    }
    let Some(paths) = std::env::var_os("PATH") else {
        return false;
    };
    std::env::split_paths(&paths).any(|dir| dir.join(name).is_file())
}

fn row_layer_status(
    tools: &ToolAvailability,
    scripts: &[ScriptReport],
    trex_reports: &[TrexRunReport],
) -> RowLayerStatus {
    if let Some(tool) = tools.required_unavailable.first() {
        return RowLayerStatus {
            status: "failed".to_string(),
            failed_layer: Some(format!("required_tool:{tool}")),
        };
    }
    for script in scripts {
        if !script.success {
            return RowLayerStatus {
                status: "failed".to_string(),
                failed_layer: Some(format!("script:{}", script.path)),
            };
        }
    }
    for report in trex_reports {
        if !report.success {
            return RowLayerStatus {
                status: "failed".to_string(),
                failed_layer: Some("trex".to_string()),
            };
        }
        if report
            .quast
            .as_ref()
            .map(|quast| !quast.success)
            .unwrap_or(false)
        {
            return RowLayerStatus {
                status: "failed".to_string(),
                failed_layer: Some("quast".to_string()),
            };
        }
    }
    RowLayerStatus {
        status: "passed".to_string(),
        failed_layer: None,
    }
}

fn run_gate(root: &Path, tier: Tier) -> DynResult<()> {
    validate_matrix(root)?;
    validate_capabilities(root)?;
    run_command_passthrough(
        root,
        "cargo",
        &[
            "clippy",
            "--workspace",
            "--all-features",
            "--",
            "-D",
            "warnings",
        ],
    )?;
    let bench_out = format!("target/benchmarks/{}.json", tier_name(tier));
    run_bench(root, tier, None, Path::new(&bench_out))?;
    if tier == Tier::Pr {
        run_command_passthrough(root, "bash", &["scripts/pr_smoke.sh"])?;
    }
    Ok(())
}

fn required_str<'a>(value: Option<&'a str>, idx: usize, key: &str) -> DynResult<&'a str> {
    match value {
        Some(value) if !value.trim().is_empty() => Ok(value),
        _ => Err(format!("row {}: missing required key {key}", idx + 1).into()),
    }
}

fn required_list<'a>(
    value: Option<&'a Vec<String>>,
    row_id: &str,
    key: &str,
) -> DynResult<&'a [String]> {
    match value {
        Some(value) if !value.is_empty() && value.iter().all(|item| !item.trim().is_empty()) => {
            Ok(value.as_slice())
        }
        _ => Err(format!("{row_id}: {key} must be a non-empty list of strings").into()),
    }
}

fn validate_claim_level(row_id: &str, value: &str) -> DynResult<()> {
    if matches!(
        value,
        "observed" | "tested" | "benchmarked" | "production-gated"
    ) {
        Ok(())
    } else {
        Err(format!("{row_id}: unsupported claim_level {value:?}").into())
    }
}

fn validate_reference_availability(row_id: &str, value: &str) -> DynResult<()> {
    if matches!(
        value,
        "none" | "single_reference" | "parental_haplotypes" | "external_reference"
    ) {
        Ok(())
    } else {
        Err(format!("{row_id}: unsupported reference_availability {value:?}").into())
    }
}

fn validate_artifact_policy(row_id: &str, value: &str) -> DynResult<()> {
    if matches!(value, "required" | "optional" | "manual_archive") {
        Ok(())
    } else {
        Err(format!("{row_id}: unsupported artifact_policy {value:?}").into())
    }
}

fn validate_tool_list(row_id: &str, key: &str, values: Option<&[String]>) -> DynResult<()> {
    for value in values.unwrap_or(&[]) {
        if value.trim().is_empty() {
            return Err(format!("{row_id}: {key} contains an empty entry").into());
        }
        if value.contains('/') || value.contains('\\') {
            return Err(format!("{row_id}: {key} must contain tool names, got {value}").into());
        }
    }
    Ok(())
}

fn require_rel_path(
    root: &Path,
    row_id: &str,
    field: &str,
    value: &str,
    must_exist: bool,
) -> DynResult<()> {
    let path = Path::new(value);
    if path.is_absolute() {
        return Err(format!("{row_id}: {field} must be repo-relative, got {value}").into());
    }
    if path
        .components()
        .any(|component| matches!(component, Component::ParentDir))
    {
        return Err(format!("{row_id}: {field} must not escape repo root: {value}").into());
    }
    if must_exist && !root.join(path).exists() {
        return Err(format!("{row_id}: {field} path does not exist: {value}").into());
    }
    Ok(())
}

fn validate_trex_bench(
    root: &Path,
    row_id: &str,
    trex: &TrexBench,
    is_external: bool,
) -> DynResult<()> {
    if trex.tiers.is_empty() {
        return Err(format!("{row_id}: [rows.trex].tiers must not be empty").into());
    }
    if trex.args.is_empty() || trex.args.iter().any(|arg| arg.trim().is_empty()) {
        return Err(format!("{row_id}: [rows.trex].args must be non-empty strings").into());
    }
    require_rel_path(root, row_id, "trex.out_dir", &trex.out_dir, false)?;
    if let Some(reference) = trex.reference.as_deref() {
        require_rel_path(root, row_id, "trex.reference", reference, !is_external)?;
    }
    let mut parent_labels = BTreeSet::new();
    for parent in &trex.parental_references {
        if parent.label.trim().is_empty() {
            return Err(format!("{row_id}: parental reference label must not be empty").into());
        }
        if !parent_labels.insert(parent.label.clone()) {
            return Err(format!(
                "{row_id}: duplicate parental reference label {:?}",
                parent.label
            )
            .into());
        }
        require_rel_path(
            root,
            row_id,
            "trex.parental_references.path",
            &parent.path,
            !is_external,
        )?;
    }
    if let Some(r1) = value_after(&trex.args, "--r1") {
        require_rel_path(root, row_id, "trex.--r1", r1, !is_external)?;
    }
    if let Some(r2) = value_after(&trex.args, "--r2") {
        require_rel_path(root, row_id, "trex.--r2", r2, !is_external)?;
    }
    match value_after(&trex.args, "--out-dir") {
        Some(out_dir) if out_dir == trex.out_dir => {}
        Some(out_dir) => {
            return Err(format!(
                "{row_id}: [rows.trex].out_dir must match --out-dir, got {out_dir}"
            )
            .into());
        }
        None => return Err(format!("{row_id}: [rows.trex].args must include --out-dir").into()),
    }
    Ok(())
}

fn scripts_for_any_tier(row: &Row) -> Vec<&String> {
    let mut scripts = Vec::new();
    for list in [
        row.pr_scripts.as_ref(),
        row.main_scripts.as_ref(),
        row.nightly_scripts.as_ref(),
        row.manual_scripts.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        scripts.extend(list.iter());
    }
    scripts
}

fn scripts_for_tier(row: &Row, tier: Tier) -> &[String] {
    match tier {
        Tier::Pr => row.pr_scripts.as_deref().unwrap_or(&[]),
        Tier::Main => row.main_scripts.as_deref().unwrap_or(&[]),
        Tier::Nightly => row.nightly_scripts.as_deref().unwrap_or(&[]),
        Tier::Manual => row.manual_scripts.as_deref().unwrap_or(&[]),
    }
}

fn artifacts_for_any_tier(row: &Row) -> Vec<&String> {
    let mut artifacts = Vec::new();
    for list in [
        row.pr_artifacts.as_ref(),
        row.main_artifacts.as_ref(),
        row.nightly_artifacts.as_ref(),
        row.manual_artifacts.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        artifacts.extend(list.iter());
    }
    artifacts
}

fn artifacts_for_tier(row: &Row, tier: Tier) -> Vec<&String> {
    let mut artifacts = Vec::new();
    if let Some(base) = row.artifacts.as_ref() {
        artifacts.extend(base.iter());
    }
    let tier_artifacts = match tier {
        Tier::Pr => row.pr_artifacts.as_ref(),
        Tier::Main => row.main_artifacts.as_ref(),
        Tier::Nightly => row.nightly_artifacts.as_ref(),
        Tier::Manual => row.manual_artifacts.as_ref(),
    };
    if let Some(tier_artifacts) = tier_artifacts {
        artifacts.extend(tier_artifacts.iter());
    }
    artifacts
}

fn run_trex_bench(root: &Path, trex: &TrexBench) -> DynResult<TrexRunReport> {
    run_command_passthrough(
        root,
        "cargo",
        &["build", "-q", "--release", "-p", "trex-cli"],
    )?;
    let mut command = vec!["target/release/trex".to_string()];
    command.extend(trex.args.clone());

    println!("xtask bench: running {}", command.join(" "));
    let start = std::time::Instant::now();
    let output = if Path::new("/usr/bin/time").exists() {
        let mut cmd = Command::new("/usr/bin/time");
        cmd.current_dir(root)
            .args(["-f", "TREX_XTASK_TIME\t%e\t%M"])
            .args(&command)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    } else {
        Command::new(root.join("target/release/trex"))
            .current_dir(root)
            .args(&command[1..])
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    };
    let wall_ms = start.elapsed().as_millis();
    let stdout = String::from_utf8_lossy(&output.stdout);
    let observed = parse_trex_observed(&stdout);
    if !stdout.is_empty() {
        print!("{stdout}");
    }
    let stderr = String::from_utf8_lossy(&output.stderr);
    let (time_seconds, max_rss_kib, passthrough_stderr) = parse_time_stderr(&stderr);
    if !passthrough_stderr.trim().is_empty() {
        eprint!("{passthrough_stderr}");
        if !passthrough_stderr.ends_with('\n') {
            eprintln!();
        }
    }
    let metrics = if output.status.success() {
        Some(assembly_metrics(root, trex)?)
    } else {
        None
    };
    let quast = if output.status.success() {
        run_optional_quast(root, trex)?
    } else {
        None
    };
    Ok(TrexRunReport {
        command,
        exit_code: output.status.code(),
        success: output.status.success(),
        wall_ms,
        time_seconds,
        max_rss_kib,
        observed,
        metrics,
        quast,
    })
}

fn assembly_metrics(root: &Path, trex: &TrexBench) -> DynResult<AssemblyMetrics> {
    let out_dir = root.join(&trex.out_dir);
    let kmer_size =
        value_after(&trex.args, "--kmer-size").and_then(|value| value.parse::<usize>().ok());
    let reads = value_after(&trex.args, "--r1")
        .map(|path| fastq_stats(&root.join(path), kmer_size))
        .transpose()?;
    let r2_reads = value_after(&trex.args, "--r2")
        .map(|path| fastq_stats(&root.join(path), kmer_size))
        .transpose()?;
    let candidate_kmers = match (
        reads.as_ref().and_then(|stats| stats.candidate_kmers),
        r2_reads.as_ref().and_then(|stats| stats.candidate_kmers),
    ) {
        (Some(r1), Some(r2)) => Some(r1 + r2),
        (Some(r1), None) => Some(r1),
        (None, Some(r2)) => Some(r2),
        (None, None) => None,
    };
    let reference = trex
        .reference
        .as_deref()
        .map(|path| fasta_stats(&root.join(path)))
        .transpose()?;
    let contigs = fasta_stats_optional(&out_dir.join("contigs.fa"))?;
    let unitigs = fasta_stats_optional(&out_dir.join("unitigs.fa"))?;
    let gfa = gfa_stats_optional(&out_dir.join("graph.gfa"))?;
    let evidence = json_summary_optional(&out_dir.join("evidence.json"))?;
    let annotations = json_summary_optional(&out_dir.join("annotations.json"))?;
    let simplification = json_summary_optional(&out_dir.join("simplification.json"))?;
    let scaffolds = json_summary_optional(&out_dir.join("scaffolds.json"))?;
    let multi_k = json_summary_optional(&out_dir.join("multi_k.json"))?;
    let fragmentation = json_summary_optional(&out_dir.join("fragmentation.json"))?;
    let audit = json_summary_optional(&out_dir.join("audit.json"))?;
    let diploid = json_summary_optional(&out_dir.join("diploid.json"))?;
    let reference_quality = match (trex.reference.as_deref(), kmer_size) {
        (Some(reference), Some(k)) => {
            let contigs_path = out_dir.join("contigs.fa");
            if contigs_path.exists() {
                Some(reference_quality(&root.join(reference), &contigs_path, k)?)
            } else {
                None
            }
        }
        _ => None,
    };
    let parental_reference_quality = match kmer_size {
        Some(k) => parental_reference_quality(root, trex, &out_dir.join("contigs.fa"), k)?,
        None => Vec::new(),
    };
    let read_assembly_quality = match kmer_size {
        Some(k) => {
            let contigs_path = out_dir.join("contigs.fa");
            if contigs_path.exists() {
                let reliable_threshold = value_after(&trex.args, "--trusted-threshold")
                    .and_then(|value| value.parse::<usize>().ok())
                    .unwrap_or(2);
                Some(read_assembly_kmer_quality(
                    root,
                    &trex.args,
                    &contigs_path,
                    k,
                    reliable_threshold,
                )?)
            } else {
                None
            }
        }
        None => None,
    };
    Ok(AssemblyMetrics {
        kmer_size,
        candidate_kmers,
        reads,
        r2_reads,
        reference,
        contigs,
        unitigs,
        gfa,
        evidence,
        annotations,
        simplification,
        scaffolds,
        multi_k,
        fragmentation,
        audit,
        diploid,
        reference_quality,
        parental_reference_quality,
        read_assembly_quality,
    })
}

fn json_summary_optional(path: &Path) -> DynResult<Option<Value>> {
    if !path.exists() {
        return Ok(None);
    }
    let bytes = fs::read(path)?;
    let value: Value = serde_json::from_slice(&bytes)?;
    Ok(Some(summarize_json_sidecar(value)))
}

fn summarize_json_sidecar(value: Value) -> Value {
    let Value::Object(object) = value else {
        return value;
    };
    let mut summary = Map::new();
    for key in ["schema_version", "summary"] {
        if let Some(value) = object.get(key) {
            summary.insert(key.to_string(), compact_json_value(value));
        }
    }
    for (key, value) in object {
        if summary.contains_key(&key) {
            continue;
        }
        match value {
            Value::Array(values) => {
                summary.insert(format!("{key}_count"), Value::from(values.len()));
            }
            Value::Object(values) => {
                summary.insert(format!("{key}_count"), Value::from(values.len()));
            }
            value => {
                summary.insert(key, value);
            }
        }
    }
    Value::Object(summary)
}

fn compact_json_value(value: &Value) -> Value {
    match value {
        Value::Array(values) => Value::from(values.len()),
        Value::Object(values) => {
            let mut compacted = Map::new();
            for (key, value) in values {
                match value {
                    Value::Array(values) => {
                        compacted.insert(format!("{key}_count"), Value::from(values.len()));
                    }
                    Value::Object(_) => {
                        compacted.insert(key.clone(), compact_json_value(value));
                    }
                    value => {
                        compacted.insert(key.clone(), value.clone());
                    }
                }
            }
            Value::Object(compacted)
        }
        value => value.clone(),
    }
}

fn run_optional_quast(root: &Path, trex: &TrexBench) -> DynResult<Option<ScriptReport>> {
    if std::env::var("TREX_RUN_QUAST").ok().as_deref() != Some("1") {
        return Ok(None);
    }
    let Some(reference) = trex.reference.as_deref() else {
        return Ok(None);
    };
    let asm = root.join(&trex.out_dir).join("contigs.fa");
    if !asm.exists() {
        return Ok(None);
    }
    let quast_out = root.join(&trex.out_dir).with_file_name(format!(
        "{}-quast",
        Path::new(&trex.out_dir)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("trex")
    ));
    let script = "scripts/reference_quast.sh";
    println!("xtask bench: running optional QUAST for {}", trex.out_dir);
    let start = std::time::Instant::now();
    let output = if Path::new("/usr/bin/time").exists() {
        Command::new("/usr/bin/time")
            .current_dir(root)
            .env("TREX_QUAST_REF", root.join(reference))
            .env("TREX_QUAST_ASM", asm)
            .env("TREX_QUAST_OUT", &quast_out)
            .args(["-f", "TREX_XTASK_TIME\t%e\t%M", "bash", script])
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    } else {
        Command::new("bash")
            .current_dir(root)
            .env("TREX_QUAST_REF", root.join(reference))
            .env("TREX_QUAST_ASM", asm)
            .env("TREX_QUAST_OUT", &quast_out)
            .arg(script)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    };
    let wall_ms = start.elapsed().as_millis();
    let stdout = String::from_utf8_lossy(&output.stdout);
    if !stdout.is_empty() {
        print!("{stdout}");
    }
    let stderr = String::from_utf8_lossy(&output.stderr);
    let (time_seconds, max_rss_kib, passthrough_stderr) = parse_time_stderr(&stderr);
    if !passthrough_stderr.trim().is_empty() {
        eprint!("{passthrough_stderr}");
        if !passthrough_stderr.ends_with('\n') {
            eprintln!();
        }
    }
    Ok(Some(ScriptReport {
        path: script.to_string(),
        exit_code: output.status.code(),
        success: output.status.success(),
        wall_ms,
        time_seconds,
        max_rss_kib,
    }))
}

fn parse_trex_observed(stdout: &str) -> Option<TrexObserved> {
    let mut observed = TrexObserved {
        reads: None,
        unique_kmers: None,
        trusted_kmers: None,
        unitigs: None,
        contigs: None,
    };
    for token in stdout.split_whitespace() {
        let Some((key, value)) = token.split_once('=') else {
            continue;
        };
        let value = value.trim_end_matches(',');
        let parsed = value.parse::<usize>().ok();
        match key {
            "reads" => observed.reads = parsed,
            "unique_kmers" => observed.unique_kmers = parsed,
            "trusted_kmers" => observed.trusted_kmers = parsed,
            "unitigs" => observed.unitigs = parsed,
            "contigs" => observed.contigs = parsed,
            _ => {}
        }
    }
    if observed.reads.is_some()
        || observed.unique_kmers.is_some()
        || observed.trusted_kmers.is_some()
        || observed.unitigs.is_some()
        || observed.contigs.is_some()
    {
        Some(observed)
    } else {
        None
    }
}

fn fasta_stats_optional(path: &Path) -> DynResult<Option<FastaStats>> {
    if path.exists() {
        Ok(Some(fasta_stats(path)?))
    } else {
        Ok(None)
    }
}

fn gfa_stats_optional(path: &Path) -> DynResult<Option<GfaStats>> {
    if !path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(path)?;
    let mut s_lines = 0;
    let mut l_lines = 0;
    let mut p_lines = 0;
    for line in text.lines() {
        if line.starts_with("S\t") {
            s_lines += 1;
        } else if line.starts_with("L\t") {
            l_lines += 1;
        } else if line.starts_with("P\t") {
            p_lines += 1;
        }
    }
    Ok(Some(GfaStats {
        s_lines,
        l_lines,
        p_lines,
        bytes: fs::metadata(path)?.len(),
    }))
}

fn fasta_stats(path: &Path) -> DynResult<FastaStats> {
    let seqs = fasta_sequences(path)?;
    let lens: Vec<usize> = seqs.iter().map(Vec::len).collect();
    Ok(length_stats(lens))
}

fn fastq_stats(path: &Path, kmer_size: Option<usize>) -> DynResult<FastqStats> {
    let seqs = fastq_sequences(path)?;
    let records = seqs.len();
    let total_bases: usize = seqs.iter().map(Vec::len).sum();
    let min_len = seqs.iter().map(Vec::len).min().unwrap_or(0);
    let max_len = seqs.iter().map(Vec::len).max().unwrap_or(0);
    let candidate_kmers = kmer_size.map(|k| {
        seqs.iter()
            .filter(|seq| seq.len() >= k)
            .map(|seq| seq.len() - k + 1)
            .sum()
    });
    Ok(FastqStats {
        records,
        total_bases,
        min_len,
        max_len,
        candidate_kmers,
    })
}

fn fastq_sequences(path: &Path) -> DynResult<Vec<Vec<u8>>> {
    let text = fs::read_to_string(path)?;
    let mut seqs = Vec::new();
    let mut lines = text.lines();
    while let Some(header) = lines.next() {
        if header.is_empty() {
            continue;
        }
        let seq = lines.next().ok_or("truncated FASTQ sequence line")?;
        let plus = lines.next().ok_or("truncated FASTQ plus line")?;
        let qual = lines.next().ok_or("truncated FASTQ quality line")?;
        if !header.starts_with('@') || !plus.starts_with('+') || qual.len() != seq.len() {
            return Err(format!("invalid FASTQ record in {}", path.display()).into());
        }
        seqs.push(seq.as_bytes().iter().map(u8::to_ascii_uppercase).collect());
    }
    Ok(seqs)
}

fn length_stats(mut lens: Vec<usize>) -> FastaStats {
    if lens.is_empty() {
        return FastaStats {
            records: 0,
            total_bases: 0,
            min_len: 0,
            max_len: 0,
            n50: 0,
        };
    }
    let records = lens.len();
    let total_bases: usize = lens.iter().sum();
    let min_len = *lens.iter().min().unwrap_or(&0);
    let max_len = *lens.iter().max().unwrap_or(&0);
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let mut acc = 0;
    let mut n50 = 0;
    for len in lens {
        acc += len;
        if acc * 2 >= total_bases {
            n50 = len;
            break;
        }
    }
    FastaStats {
        records,
        total_bases,
        min_len,
        max_len,
        n50,
    }
}

fn fasta_sequences(path: &Path) -> DynResult<Vec<Vec<u8>>> {
    let text = fs::read_to_string(path)?;
    let mut seqs = Vec::new();
    let mut cur = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !cur.is_empty() {
                seqs.push(cur);
                cur = Vec::new();
            }
        } else {
            cur.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    Ok(seqs)
}

fn fasta_sequence_for_contig(path: &Path, target: &str) -> DynResult<Vec<u8>> {
    let text = fs::read_to_string(path)?;
    let mut current_name: Option<String> = None;
    let mut current_sequence = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if current_name.as_deref() == Some(target) {
                return Ok(current_sequence);
            }
            current_name = header.split_whitespace().next().map(str::to_string);
            current_sequence = Vec::new();
        } else if current_name.is_some() {
            current_sequence.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if current_name.as_deref() == Some(target) {
        return Ok(current_sequence);
    }
    Err(format!("{}: FASTA contig {target:?} not found", path.display()).into())
}

fn parse_region(value: &str) -> DynResult<GenomicRegion> {
    let (contig, rest) = value
        .split_once(':')
        .ok_or_else(|| format!("invalid genomic region {value:?}: missing ':'"))?;
    let (start, end) = rest
        .split_once('-')
        .ok_or_else(|| format!("invalid genomic region {value:?}: missing '-'"))?;
    if contig.trim().is_empty() {
        return Err(format!("invalid genomic region {value:?}: empty contig").into());
    }
    let start = start
        .replace(',', "")
        .parse::<usize>()
        .map_err(|_| format!("invalid genomic region {value:?}: bad start"))?;
    let end = end
        .replace(',', "")
        .parse::<usize>()
        .map_err(|_| format!("invalid genomic region {value:?}: bad end"))?;
    if start == 0 || end < start {
        return Err(
            format!("invalid genomic region {value:?}: expected 1-based start <= end").into(),
        );
    }
    Ok(GenomicRegion {
        contig: contig.to_string(),
        start,
        end,
    })
}

fn write_wrapped_fasta(path: &Path, name: &str, sequence: &[u8]) -> DynResult<()> {
    let mut out = String::new();
    out.push('>');
    out.push_str(name);
    out.push('\n');
    for chunk in sequence.chunks(80) {
        out.push_str(std::str::from_utf8(chunk)?);
        out.push('\n');
    }
    fs::write(path, out)?;
    Ok(())
}

fn reference_quality(
    reference_path: &Path,
    contigs_path: &Path,
    k: usize,
) -> DynResult<ReferenceQuality> {
    let reference = fasta_sequences(reference_path)?;
    let contigs = fasta_sequences(contigs_path)?;
    Ok(reference_quality_from_sequences(&reference, &contigs, k))
}

fn parental_reference_quality(
    root: &Path,
    trex: &TrexBench,
    contigs_path: &Path,
    k: usize,
) -> DynResult<Vec<NamedReferenceQuality>> {
    if trex.parental_references.is_empty() || !contigs_path.exists() {
        return Ok(Vec::new());
    }
    let contigs = fasta_sequences(contigs_path)?;
    let mut out = Vec::with_capacity(trex.parental_references.len());
    for parent in &trex.parental_references {
        let path = root.join(&parent.path);
        let reference = fasta_sequences(&path)?;
        let reference_stats = length_stats(reference.iter().map(Vec::len).collect());
        let quality = reference_quality_from_sequences(&reference, &contigs, k);
        out.push(NamedReferenceQuality {
            label: parent.label.clone(),
            path: parent.path.clone(),
            reference: reference_stats,
            quality,
        });
    }
    Ok(out)
}

fn reference_quality_from_sequences(
    reference: &[Vec<u8>],
    contigs: &[Vec<u8>],
    k: usize,
) -> ReferenceQuality {
    let reference_kmers = canonical_kmer_set(reference, k);
    let mut contig_kmers = 0usize;
    let mut contig_kmers_in_reference = 0usize;
    let mut contigs_with_reference_hit = 0usize;
    let mut contig_bases = 0usize;
    let mut contig_bases_with_reference_hit = 0usize;

    for contig in contigs {
        contig_bases += contig.len();
        let mut has_hit = false;
        for window in acgt_windows(contig, k) {
            contig_kmers += 1;
            if reference_kmers.contains(&canonical_kmer(window)) {
                contig_kmers_in_reference += 1;
                has_hit = true;
            }
        }
        if has_hit {
            contigs_with_reference_hit += 1;
            contig_bases_with_reference_hit += contig.len();
        }
    }

    let contig_kmer_reference_fraction = if contig_kmers == 0 {
        0.0
    } else {
        contig_kmers_in_reference as f64 / contig_kmers as f64
    };

    ReferenceQuality {
        kmer_size: k,
        contig_kmers,
        contig_kmers_in_reference,
        contig_kmer_reference_fraction,
        contigs_with_reference_hit,
        contig_records: contigs.len(),
        contig_bases,
        contig_bases_with_reference_hit,
    }
}

fn read_assembly_kmer_quality(
    root: &Path,
    args: &[String],
    contigs_path: &Path,
    k: usize,
    reliable_threshold: usize,
) -> DynResult<ReadAssemblyKmerQuality> {
    let r1 = value_after(args, "--r1").ok_or("Trex benchmark args missing --r1")?;
    let mut reads = sequence_records(&root.join(r1))?;
    if let Some(r2) = value_after(args, "--r2") {
        reads.extend(sequence_records(&root.join(r2))?);
    }
    let assembly = fasta_sequences(contigs_path)?;
    Ok(read_assembly_kmer_quality_from_sequences(
        &reads,
        &assembly,
        k,
        reliable_threshold,
    ))
}

fn sequence_records(path: &Path) -> DynResult<Vec<Vec<u8>>> {
    if path_is_fasta(path) {
        fasta_sequences(path)
    } else {
        fastq_sequences(path)
    }
}

fn path_is_fasta(path: &Path) -> bool {
    let s = path.to_string_lossy().to_ascii_lowercase();
    s.ends_with(".fa")
        || s.ends_with(".fasta")
        || s.ends_with(".fna")
        || s.ends_with(".fa.gz")
        || s.ends_with(".fasta.gz")
        || s.ends_with(".fna.gz")
}

fn read_assembly_kmer_quality_from_sequences(
    reads: &[Vec<u8>],
    assembly: &[Vec<u8>],
    k: usize,
    reliable_threshold: usize,
) -> ReadAssemblyKmerQuality {
    let read_counts = canonical_kmer_counts(reads, k);
    let assembly_counts = canonical_kmer_counts(assembly, k);
    let reliable_read_kmers = read_counts
        .values()
        .filter(|count| **count >= reliable_threshold)
        .count();
    let reliable_read_kmers_in_assembly = read_counts
        .iter()
        .filter(|(kmer, count)| {
            **count >= reliable_threshold && assembly_counts.contains_key(*kmer)
        })
        .count();
    let assembly_only_kmers = assembly_counts
        .keys()
        .filter(|kmer| {
            read_counts
                .get(*kmer)
                .map(|count| *count < reliable_threshold)
                .unwrap_or(true)
        })
        .count();
    let assembly_total_kmers: usize = assembly_counts.values().sum();
    let assembly_only_total_kmers: usize = assembly_counts
        .iter()
        .filter(|(kmer, _count)| {
            read_counts
                .get(*kmer)
                .map(|count| *count < reliable_threshold)
                .unwrap_or(true)
        })
        .map(|(_kmer, count)| *count)
        .sum();
    let reliable_read_containment_fraction =
        fraction(reliable_read_kmers_in_assembly, reliable_read_kmers);
    let assembly_only_fraction = fraction(assembly_only_kmers, assembly_counts.len());
    let assembly_only_total_fraction = fraction(assembly_only_total_kmers, assembly_total_kmers);
    let approximate_qv = if assembly_only_fraction > 0.0 {
        Some(-10.0 * assembly_only_fraction.log10())
    } else {
        None
    };

    ReadAssemblyKmerQuality {
        kmer_size: k,
        reliable_threshold,
        read_distinct_kmers: read_counts.len(),
        reliable_read_kmers,
        assembly_distinct_kmers: assembly_counts.len(),
        assembly_total_kmers,
        reliable_read_kmers_in_assembly,
        reliable_read_containment_fraction,
        assembly_only_kmers,
        assembly_only_total_kmers,
        assembly_only_fraction,
        assembly_only_total_fraction,
        qv_error_rate: assembly_only_fraction,
        approximate_qv,
    }
}

fn canonical_kmer_counts(seqs: &[Vec<u8>], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut out = HashMap::new();
    if k == 0 {
        return out;
    }
    for seq in seqs {
        for window in acgt_windows(seq, k) {
            *out.entry(canonical_kmer(window)).or_insert(0) += 1;
        }
    }
    out
}

fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn canonical_kmer_set(seqs: &[Vec<u8>], k: usize) -> HashSet<Vec<u8>> {
    let mut out = HashSet::new();
    if k == 0 {
        return out;
    }
    for seq in seqs {
        for window in acgt_windows(seq, k) {
            out.insert(canonical_kmer(window));
        }
    }
    out
}

fn acgt_windows(seq: &[u8], k: usize) -> impl Iterator<Item = &[u8]> {
    seq.windows(k).filter(|window| {
        window
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    })
}

fn canonical_kmer(window: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(window);
    if window <= rc.as_slice() {
        window.to_vec()
    } else {
        rc
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn fetch_data(root: &Path, dataset_filter: Option<&str>) -> DynResult<()> {
    validate_data_catalog(root)?;
    let catalog = load_data_catalog(root)?;
    let mut matched = 0usize;
    for dataset in &catalog.datasets {
        if let Some(filter) = dataset_filter {
            if dataset.id != filter {
                continue;
            }
        }
        matched += 1;
        fetch_dataset(root, dataset)?;
    }
    if matched == 0 {
        return Err(format!(
            "xtask fetch-data: no dataset matched {}",
            dataset_filter.unwrap_or("<all>")
        )
        .into());
    }
    Ok(())
}

fn fetch_dataset(root: &Path, dataset: &DataSet) -> DynResult<()> {
    for local in &dataset.local_files {
        link_local_data_file(root, dataset, local)?;
    }
    for reference in &dataset.references {
        fetch_reference(root, dataset, reference)?;
    }
    let references: BTreeMap<&str, &DataReference> = dataset
        .references
        .iter()
        .map(|reference| (reference.role.as_str(), reference))
        .collect();
    for slice in &dataset.reference_slices {
        let source = references.get(slice.source_role.as_str()).ok_or_else(|| {
            format!(
                "{}: reference slice {} references unknown reference role {}",
                dataset.id, slice.role, slice.source_role
            )
        })?;
        fetch_reference_slice(root, dataset, source, slice)?;
    }
    for truth in &dataset.truth {
        fetch_truth_file(root, dataset, truth)?;
    }

    let alignments: BTreeMap<&str, &DataAlignmentFile> = dataset
        .alignments
        .iter()
        .map(|alignment| (alignment.role.as_str(), alignment))
        .collect();
    for pair in &dataset.interval_pairs {
        let source = alignments.get(pair.source_role.as_str()).ok_or_else(|| {
            format!(
                "{}: interval pair {} references unknown alignment role {}",
                dataset.id, pair.role, pair.source_role
            )
        })?;
        fetch_interval_pair(root, dataset, source, pair)?;
    }

    let files: BTreeMap<&str, &DataFile> = dataset
        .files
        .iter()
        .map(|file| (file.role.as_str(), file))
        .collect();
    for prepared in &dataset.prepared {
        let source = files.get(prepared.source_role.as_str()).ok_or_else(|| {
            format!(
                "{}: prepared {} references unknown source role {}",
                dataset.id, prepared.role, prepared.source_role
            )
        })?;
        let out_path = root.join(&prepared.path);
        let expected = prepared.sha256.as_deref();
        if out_path.exists() {
            let got = sha256_file(&out_path)?;
            if expected.map(|sha| sha == got).unwrap_or(false) {
                println!(
                    "xtask fetch-data: {} {} already present ({got})",
                    dataset.id, prepared.role
                );
                continue;
            }
        }
        stream_fastq_prefix(source, &out_path, prepared.records)?;
        validate_fastq_record_count(&out_path, prepared.records)?;
        let got = sha256_file(&out_path)?;
        match expected {
            Some(expected) if expected == got => {
                println!(
                    "xtask fetch-data: {} {} OK ({got})",
                    dataset.id, prepared.role
                );
            }
            Some(expected) => {
                return Err(format!(
                    "xtask fetch-data: {} {} digest mismatch: got {}, expected {}",
                    dataset.id, prepared.role, got, expected
                )
                .into());
            }
            None => {
                println!(
                    "xtask fetch-data: {} {} wrote {} records sha256={got}",
                    dataset.id, prepared.role, prepared.records
                );
            }
        }
    }
    Ok(())
}

fn link_local_data_file(root: &Path, dataset: &DataSet, local: &LocalDataFile) -> DynResult<()> {
    let source_path = Path::new(&local.source_path);
    if !source_path.exists() {
        return Err(format!(
            "xtask fetch-data: {} local file {} source missing: {}",
            dataset.id, local.role, local.source_path
        )
        .into());
    }

    let out_path = root.join(&local.path);
    let expected = local.sha256.as_deref();
    if out_path.exists() {
        if let Some(records) = local.records {
            validate_fastq_record_count(&out_path, records)?;
        }
        let got = sha256_file(&out_path)?;
        if expected.map(|sha| sha == got).unwrap_or(false) {
            println!(
                "xtask fetch-data: {} local file {} already linked ({got})",
                dataset.id, local.role
            );
            return Ok(());
        }
        return Err(format!(
            "xtask fetch-data: {} local file {} exists with unexpected digest: got {}, expected {}",
            dataset.id,
            local.role,
            got,
            expected.unwrap_or("<unrecorded>")
        )
        .into());
    }

    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    symlink_or_copy_file(source_path, &out_path)?;
    if let Some(records) = local.records {
        validate_fastq_record_count(&out_path, records)?;
    }
    let got = sha256_file(&out_path)?;
    match expected {
        Some(expected) if expected == got => {
            println!(
                "xtask fetch-data: {} local file {} linked ({got})",
                dataset.id, local.role
            );
        }
        Some(expected) => {
            let _ = fs::remove_file(&out_path);
            return Err(format!(
                "xtask fetch-data: {} local file {} digest mismatch: got {}, expected {}",
                dataset.id, local.role, got, expected
            )
            .into());
        }
        None => {
            println!(
                "xtask fetch-data: {} local file {} linked sha256={got}",
                dataset.id, local.role
            );
        }
    }
    Ok(())
}

fn symlink_or_copy_file(source: &Path, target: &Path) -> DynResult<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(source, target)?;
        Ok(())
    }
    #[cfg(not(unix))]
    {
        fs::copy(source, target)?;
        Ok(())
    }
}

fn fetch_reference_slice(
    root: &Path,
    dataset: &DataSet,
    source: &DataReference,
    slice: &DataReferenceSlice,
) -> DynResult<()> {
    let out_path = root.join(&slice.path);
    let expected = slice.sha256.as_deref();
    if out_path.exists() {
        let got = sha256_file(&out_path)?;
        if expected.map(|sha| sha == got).unwrap_or(false) {
            println!(
                "xtask fetch-data: {} reference slice {} already present ({got})",
                dataset.id, slice.role
            );
            return Ok(());
        }
    }

    let region = parse_region(&slice.region)?;
    let source_path = root.join(&source.path);
    let sequence = fasta_sequence_for_contig(&source_path, &region.contig)?;
    if region.end > sequence.len() {
        return Err(format!(
            "{}: reference slice {} end {} exceeds {} length {}",
            dataset.id,
            slice.role,
            region.end,
            region.contig,
            sequence.len()
        )
        .into());
    }
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = out_path.with_extension("tmp");
    let subseq = &sequence[(region.start - 1)..region.end];
    write_wrapped_fasta(&tmp, &slice.region, subseq)?;
    let got = sha256_file(&tmp)?;
    match expected {
        Some(expected) if expected == got => {}
        Some(expected) => {
            let _ = fs::remove_file(&tmp);
            return Err(format!(
                "xtask fetch-data: {} reference slice {} digest mismatch: got {}, expected {}",
                dataset.id, slice.role, got, expected
            )
            .into());
        }
        None => {}
    }
    fs::rename(tmp, &out_path)?;
    println!(
        "xtask fetch-data: {} reference slice {} wrote {} bp sha256={got}",
        dataset.id,
        slice.role,
        subseq.len()
    );
    Ok(())
}

fn fetch_interval_pair(
    root: &Path,
    dataset: &DataSet,
    source: &DataAlignmentFile,
    pair: &IntervalPairDataFile,
) -> DynResult<()> {
    let r1_path = root.join(&pair.r1_path);
    let r2_path = root.join(&pair.r2_path);
    if r1_path.exists() && r2_path.exists() {
        validate_paired_fastq_record_count(&r1_path, &r2_path, pair.records)?;
        let r1 = sha256_file(&r1_path)?;
        let r2 = sha256_file(&r2_path)?;
        let r1_ok = pair
            .r1_sha256
            .as_deref()
            .map(|expected| expected == r1)
            .unwrap_or(false);
        let r2_ok = pair
            .r2_sha256
            .as_deref()
            .map(|expected| expected == r2)
            .unwrap_or(false);
        if r1_ok && r2_ok {
            println!(
                "xtask fetch-data: {} {} already present r1_sha256={} r2_sha256={}",
                dataset.id, pair.role, r1, r2
            );
            return Ok(());
        }
    }

    if let Some(parent) = r1_path.parent() {
        fs::create_dir_all(parent)?;
    }
    if let Some(parent) = r2_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let scratch = root
        .join("target/benchmarks/fetch-tmp")
        .join(&dataset.id)
        .join(&pair.role);
    fs::create_dir_all(&scratch)?;
    let r1_all = scratch.join("r1.all.fq");
    let r2_all = scratch.join("r2.all.fq");
    let singleton = scratch.join("singletons.fq");
    let other = scratch.join("other.fq");
    for path in [&r1_all, &r2_all, &singleton, &other] {
        let _ = fs::remove_file(path);
    }

    let script = format!(
        "samtools view -u -f 3 -F 2304 '{}' '{}' | samtools collate -u -O - | samtools fastq -n -1 '{}' -2 '{}' -0 '{}' -s '{}' -",
        source.url,
        pair.region,
        r1_all.display(),
        r2_all.display(),
        other.display(),
        singleton.display()
    );
    let status = Command::new("bash").arg("-c").arg(&script).status()?;
    if !status.success() {
        return Err(format!(
            "xtask fetch-data: {} {} interval extraction failed",
            dataset.id, pair.role
        )
        .into());
    }

    let r1_tmp = r1_path.with_extension("tmp");
    let r2_tmp = r2_path.with_extension("tmp");
    write_fastq_prefix(&r1_all, &r1_tmp, pair.records)?;
    write_fastq_prefix(&r2_all, &r2_tmp, pair.records)?;
    validate_paired_fastq_record_count(&r1_tmp, &r2_tmp, pair.records)?;
    let r1 = sha256_file(&r1_tmp)?;
    let r2 = sha256_file(&r2_tmp)?;
    match (pair.r1_sha256.as_deref(), pair.r2_sha256.as_deref()) {
        (Some(expected_r1), Some(expected_r2)) if expected_r1 == r1 && expected_r2 == r2 => {}
        (Some(expected_r1), Some(expected_r2)) => {
            let _ = fs::remove_file(&r1_tmp);
            let _ = fs::remove_file(&r2_tmp);
            return Err(format!(
                "xtask fetch-data: {} {} digest mismatch: got r1={} r2={}, expected r1={} r2={}",
                dataset.id, pair.role, r1, r2, expected_r1, expected_r2
            )
            .into());
        }
        _ => {}
    }
    fs::rename(r1_tmp, &r1_path)?;
    fs::rename(r2_tmp, &r2_path)?;
    println!(
        "xtask fetch-data: {} {} wrote {} read pairs from {} r1_sha256={} r2_sha256={}",
        dataset.id, pair.role, pair.records, pair.region, r1, r2
    );
    Ok(())
}

fn fetch_truth_file(root: &Path, dataset: &DataSet, truth: &DataTruthFile) -> DynResult<()> {
    let out_path = root.join(&truth.path);
    if out_path.exists() {
        if let Some(bytes) = truth.bytes {
            let got = fs::metadata(&out_path)?.len();
            if got == bytes {
                println!(
                    "xtask fetch-data: {} truth {} already present ({} bytes)",
                    dataset.id, truth.role, got
                );
                return Ok(());
            }
        } else {
            println!(
                "xtask fetch-data: {} truth {} already present",
                dataset.id, truth.role
            );
            return Ok(());
        }
    }
    stream_opaque_file(&truth.url, &out_path)?;
    let got = fs::metadata(&out_path)?.len();
    if let Some(expected) = truth.bytes {
        if got != expected {
            return Err(format!(
                "xtask fetch-data: {} truth {} byte mismatch: got {}, expected {}",
                dataset.id, truth.role, got, expected
            )
            .into());
        }
    }
    println!(
        "xtask fetch-data: {} truth {} wrote {} bytes",
        dataset.id, truth.role, got
    );
    Ok(())
}

fn fetch_reference(root: &Path, dataset: &DataSet, reference: &DataReference) -> DynResult<()> {
    let out_path = root.join(&reference.path);
    let expected = reference.sha256.as_deref();
    if out_path.exists() {
        let got = sha256_file(&out_path)?;
        if expected.map(|sha| sha == got).unwrap_or(false) {
            println!(
                "xtask fetch-data: {} reference {} already present ({got})",
                dataset.id, reference.role
            );
            return Ok(());
        }
    }
    stream_reference_fasta(reference, &out_path)?;
    let stats = fasta_stats(&out_path)?;
    if stats.records == 0 || stats.total_bases == 0 {
        return Err(format!(
            "xtask fetch-data: {} reference {} produced empty FASTA",
            dataset.id, reference.role
        )
        .into());
    }
    let got = sha256_file(&out_path)?;
    match expected {
        Some(expected) if expected == got => {
            println!(
                "xtask fetch-data: {} reference {} OK ({got})",
                dataset.id, reference.role
            );
        }
        Some(expected) => {
            return Err(format!(
                "xtask fetch-data: {} reference {} digest mismatch: got {}, expected {}",
                dataset.id, reference.role, got, expected
            )
            .into());
        }
        None => {
            println!(
                "xtask fetch-data: {} reference {} wrote {} bp sha256={got}",
                dataset.id, reference.role, stats.total_bases
            );
        }
    }
    Ok(())
}

fn stream_reference_fasta(reference: &DataReference, out_path: &Path) -> DynResult<()> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = out_path.with_extension("tmp");
    let script = if reference.url.ends_with(".gz") {
        format!(
            "curl -fsSL '{}' | gzip -dc > '{}'",
            reference.url,
            tmp.display()
        )
    } else {
        format!("curl -fsSL '{}' > '{}'", reference.url, tmp.display())
    };
    let status = Command::new("bash").arg("-c").arg(&script).status()?;
    if !status.success() {
        let _ = fs::remove_file(&tmp);
        return Err(format!("fetch failed for {}", reference.url).into());
    }
    fs::rename(tmp, out_path)?;
    Ok(())
}

fn stream_opaque_file(url: &str, out_path: &Path) -> DynResult<()> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = out_path.with_extension("tmp");
    let script = format!("curl -fsSL '{}' > '{}'", url, tmp.display());
    let status = Command::new("bash").arg("-c").arg(&script).status()?;
    if !status.success() {
        let _ = fs::remove_file(&tmp);
        return Err(format!("fetch failed for {url}").into());
    }
    fs::rename(tmp, out_path)?;
    Ok(())
}

fn stream_fastq_prefix(source: &DataFile, out_path: &Path, records: usize) -> DynResult<()> {
    let lines = records
        .checked_mul(4)
        .ok_or("FASTQ record count overflow")?;
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = out_path.with_extension("tmp");
    let script = format!(
        "curl -fsSL '{}' 2> >(grep -v 'Failure writing output to destination' >&2) | gzip -dc | awk 'NR <= {} {{ print }} NR == {} {{ exit }}' > '{}'",
        source.url,
        lines,
        lines,
        tmp.display()
    );
    let status = Command::new("bash").arg("-c").arg(&script).status()?;
    if !status.success() {
        let _ = fs::remove_file(&tmp);
        return Err(format!("fetch failed for {}", source.url).into());
    }
    fs::rename(tmp, out_path)?;
    Ok(())
}

fn validate_fastq_record_count(path: &Path, expected_records: usize) -> DynResult<()> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let mut records = 0usize;
    while let Some(header) = lines.next() {
        let seq = lines.next().ok_or("truncated FASTQ sequence line")?;
        let plus = lines.next().ok_or("truncated FASTQ plus line")?;
        let qual = lines.next().ok_or("truncated FASTQ quality line")?;
        if !header.starts_with('@') || !plus.starts_with('+') || seq.len() != qual.len() {
            return Err(format!("invalid FASTQ record in {}", path.display()).into());
        }
        records += 1;
    }
    if records != expected_records {
        return Err(format!(
            "{}: expected {} FASTQ records, found {}",
            path.display(),
            expected_records,
            records
        )
        .into());
    }
    Ok(())
}

fn validate_paired_fastq_record_count(
    r1_path: &Path,
    r2_path: &Path,
    expected_records: usize,
) -> DynResult<()> {
    let r1 = fastq_record_names(r1_path)?;
    let r2 = fastq_record_names(r2_path)?;
    if r1.len() != expected_records || r2.len() != expected_records {
        return Err(format!(
            "expected {} paired FASTQ records, found {} in {} and {} in {}",
            expected_records,
            r1.len(),
            r1_path.display(),
            r2.len(),
            r2_path.display()
        )
        .into());
    }
    for (idx, (left, right)) in r1.iter().zip(r2.iter()).enumerate() {
        if left != right {
            return Err(format!(
                "paired FASTQ mate name mismatch at record {}: {} vs {}",
                idx + 1,
                left,
                right
            )
            .into());
        }
    }
    Ok(())
}

fn fastq_record_names(path: &Path) -> DynResult<Vec<String>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let mut names = Vec::new();
    while let Some(header) = lines.next() {
        let seq = lines.next().ok_or("truncated FASTQ sequence line")?;
        let plus = lines.next().ok_or("truncated FASTQ plus line")?;
        let qual = lines.next().ok_or("truncated FASTQ quality line")?;
        if !header.starts_with('@') || !plus.starts_with('+') || seq.len() != qual.len() {
            return Err(format!("invalid FASTQ record in {}", path.display()).into());
        }
        names.push(fastq_core_name(header).to_string());
    }
    Ok(names)
}

fn fastq_core_name(header: &str) -> &str {
    let name = header
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or("");
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}

fn write_fastq_prefix(input: &Path, output: &Path, records: usize) -> DynResult<()> {
    let text = fs::read_to_string(input)?;
    let mut lines = text.lines();
    let mut out = String::new();
    for _ in 0..records {
        let header = lines.next().ok_or_else(|| {
            format!(
                "{}: fewer than {} FASTQ records available",
                input.display(),
                records
            )
        })?;
        let seq = lines.next().ok_or("truncated FASTQ sequence line")?;
        let plus = lines.next().ok_or("truncated FASTQ plus line")?;
        let qual = lines.next().ok_or("truncated FASTQ quality line")?;
        if !header.starts_with('@') || !plus.starts_with('+') || seq.len() != qual.len() {
            return Err(format!("invalid FASTQ record in {}", input.display()).into());
        }
        out.push_str(header);
        out.push('\n');
        out.push_str(seq);
        out.push('\n');
        out.push_str(plus);
        out.push('\n');
        out.push_str(qual);
        out.push('\n');
    }
    fs::write(output, out)?;
    Ok(())
}

fn generate_reads(
    root: &Path,
    reference: &Path,
    out: &Path,
    read_len: usize,
    step: usize,
    circular: bool,
) -> DynResult<()> {
    if read_len == 0 || step == 0 {
        return Err("read_len and step must be greater than zero".into());
    }
    let reference_path = if reference.is_absolute() {
        reference.to_path_buf()
    } else {
        root.join(reference)
    };
    let out_path = if out.is_absolute() {
        out.to_path_buf()
    } else {
        root.join(out)
    };
    let seqs = fasta_sequences(&reference_path)?;
    if seqs.len() != 1 {
        return Err("generate-reads expects exactly one FASTA sequence".into());
    }
    let seq = &seqs[0];
    if seq.len() < read_len {
        return Err("reference shorter than requested read length".into());
    }
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut body = String::new();
    let last_start = if circular {
        seq.len()
    } else {
        seq.len() - read_len + 1
    };
    let mut read_idx = 0usize;
    let mut start = 0usize;
    while start < last_start {
        let mut read = Vec::with_capacity(read_len);
        for offset in 0..read_len {
            let pos = start + offset;
            let base = if circular {
                seq[pos % seq.len()]
            } else {
                seq[pos]
            };
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                return Err(format!("reference contains non-ACGT base {}", base as char).into());
            }
            read.push(base);
        }
        read_idx += 1;
        body.push_str(&format!("@read_{read_idx:06}_start_{start}\n"));
        body.push_str(std::str::from_utf8(&read)?);
        body.push_str("\n+\n");
        body.push_str(&"I".repeat(read_len));
        body.push('\n');
        start += step;
    }
    fs::write(&out_path, body)?;
    println!(
        "xtask generate-reads: wrote {} reads to {}",
        read_idx,
        out_path.display()
    );
    Ok(())
}

fn value_after<'a>(args: &'a [String], key: &str) -> Option<&'a str> {
    args.windows(2)
        .find_map(|pair| (pair[0] == key).then_some(pair[1].as_str()))
}

fn tier_name(tier: Tier) -> &'static str {
    match tier {
        Tier::Pr => "pr",
        Tier::Main => "main",
        Tier::Nightly => "nightly",
        Tier::Manual => "manual",
    }
}

fn run_command_passthrough(root: &Path, program: &str, args: &[&str]) -> DynResult<()> {
    let status = Command::new(program)
        .current_dir(root)
        .args(args)
        .status()?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{program} {} failed with {status}", args.join(" ")).into())
    }
}

fn extract_assemble_flags(path: &Path) -> DynResult<BTreeSet<String>> {
    let text = fs::read_to_string(path)?;
    let mut flags = BTreeSet::new();
    let mut pending_arg: Option<String> = None;
    let mut in_assemble = false;
    let mut brace_depth = 0i32;

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("Assemble {") {
            in_assemble = true;
            brace_depth = 1;
            continue;
        }
        if !in_assemble {
            continue;
        }

        brace_depth += trimmed.matches('{').count() as i32;
        brace_depth -= trimmed.matches('}').count() as i32;
        if brace_depth <= 0 {
            break;
        }

        if trimmed.starts_with("#[arg(") {
            if let Some(long) = extract_long_name(trimmed) {
                flags.insert(format!("--{long}"));
            }
            pending_arg = Some(trimmed.to_string());
            continue;
        }

        if let Some(arg) = pending_arg.as_deref() {
            if !trimmed.starts_with("///") && trimmed.contains(':') {
                if arg.contains("long") && !arg.contains("long =") {
                    let field = trimmed
                        .split(':')
                        .next()
                        .ok_or("invalid field line")?
                        .trim();
                    flags.insert(format!("--{}", field.replace('_', "-")));
                }
                pending_arg = None;
            }
        }
    }

    if flags.is_empty() {
        return Err("no Assemble CLI flags found".into());
    }
    Ok(flags)
}

fn extract_long_name(arg: &str) -> Option<String> {
    let marker = "long = \"";
    let start = arg.find(marker)? + marker.len();
    let rest = &arg[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

fn run_script(root: &Path, script: &str) -> DynResult<ScriptReport> {
    println!("xtask bench: running {script}");
    let script_path = root.join(script);
    let start = std::time::Instant::now();
    let time_path = Path::new("/usr/bin/time");
    let output = if time_path.exists() {
        Command::new(time_path)
            .current_dir(root)
            .args(["-f", "TREX_XTASK_TIME\t%e\t%M", "bash"])
            .arg(&script_path)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    } else {
        Command::new("bash")
            .current_dir(root)
            .arg(&script_path)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?
    };
    let wall_ms = start.elapsed().as_millis();

    let stdout = String::from_utf8_lossy(&output.stdout);
    if !stdout.is_empty() {
        print!("{stdout}");
    }

    let stderr = String::from_utf8_lossy(&output.stderr);
    let (time_seconds, max_rss_kib, passthrough_stderr) = parse_time_stderr(&stderr);
    if !passthrough_stderr.trim().is_empty() {
        eprint!("{passthrough_stderr}");
        if !passthrough_stderr.ends_with('\n') {
            eprintln!();
        }
    }

    Ok(ScriptReport {
        path: script.to_string(),
        exit_code: output.status.code(),
        success: output.status.success(),
        wall_ms,
        time_seconds,
        max_rss_kib,
    })
}

fn parse_time_stderr(stderr: &str) -> (Option<f64>, Option<u64>, String) {
    let mut time_seconds = None;
    let mut max_rss_kib = None;
    let mut passthrough = String::new();

    for line in stderr.lines() {
        if let Some(rest) = line.strip_prefix("TREX_XTASK_TIME\t") {
            let mut parts = rest.split('\t');
            let secs = parts.next().and_then(|s| s.parse::<f64>().ok());
            let rss = parts.next().and_then(|s| s.parse::<u64>().ok());
            time_seconds = secs;
            max_rss_kib = rss;
        } else {
            passthrough.push_str(line);
            passthrough.push('\n');
        }
    }

    (time_seconds, max_rss_kib, passthrough)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_root(name: &str) -> PathBuf {
        let millis = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_millis();
        let root =
            std::env::temp_dir().join(format!("trex-xtask-{name}-{}-{millis}", std::process::id()));
        fs::create_dir_all(&root).expect("create temp root");
        root
    }

    #[test]
    fn reference_quality_counts_forward_and_reverse_kmers() {
        let reference = vec![b"ACCGTTAA".to_vec()];
        let contigs = vec![b"AACGGT".to_vec()];
        let q = reference_quality_from_sequences(&reference, &contigs, 4);

        assert_eq!(q.contig_records, 1);
        assert_eq!(q.contig_kmers, 3);
        assert_eq!(q.contig_kmers_in_reference, 3);
        assert_eq!(q.contigs_with_reference_hit, 1);
        assert_eq!(q.contig_bases_with_reference_hit, 6);
        assert_eq!(q.contig_kmer_reference_fraction, 1.0);
    }

    #[test]
    fn reference_quality_ignores_non_acgt_windows() {
        let reference = vec![b"ACGTACGT".to_vec()];
        let contigs = vec![b"ACGTNNNN".to_vec()];
        let q = reference_quality_from_sequences(&reference, &contigs, 4);

        assert_eq!(q.contig_kmers, 1);
        assert_eq!(q.contig_kmers_in_reference, 1);
        assert_eq!(q.contigs_with_reference_hit, 1);
    }

    #[test]
    fn matrix_deserializes_parental_references() {
        let text = r#"
schema_version = 1

[[rows]]
id = "phase2"
technology = "illumina_pe"
organism = "synthetic"
license = "in_repo"
provenance = "fixtures"
depth_class = "smoke"
ci_tier = "manual"
product_claim = "observed parent-specific metric"
claim_level = "observed"
reference_availability = "parental_haplotypes"
artifact_policy = "required"
fixtures = ["fixtures/phase2_synthetic/reads.fq"]
required_tools = ["cargo"]
artifacts = ["target/x/contigs.fa"]

[rows.trex]
tiers = ["manual"]
out_dir = "target/x"
args = ["illumina", "assemble", "--r1", "fixtures/phase2_synthetic/reads.fq", "--kmer-size", "3", "--out-dir", "target/x"]
parental_references = [
  { label = "parent1", path = "fixtures/phase2_synthetic/parent1.fa" },
  { label = "parent2", path = "fixtures/phase2_synthetic/parent2.fa" },
]
"#;
        let matrix: Matrix = toml::from_str(text).expect("matrix");
        let trex = matrix.rows[0].trex.as_ref().expect("trex");
        assert_eq!(trex.parental_references.len(), 2);
        assert_eq!(trex.parental_references[0].label, "parent1");
        assert_eq!(
            trex.parental_references[1].path,
            "fixtures/phase2_synthetic/parent2.fa"
        );
    }

    #[test]
    fn parental_reference_quality_reports_each_parent() {
        let root = temp_root("parental-quality");
        fs::create_dir_all(root.join("parents")).expect("parents dir");
        fs::create_dir_all(root.join("out")).expect("out dir");
        fs::write(root.join("parents/p1.fa"), b">p1\nACGTACGT\n").expect("p1");
        fs::write(root.join("parents/p2.fa"), b">p2\nTTTTTTTT\n").expect("p2");
        fs::write(root.join("out/contigs.fa"), b">ctg\nACGTAC\n").expect("contigs");

        let trex = TrexBench {
            tiers: vec![Tier::Manual],
            args: vec![],
            out_dir: "out".to_string(),
            reference: None,
            parental_references: vec![
                ParentReference {
                    label: "p1".to_string(),
                    path: "parents/p1.fa".to_string(),
                },
                ParentReference {
                    label: "p2".to_string(),
                    path: "parents/p2.fa".to_string(),
                },
            ],
        };
        let q = parental_reference_quality(&root, &trex, &root.join("out/contigs.fa"), 3)
            .expect("parental quality");
        assert_eq!(q.len(), 2);
        assert_eq!(q[0].label, "p1");
        assert_eq!(q[0].quality.contig_kmer_reference_fraction, 1.0);
        assert_eq!(q[1].label, "p2");
        assert_eq!(q[1].quality.contig_kmer_reference_fraction, 0.0);

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn read_assembly_quality_reports_containment_and_assembly_only_kmers() {
        let reads = vec![b"ACGTACGT".to_vec()];
        let assembly = vec![b"ACGTTCGT".to_vec()];
        let q = read_assembly_kmer_quality_from_sequences(&reads, &assembly, 3, 2);

        assert_eq!(q.read_distinct_kmers, 2);
        assert_eq!(q.reliable_read_kmers, 2);
        assert_eq!(q.assembly_distinct_kmers, 4);
        assert_eq!(q.assembly_total_kmers, 6);
        assert_eq!(q.reliable_read_kmers_in_assembly, 1);
        assert_eq!(q.reliable_read_containment_fraction, 0.5);
        assert_eq!(q.assembly_only_kmers, 3);
        assert_eq!(q.assembly_only_total_kmers, 3);
        assert_eq!(q.assembly_only_fraction, 0.75);
        assert_eq!(q.assembly_only_total_fraction, 0.5);
        let qv = q
            .approximate_qv
            .expect("nonzero assembly-only fraction has QV");
        assert!((qv - 1.249_387).abs() < 0.000_001);
    }

    #[test]
    fn read_assembly_quality_uses_none_for_zero_error_qv() {
        let reads = vec![b"ACGTACGT".to_vec()];
        let assembly = vec![b"ACGTACGT".to_vec()];
        let q = read_assembly_kmer_quality_from_sequences(&reads, &assembly, 3, 1);

        assert_eq!(q.assembly_only_kmers, 0);
        assert_eq!(q.assembly_only_total_kmers, 0);
        assert_eq!(q.assembly_only_total_fraction, 0.0);
        assert_eq!(q.qv_error_rate, 0.0);
        assert_eq!(q.approximate_qv, None);
    }

    #[test]
    fn layer_status_reports_missing_required_tool_first() {
        let tools = ToolAvailability {
            required_available: Vec::new(),
            required_unavailable: vec!["definitely_missing_trex_tool".to_string()],
            optional_available: Vec::new(),
            optional_unavailable: Vec::new(),
        };

        let status = row_layer_status(&tools, &[], &[]);

        assert_eq!(status.status, "failed");
        assert_eq!(
            status.failed_layer.as_deref(),
            Some("required_tool:definitely_missing_trex_tool")
        );
    }

    #[test]
    fn fastq_core_name_normalizes_pair_suffixes() {
        assert_eq!(fastq_core_name("@read/1 extra"), "read");
        assert_eq!(fastq_core_name("@read/2"), "read");
        assert_eq!(
            fastq_core_name("@HISEQ1:26:HA2RRADXX:1:2203"),
            "HISEQ1:26:HA2RRADXX:1:2203"
        );
    }

    #[test]
    fn parse_region_accepts_1_based_closed_coordinates() {
        assert_eq!(
            parse_region("chr20:10,000,000-10,010,000").expect("region"),
            GenomicRegion {
                contig: "chr20".to_string(),
                start: 10_000_000,
                end: 10_010_000,
            }
        );
        assert!(parse_region("chr20:0-10").is_err());
        assert!(parse_region("chr20:20-10").is_err());
    }
}
