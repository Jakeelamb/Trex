use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::error::Error;
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
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
    fixtures: Option<Vec<String>>,
    external_data: Option<String>,
    digest_manifest: Option<String>,
    manifest_table: Option<String>,
    pr_scripts: Option<Vec<String>>,
    main_scripts: Option<Vec<String>>,
    nightly_scripts: Option<Vec<String>>,
    manual_scripts: Option<Vec<String>>,
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
    files: Vec<DataFile>,
    #[serde(default)]
    references: Vec<DataReference>,
    #[serde(default)]
    prepared: Vec<PreparedDataFile>,
}

#[derive(Debug, Deserialize)]
struct DataFile {
    role: String,
    url: String,
    md5: String,
    bytes: u64,
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
struct DataReference {
    role: String,
    url: String,
    path: String,
    sha256: Option<String>,
    notes: Option<String>,
}

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
    scripts: Vec<ScriptReport>,
    trex_runs: Vec<TrexRunReport>,
    artifacts: Vec<ArtifactReport>,
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
    reference_quality: Option<ReferenceQuality>,
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
        }
        Cmd::ValidateMatrix => validate_matrix(&root)?,
        Cmd::ValidateCapabilities => validate_capabilities(&root)?,
        Cmd::ValidateData => validate_data_catalog(&root)?,
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

        for artifact in row.artifacts.as_deref().unwrap_or(&[]) {
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

        for optional_tool in row.optional_tools.as_deref().unwrap_or(&[]) {
            if optional_tool.trim().is_empty() {
                return Err(format!("{row_id}: optional_tools contains an empty entry").into());
            }
        }
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
            if file.md5.len() != 32 || !file.md5.chars().all(|c| c.is_ascii_hexdigit()) {
                return Err(format!("{}: invalid md5 for {}", dataset.id, file.role).into());
            }
            if file.bytes == 0 {
                return Err(
                    format!("{}: bytes must be non-zero for {}", dataset.id, file.role).into(),
                );
            }
        }
        if roles.is_empty() {
            return Err(format!("{}: at least one source file is required", dataset.id).into());
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
        for script in scripts {
            let report = run_script(root, script)?;
            if !report.success {
                failed = true;
            }
            script_reports.push(report);
        }

        let mut trex_reports = Vec::new();
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
        return Err("xtask bench: one or more scripts failed".into());
    }
    Ok(())
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
    Ok(AssemblyMetrics {
        kmer_size,
        candidate_kmers,
        reads,
        r2_reads,
        reference,
        contigs,
        unitigs,
        gfa,
        reference_quality,
    })
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
    let text = fs::read_to_string(path)?;
    let mut records = 0;
    let mut total_bases = 0;
    let mut min_len = usize::MAX;
    let mut max_len = 0;
    let mut candidate_kmers = kmer_size.map(|_| 0usize);
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
        records += 1;
        total_bases += seq.len();
        min_len = min_len.min(seq.len());
        max_len = max_len.max(seq.len());
        if let (Some(k), Some(total)) = (kmer_size, candidate_kmers.as_mut()) {
            if seq.len() >= k {
                *total += seq.len() - k + 1;
            }
        }
    }
    if records == 0 {
        min_len = 0;
    }
    Ok(FastqStats {
        records,
        total_bases,
        min_len,
        max_len,
        candidate_kmers,
    })
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

fn reference_quality(
    reference_path: &Path,
    contigs_path: &Path,
    k: usize,
) -> DynResult<ReferenceQuality> {
    let reference = fasta_sequences(reference_path)?;
    let contigs = fasta_sequences(contigs_path)?;
    Ok(reference_quality_from_sequences(&reference, &contigs, k))
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
    for reference in &dataset.references {
        fetch_reference(root, dataset, reference)?;
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
}
