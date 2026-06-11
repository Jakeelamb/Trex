use std::collections::BTreeSet;
use std::error::Error;
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};

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
    /// Run benchmark matrix scripts for one tier and write a JSON artifact.
    Bench {
        #[arg(long, value_enum, default_value_t = Tier::Pr)]
        tier: Tier,
        #[arg(long, default_value = "target/benchmarks/matrix.json")]
        out: PathBuf,
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
        }
        Cmd::ValidateMatrix => validate_matrix(&root)?,
        Cmd::ValidateCapabilities => validate_capabilities(&root)?,
        Cmd::Bench { tier, out } => run_bench(&root, tier, &out)?,
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
    let text = fs::read_to_string(&path)?;
    Ok(toml::from_str(&text)?)
}

fn validate_matrix(root: &Path) -> DynResult<()> {
    let matrix = load_matrix(root)?;
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

        let fixtures = required_list(row.fixtures.as_ref(), row_id, "fixtures")?;
        for fixture in fixtures {
            require_rel_path(root, row_id, "fixture", fixture, true)?;
        }

        let script_paths = scripts_for_any_tier(row);
        if script_paths.is_empty() {
            return Err(format!("{row_id}: at least one *_scripts list is required").into());
        }
        if ci_tier == Tier::Pr && row.pr_scripts.as_deref().unwrap_or(&[]).is_empty() {
            return Err(format!("{row_id}: ci_tier=pr requires pr_scripts").into());
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
    }

    println!("xtask validate-matrix: OK ({} rows)", matrix.rows.len());
    Ok(())
}

fn validate_capabilities(root: &Path) -> DynResult<()> {
    let doc_path = root.join("docs/CAPABILITIES.md");
    let text = fs::read_to_string(&doc_path)?;
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

fn run_bench(root: &Path, tier: Tier, out: &Path) -> DynResult<()> {
    validate_matrix(root)?;
    let matrix = load_matrix(root)?;
    let mut rows = Vec::new();
    let mut failed = false;

    for row in &matrix.rows {
        let row_id = row.id.clone().ok_or("validated row missing id")?;
        let ci_tier = row.ci_tier.ok_or("validated row missing ci_tier")?;
        let scripts = scripts_for_tier(row, tier);
        if scripts.is_empty() {
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
            artifacts,
        });
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

fn scripts_for_any_tier(row: &Row) -> Vec<&String> {
    let mut scripts = Vec::new();
    for list in [
        row.pr_scripts.as_ref(),
        row.main_scripts.as_ref(),
        row.nightly_scripts.as_ref(),
        row.manual_scripts.as_ref(),
    ] {
        if let Some(list) = list {
            scripts.extend(list.iter());
        }
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
    ] {
        if let Some(list) = list {
            artifacts.extend(list.iter());
        }
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
