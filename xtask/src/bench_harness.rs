use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::{
    artifacts_for_tier, load_matrix, run_script, run_trex_bench, scripts_for_tier, validate_matrix,
    ArtifactReport, BenchReport, BenchRowReport, DynResult, RowLayerStatus, ScriptReport, Tier,
    ToolAvailability, TrexRunReport,
};

pub(crate) fn run_bench(
    root: &Path,
    tier: Tier,
    row_filter: Option<&str>,
    out: &Path,
) -> DynResult<()> {
    validate_matrix(root)?;
    let matrix = load_matrix(root)?;
    let harness = BenchmarkHarness::new(root);
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
        let comparison_threads = row.comparison_threads;
        let tools = harness.tool_availability(
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
        let layer_status = harness.row_layer_status(&tools, &script_reports, &trex_reports);
        let artifacts = harness.artifact_reports(artifacts_for_tier(row, tier));

        rows.push(BenchRowReport {
            id: row_id,
            ci_tier,
            product_claim,
            claim_level,
            reference_availability,
            artifact_policy,
            comparison_threads,
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

    let out_path = harness.write_report(tier, rows, out)?;
    println!("xtask bench: wrote {}", out_path.display());

    if failed {
        return Err("xtask bench: one or more rows failed".into());
    }
    Ok(())
}

pub(crate) struct BenchmarkHarness<'a> {
    root: &'a Path,
}

impl<'a> BenchmarkHarness<'a> {
    pub(crate) fn new(root: &'a Path) -> Self {
        Self { root }
    }

    pub(crate) fn tool_availability(
        &self,
        required_tools: &[String],
        optional_tools: &[String],
    ) -> ToolAvailability {
        let (required_available, required_unavailable) = partition_tools(required_tools);
        let (optional_available, optional_unavailable) = partition_tools(optional_tools);
        ToolAvailability {
            required_available,
            required_unavailable,
            optional_available,
            optional_unavailable,
        }
    }

    pub(crate) fn row_layer_status(
        &self,
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

    pub(crate) fn artifact_reports<'p>(
        &self,
        paths: impl IntoIterator<Item = &'p String>,
    ) -> Vec<ArtifactReport> {
        paths
            .into_iter()
            .map(|path| {
                let meta = fs::metadata(self.root.join(path)).ok();
                ArtifactReport {
                    path: path.to_string(),
                    exists: meta.is_some(),
                    bytes: meta.map(|m| m.len()),
                }
            })
            .collect()
    }

    pub(crate) fn write_report(
        &self,
        tier: Tier,
        rows: Vec<BenchRowReport>,
        out: &Path,
    ) -> DynResult<PathBuf> {
        let report = BenchReport {
            schema_version: 1,
            tier,
            generated_unix_ms: SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
            rows,
        };
        let out_path = if out.is_absolute() {
            out.to_path_buf()
        } else {
            self.root.join(out)
        };
        if let Some(parent) = out_path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(&out_path, serde_json::to_string_pretty(&report)? + "\n")?;
        Ok(out_path)
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
