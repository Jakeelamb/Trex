use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::{
    ArtifactReport, BenchReport, BenchRowReport, DynResult, RowLayerStatus, ScriptReport, Tier,
    ToolAvailability, TrexRunReport,
};

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
