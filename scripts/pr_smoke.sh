#!/usr/bin/env bash
# **Phase-1 PR smoke** + lightweight **Phase-2 Illumina** layers (no minimap2 / no full matrix).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
bash scripts/ref_free_smoke.sh
bash scripts/phase2_illumina_diploid_reference_layer.sh
bash scripts/phase2_illumina_graph_summaries.sh
bash scripts/phase2_illumina_haplotype_metrics.sh
cargo test --workspace --all-features -q
echo "pr_smoke: OK"
