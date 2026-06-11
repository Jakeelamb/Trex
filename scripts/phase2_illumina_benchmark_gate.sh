#!/usr/bin/env bash
# **Phase-2 Illumina benchmark gate** (CONTEXT): layer 1 = Phase-1 gate; layer 2+ = synthetic diploid + graph + haplotype metrics.
# Set `TREX_RUN_QUAST=1` to opt into **QUAST** when installed. Local devs keep full Phase-1 layer via `benchmark_gate.sh`;
# PR CI uses a slimmer path in `.github/workflows/ci.yml` (see CONTEXT: Phase-2 Illumina scope and CI cadence).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "${ROOT}/scripts/benchmark_gate.sh"
bash "${ROOT}/scripts/phase2_illumina_diploid_reference_layer.sh"
bash "${ROOT}/scripts/phase2_illumina_graph_summaries.sh"
bash "${ROOT}/scripts/phase2_illumina_haplotype_metrics.sh"
# Optional **QUAST**-class row (see `scripts/reference_quast.sh`); skipped when quast not installed.
if [[ "${TREX_RUN_QUAST:-0}" == "1" ]]; then
  bash "${ROOT}/scripts/reference_quast.sh"
fi
echo "phase2_illumina_benchmark_gate: all layers OK"
