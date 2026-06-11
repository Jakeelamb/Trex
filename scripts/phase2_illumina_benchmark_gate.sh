#!/usr/bin/env bash
# **Phase-2 Illumina benchmark gate** (CONTEXT): layer 1 = Phase-1 gate; layer 2+ = synthetic diploid + graph + haplotype metrics.
# Set `TREX_RUN_QUAST=1` to opt into **QUAST** when installed. Local devs keep full Phase-1 layer via `benchmark_gate.sh`;
# PR CI uses a slimmer path in `.github/workflows/ci.yml` (see CONTEXT: Phase-2 Illumina scope and CI cadence).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

run_layer() {
  local code="$1"
  local name="$2"
  shift 2
  echo "phase2_illumina_benchmark_gate: layer ${name}"
  if ! "$@"; then
    echo "phase2_illumina_benchmark_gate: layer ${name} failed" >&2
    exit "$code"
  fi
}

run_layer 10 "phase1_benchmark_gate" bash "${ROOT}/scripts/benchmark_gate.sh"
run_layer 20 "phase2_diploid_reference" bash "${ROOT}/scripts/phase2_illumina_diploid_reference_layer.sh"
run_layer 30 "phase2_graph_summaries" bash "${ROOT}/scripts/phase2_illumina_graph_summaries.sh"
run_layer 40 "phase2_haplotype_metrics" bash "${ROOT}/scripts/phase2_illumina_haplotype_metrics.sh"
# Optional **QUAST**-class row (see `scripts/reference_quast.sh`); skipped when quast not installed.
if [[ "${TREX_RUN_QUAST:-0}" == "1" ]]; then
  run_layer 50 "optional_quast" bash "${ROOT}/scripts/reference_quast.sh"
fi
echo "phase2_illumina_benchmark_gate: all layers OK"
