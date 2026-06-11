#!/usr/bin/env bash
# **Phase-1 reference metrics toolchain** hook: run MetaQuast / QUAST when available (CONTEXT).
# Exits 0 when skipped; exits non-zero only when QUAST is invoked and fails.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${ROOT}/fixtures/phase2_synthetic/parent1.fa"
ASM="${ROOT}/target/phase2-graph-summaries/contigs.fa"
OUT="${ROOT}/target/quast-phase2-synthetic"

if [[ ! -f "$REF" || ! -f "$ASM" ]]; then
  echo "reference_quast: skip (missing ref or assembly; run phase2 graph summaries first)"
  exit 0
fi

# Default **1** so smoke-scale synthetic contigs (e.g. Phase-2 fixture) are not dropped (QUAST default is 500).
MIN_CONTIG="${TREX_QUAST_MIN_CONTIG:-1}"
MIN_ALIGN="${TREX_QUAST_MIN_ALIGNMENT:-1}"

if command -v metaquast.py >/dev/null 2>&1; then
  mkdir -p "$OUT"
  metaquast.py "$ASM" -r "$REF" -o "$OUT" --no-html \
    --min-contig "$MIN_CONTIG" --min-alignment "$MIN_ALIGN"
  echo "reference_quast: metaquast.py finished -> $OUT"
  exit 0
fi

if command -v quast.py >/dev/null 2>&1; then
  mkdir -p "$OUT"
  quast.py "$ASM" -r "$REF" -o "$OUT" --no-html \
    --min-contig "$MIN_CONTIG" --min-alignment "$MIN_ALIGN"
  echo "reference_quast: quast.py finished -> $OUT"
  exit 0
fi

echo "reference_quast: QUAST not installed; skipping (install quast for CI/local matrix row)"
