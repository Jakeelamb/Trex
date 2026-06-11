#!/usr/bin/env bash
# **Phase-2 Illumina graph summaries** (CONTEXT): GFA topology checks + **Phase-1 reference-free metrics** on **primary** contigs.fa only (not labeled as extending Phase-1 gate vocabulary).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FQ="${ROOT}/fixtures/phase2_synthetic/reads.fq"
OUT="${ROOT}/target/phase2-graph-summaries"
mkdir -p "${OUT}"

cargo run -q -p trex-cli -- illumina assemble \
  --r1 "${FQ}" \
  --kmer-size 4 \
  --trusted-threshold 1 \
  --diploid \
  --out-dir "${OUT}"

python3 - "${OUT}" <<'PY'
import os
import re
import sys
from pathlib import Path

out = Path(sys.argv[1])
contigs = out / "contigs.fa"
gfa = out / "graph.gfa"
if not contigs.is_file() or not gfa.is_file():
    sys.stderr.write("phase2_illumina_graph_summaries: missing contigs.fa or graph.gfa\n")
    sys.exit(1)

# Phase-1 reference-free metrics on primary FASTA only (naming kept for operator clarity).
def fasta_stats(path: Path) -> tuple[int, int, int]:
    seqs: list[bytes] = []
    cur: list[bytes] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append(b"".join(cur))
            cur = []
        else:
            cur.append(line.strip().encode())
    if cur:
        seqs.append(b"".join(cur))
    lens = [len(s) for s in seqs]
    if not lens:
        return 0, 0, 0
    total = sum(lens)
    lens_sorted = sorted(lens, reverse=True)
    acc = 0
    n50 = 0
    for L in lens_sorted:
        acc += L
        if acc * 2 >= total:
            n50 = L
            break
    return len(lens), total, n50

ct, total, n50 = fasta_stats(contigs)
print(
    f"phase2_primary_contigs.fa (Phase-1-style ref-free stats): contigs={ct} total_bases={total} n50={n50}"
)

text = gfa.read_text()
if "trex-phase2-illumina" not in text:
    sys.stderr.write(
        "phase2_illumina_graph_summaries: expected Phase-2 GFA header tag trex-phase2-illumina\n"
    )
    sys.exit(1)
s_lines = len(re.findall(r"(?m)^S\t", text))
l_lines = len(re.findall(r"(?m)^L\t", text))
p_lines = len(re.findall(r"(?m)^P\t", text))
print(
    f"phase2_illumina_graph_summaries: GFA S_lines={s_lines} L_lines={l_lines} P_lines={p_lines}"
)

if s_lines < 1:
    sys.stderr.write("phase2_illumina_graph_summaries: expected at least one GFA S line\n")
    sys.exit(1)
PY

echo "phase2_illumina_graph_summaries: OK"
