#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/target/ref-free-smoke"
mkdir -p "${OUT}"
export REF_FREE_SMOKE_OUT="${OUT}"

cargo run -q -p trex-cli -- illumina assemble \
  --r1 "${ROOT}/fixtures/tiny.fq" \
  --kmer-size 4 \
  --trusted-threshold 1 \
  --out-dir "${OUT}"

python3 - <<'PY'
import os
from pathlib import Path

def fasta_stats(path: Path):
    seqs = []
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
        print(f"{path.name}: (empty)")
        return
    total = sum(lens)
    lens_sorted = sorted(lens, reverse=True)
    acc = 0
    n50 = 0
    for L in lens_sorted:
        acc += L
        if acc * 2 >= total:
            n50 = L
            break
    print(f"{path.name}: contigs={len(lens)} total_bases={total} n50={n50}")

out = Path(os.environ["REF_FREE_SMOKE_OUT"])
fasta_stats(out / "contigs.fa")
fasta_stats(out / "unitigs.fa")
PY

EXP="${ROOT}/fixtures/expected/ref_free_smoke"
for f in contigs.fa unitigs.fa graph.gfa; do
  if [[ ! -f "${EXP}/${f}" ]]; then
    echo "ref_free_smoke: missing golden file ${EXP}/${f}" >&2
    exit 1
  fi
  if ! cmp -s "${OUT}/${f}" "${EXP}/${f}"; then
    echo "ref_free_smoke: ${f} differs from fixtures/expected (run from repo root; refresh with fixtures/README.md)" >&2
    diff -u "${EXP}/${f}" "${OUT}/${f}" >&2 || true
    exit 1
  fi
done
echo "ref_free_smoke: output matches fixtures/expected/ref_free_smoke"
