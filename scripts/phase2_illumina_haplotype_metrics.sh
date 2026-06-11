#!/usr/bin/env bash
# **Phase-2 Illumina** haplotype-aware reference leg: best parent + Hamming (synthetic fixtures).
# Requires `contigs.fa` from a diploid `trex illumina assemble` run (e.g. target/phase2-graph-summaries).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-${ROOT}/target/phase2-graph-summaries}"
P1="${ROOT}/fixtures/phase2_synthetic/parent1.fa"
P2="${ROOT}/fixtures/phase2_synthetic/parent2.fa"
CF="${OUT}/contigs.fa"

for f in "$P1" "$P2" "$CF"; do
  if [[ ! -f "$f" ]]; then
    echo "phase2_illumina_haplotype_metrics: missing $f (run graph summaries / assemble first)" >&2
    exit 1
  fi
done

python3 - "$P1" "$P2" "$CF" <<'PY'
import sys
from pathlib import Path

def load_fasta(path: Path) -> list[bytes]:
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
    return seqs

def ham(a: bytes, b: bytes) -> int:
    return sum(x != y for x, y in zip(a, b))

def best_substring_ham(query: bytes, ref: bytes) -> int:
    L, R = len(query), len(ref)
    if L > R or L == 0:
        return R + 1
    best = R + 1
    for j in range(R - L + 1):
        d = ham(query, ref[j : j + L])
        if d < best:
            best = d
    return best

p1 = load_fasta(Path(sys.argv[1]))[0]
p2 = load_fasta(Path(sys.argv[2]))[0]
contigs = load_fasta(Path(sys.argv[3]))

if len(p1) != len(p2):
    sys.exit("parents must be same length")

for i, c in enumerate(contigs):
    if len(c) == len(p1):
        d1, d2 = ham(c, p1), ham(c, p2)
        parent = "p1" if d1 <= d2 else "p2"
        dist = min(d1, d2)
    else:
        d1 = best_substring_ham(c, p1)
        d2 = best_substring_ham(c, p2)
        parent = "p1" if d1 <= d2 else "p2"
        dist = min(d1, d2)
    print(f"  ctg{i+1}: len={len(c)} best_parent={parent} hamming={dist}")

print("phase2_illumina_haplotype_metrics: OK")
PY
