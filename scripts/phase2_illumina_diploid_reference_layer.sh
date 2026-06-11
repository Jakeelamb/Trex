#!/usr/bin/env bash
# Phase-2 Illumina benchmark gate — diploid **reference-aligned** layer (synthetic fixture).
# Validates two parental haplotypes + reads per **Phase-2 Illumina benchmark gate** until Trex diploid assembly ships.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
P1="${ROOT}/fixtures/phase2_synthetic/parent1.fa"
P2="${ROOT}/fixtures/phase2_synthetic/parent2.fa"
FQ="${ROOT}/fixtures/phase2_synthetic/reads.fq"
MANIFEST="${ROOT}/tools/manifest.toml"

for f in "$P1" "$P2" "$FQ" "$MANIFEST"; do
  if [[ ! -f "$f" ]]; then
    echo "phase2_illumina_diploid_reference_layer: missing $f" >&2
    exit 1
  fi
done

digest_field() {
  local key="$1"
  awk -F '"' -v "k=${key}" '
    $0 ~ "^" k " *= *" { print $2; exit }
  ' "${MANIFEST}"
}

verify_sha256() {
  local file="$1" want="$2"
  local got
  got=$(sha256sum "$file" | awk '{print $1}')
  if [[ "$got" != "$want" ]]; then
    echo "phase2: digest mismatch ${file}: got ${got} expected ${want}" >&2
    exit 1
  fi
}

verify_sha256 "$P1" "$(digest_field parent1_fa_sha256)"
verify_sha256 "$P2" "$(digest_field parent2_fa_sha256)"
verify_sha256 "$FQ" "$(digest_field reads_fq_sha256)"
echo "phase2_illumina_diploid_reference_layer: fixture digests OK (tools/manifest.toml)"

python3 - "$P1" "$P2" "$FQ" <<'PY'
import sys
from pathlib import Path

def load_fasta_seq(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            continue
        parts.append(line.strip())
    return "".join(parts)

def load_fastq_seqs(path: Path) -> list[str]:
    seqs: list[str] = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("@"):
            if i + 1 >= len(lines):
                raise SystemExit(f"malformed FASTQ at line {i+1}")
            seqs.append(lines[i + 1].strip())
            i += 4
        else:
            i += 1
    return seqs

p1 = load_fasta_seq(Path(sys.argv[1]))
p2 = load_fasta_seq(Path(sys.argv[2]))
reads = load_fastq_seqs(Path(sys.argv[3]))

if len(p1) != len(p2):
    sys.stderr.write("phase2: parental sequences must be same length\n")
    sys.exit(1)
if p1 == p2:
    sys.stderr.write("phase2: parental sequences must differ\n")
    sys.exit(1)

for i, seq in enumerate(reads):
    if not seq:
        sys.stderr.write(f"phase2: empty read sequence at index {i}\n")
        sys.exit(1)
    if not all(c in "ACGTacgt" for c in seq):
        sys.stderr.write(f"phase2: non-ACGT read at index {i}\n")
        sys.exit(1)
    u = seq.upper()
    if u not in p1 and u not in p2:
        sys.stderr.write(
            f"phase2: read {i+1} not substring of either parent\n"
        )
        sys.exit(1)

print(f"phase2_illumina_diploid_reference_layer: OK ({len(reads)} reads, parents len={len(p1)})")
PY

if command -v minimap2 >/dev/null 2>&1; then
  for ref in "$P1" "$P2"; do
    PAF="${ROOT}/target/phase2-synthetic-$(basename "$ref" .fa).paf"
    mkdir -p "${ROOT}/target/phase2-synthetic"
    if minimap2 -x sr "$ref" "$FQ" >"$PAF" 2>/dev/null && [[ -s "$PAF" ]]; then
      echo "phase2_illumina_diploid_reference_layer: minimap2 $(basename "$ref") -> $(wc -l <"$PAF") PAF lines"
    else
      echo "phase2_illumina_diploid_reference_layer: minimap2 empty PAF for $(basename "$ref"); substring check already passed"
    fi
  done
fi
