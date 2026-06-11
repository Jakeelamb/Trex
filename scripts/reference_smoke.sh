#!/usr/bin/env bash
# Reference-aligned smoke (**Phase-1 benchmark gate**): verify contigs against `fixtures/tiny_ref.fa`.
# Uses **minimap2** when it emits PAF; otherwise falls back to substring check (ultra-short contigs).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${REF_FREE_SMOKE_OUT:-${ROOT}/target/ref-free-smoke}"
REF="${ROOT}/fixtures/tiny_ref.fa"
CONTIGS="${OUT}/contigs.fa"

if [[ ! -f "${CONTIGS}" ]]; then
  echo "reference_smoke: missing ${CONTIGS} (run ref_free_smoke first)" >&2
  exit 1
fi

verify_substring() {
  python3 - "$REF" "$CONTIGS" <<'PY'
import sys
from pathlib import Path

def load_seq(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            continue
        parts.append(line.strip())
    return "".join(parts)

ref = load_seq(Path(sys.argv[1]))
ctg = load_seq(Path(sys.argv[2]))
if ctg not in ref:
    sys.stderr.write(
        f"reference_smoke: contig not substring of ref (len ctg={len(ctg)} ref={len(ref)})\n"
    )
    sys.exit(1)
PY
}

if command -v minimap2 >/dev/null 2>&1; then
  PAF="${OUT}/tiny_ref.paf"
  if minimap2 -x sr "${REF}" "${CONTIGS}" > "${PAF}" 2>/dev/null; then
    if [[ -s "${PAF}" ]]; then
      echo "reference_smoke: minimap2 PAF ${PAF} ($(wc -l < "${PAF}") lines)"
      exit 0
    fi
  fi
  echo "reference_smoke: minimap2 produced no PAF (short query); using substring check"
fi

verify_substring
echo "reference_smoke: substring check OK"
