#!/usr/bin/env bash
# **Phase-1 benchmark gate**: reference-free stats + optional reference alignment (see `CONTEXT.md`).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "${ROOT}/scripts/ref_free_smoke.sh"
bash "${ROOT}/scripts/reference_smoke.sh"
