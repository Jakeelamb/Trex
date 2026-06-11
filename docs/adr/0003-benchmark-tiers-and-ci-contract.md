# ADR 0003 — Benchmark tiers and CI contract

## Status

Accepted.

## Context

Trex needs architecture and product claims to have local or CI pass/fail loops. The repository already has Phase-1 smoke scripts, Phase-2 Illumina synthetic scripts, and a benchmark matrix file, but the matrix must be enforced and CI tiers must be explicit so benchmark scope does not drift through prose alone.

## Decisions

1. **CI tiers**
   - Pull requests run contract validators plus `scripts/pr_smoke.sh` on MSRV and stable Rust.
   - Pushes to `main` / `master`, tags, schedules, and manual dispatch run the full `scripts/phase2_illumina_benchmark_gate.sh` after installing minimap2.
   - QUAST remains opt-in with `TREX_RUN_QUAST=1`; it is not required for the default PR or main gate until the tool installation and artifact policy are made mandatory.

2. **Matrix-as-code**
   - `tools/benchmark_matrix.toml` is a schema-validated governed matrix.
   - Every row names its minimum CI tier, fixtures, scripts, and optional artifacts or external tools.
   - `tools/validate_benchmark_matrix.py` fails CI when a row omits required fields or points at missing fixtures/scripts.

3. **Capabilities sync**
   - `docs/CAPABILITIES.md` is the single operator capability table for CLI flags, outputs, checkpoints, benchmark scripts, stability, and deferred work.
   - `tools/validate_capabilities_doc.py` fails CI when current `trex illumina assemble` flags or benchmark scripts are missing from that page.

4. **Layer failure semantics**
   - `scripts/phase2_illumina_benchmark_gate.sh` uses distinct exit codes for its layers: 10 for Phase-1 gate, 20 for Phase-2 reference, 30 for graph summaries, 40 for haplotype metrics, and 50 for optional QUAST.
   - CI still runs most PR layers through `scripts/pr_smoke.sh`; the full layered gate is the main/schedule/manual triage surface.

5. **Artifacts**
   - Smoke artifacts remain under `target/`.
   - QUAST outputs are written under `target/quast-phase2-synthetic/` when the optional row runs. Uploading those artifacts from nightly CI is allowed but not mandatory in this ADR.

## Consequences

Adding a benchmark row, CLI flag, or benchmark script now requires updating the matrix or capabilities page. This is deliberate: product claims should fail fast when they are not represented in the repo's runnable contract.
