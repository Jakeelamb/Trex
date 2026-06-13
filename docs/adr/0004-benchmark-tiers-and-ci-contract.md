# ADR 0004 — Benchmark tiers and CI contract

## Status

Accepted.

## Context

Trex needs architecture and product claims to have local or CI pass/fail loops. The repository already has Phase-1 smoke scripts, Phase-2 Illumina synthetic scripts, and a benchmark matrix file, but the matrix must be enforced and CI tiers must be explicit so benchmark scope does not drift through prose alone.

## Decisions

1. **CI tiers**
   - Pull requests run `cargo run -p xtask -- gate --tier pr` on MSRV, stable, and nightly Rust. The Rust gate owns validators, Clippy, the PR benchmark artifact, and `scripts/pr_smoke.sh`.
   - Pushes to `main` / `master`, tags, schedules, and manual dispatch run the full `scripts/phase2_illumina_benchmark_gate.sh` after installing minimap2.
   - QUAST remains opt-in with `TREX_RUN_QUAST=1`; it is not required for the default PR or main gate until the tool installation and artifact policy are made mandatory.

2. **Matrix-as-code**
   - `tools/benchmark_matrix.toml` is a schema-validated governed matrix.
   - Every row names its minimum CI tier, fixtures, scripts, and optional artifacts or external tools.
   - `cargo run -p xtask -- validate-matrix` fails CI when a row omits required fields, points at missing fixtures/scripts, drifts from declared fixture SHA-256 digests, or declares a direct Trex row whose `out_dir` and CLI `--out-dir` disagree.

3. **Capabilities sync**
   - `docs/CAPABILITIES.md` is the single operator capability table for CLI flags, outputs, checkpoints, benchmark scripts, stability, and deferred work.
   - `cargo run -p xtask -- validate-capabilities` fails CI when current `trex illumina assemble` flags or benchmark scripts are missing from that page.

4. **Layer failure semantics**
   - `scripts/phase2_illumina_benchmark_gate.sh` uses distinct exit codes for its layers: 10 for Phase-1 gate, 20 for Phase-2 reference, 30 for graph summaries, 40 for haplotype metrics, and 50 for optional QUAST.
   - CI still runs most PR layers through `scripts/pr_smoke.sh`; the full layered gate is the main/schedule/manual triage surface.

5. **Artifacts**
   - Smoke artifacts remain under `target/`.
   - `cargo run -p xtask -- bench` writes machine-readable row reports under `target/benchmarks/`, including wall-clock time, script exit codes, max RSS where GNU time is available, observed Trex counters, FASTQ/FASTA/GFA assembly metrics, stage timing summaries, and declared base plus tier-specific artifact sizes.
   - Matrix rows may declare `[rows.trex]` for direct release `trex` executions when the row should measure the product binary rather than only shell script layers. The first direct row is the nightly/manual PhiX174 real-reference micro benchmark.
   - QUAST outputs are written under `target/quast-phase2-synthetic/` when the optional row runs. Uploading those artifacts from nightly CI is allowed but not mandatory in this ADR.

## Consequences

Adding a benchmark row, CLI flag, or benchmark script now requires updating the matrix or capabilities page. This is deliberate: product claims should fail fast when they are not represented in the repo's runnable contract.
