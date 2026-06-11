# Trex scripts

| Script | Purpose |
|--------|---------|
| [`pr_smoke.sh`](pr_smoke.sh) | Local/CI **PR smoke**: `ref_free_smoke.sh`, Phase-2 diploid reference layer, graph summaries, haplotype metrics, and `cargo test --workspace --all-features -q`. |
| [`ref_free_smoke.sh`](ref_free_smoke.sh) | **Phase-1 PR smoke**: assemble `fixtures/tiny.fq`, print reference-free stats on FASTA, compare `contigs.fa` / `unitigs.fa` / `graph.gfa` to `fixtures/expected/ref_free_smoke/`. |
| [`reference_smoke.sh`](reference_smoke.sh) | **Phase-1 benchmark gate** (reference leg): check contigs from `ref_free_smoke` against `fixtures/tiny_ref.fa` via **minimap2** PAF or substring fallback. Expects `target/ref-free-smoke/contigs.fa` (or `REF_FREE_SMOKE_OUT`). |
| [`benchmark_gate.sh`](benchmark_gate.sh) | Full **Phase-1 benchmark gate**: `ref_free_smoke.sh` then `reference_smoke.sh`. |
| [`phase2_illumina_diploid_reference_layer.sh`](phase2_illumina_diploid_reference_layer.sh) | **Phase-2 Illumina** diploid fixture leg: SHA-256 vs [`tools/manifest.toml`](../tools/manifest.toml), parental/read consistency checks, optional minimap2. |
| [`phase2_illumina_graph_summaries.sh`](phase2_illumina_graph_summaries.sh) | **Phase-2 Illumina graph summaries**: `trex illumina assemble --diploid` on `fixtures/phase2_synthetic/reads.fq`, stats on `contigs.fa`, GFA record counts + `trex-phase2-illumina` header tag. |
| [`phase2_illumina_haplotype_metrics.sh`](phase2_illumina_haplotype_metrics.sh) | **Phase-2 Illumina haplotype metrics**: compare `target/phase2-graph-summaries/contigs.fa` to both synthetic parents and print best-parent Hamming-style distances. |
| [`reference_quast.sh`](reference_quast.sh) | Optional **QUAST / MetaQUAST** hook for the Phase-2 synthetic assembly. Runs only when called and QUAST is installed; use `TREX_QUAST_MIN_CONTIG` and `TREX_QUAST_MIN_ALIGNMENT` for smoke-scale thresholds. |
| [`phase2_illumina_benchmark_gate.sh`](phase2_illumina_benchmark_gate.sh) | Layered **Phase-2 Illumina benchmark gate**: Phase-1 gate → diploid reference layer → graph summaries → haplotype metrics → optional QUAST when `TREX_RUN_QUAST=1`. |

## CI tiers

- **Pull request**: `cargo run -p xtask -- gate --tier pr` on MSRV, stable, and nightly Rust. The gate runs `cargo run -p xtask -- validate`, `cargo clippy --workspace --all-features -- -D warnings`, `cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json`, and `pr_smoke.sh`.
- **Main / master / tags / schedule / workflow_dispatch**: install minimap2, then run `phase2_illumina_benchmark_gate.sh`.
- **Nightly / manual benchmark artifact**: `cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json` includes the direct release Trex PhiX174 row under `fixtures/phix174/`.
- **Biological manual rows**: `cargo run -p xtask -- fetch-data` prepares ignored ENA subsets from `tools/benchmark_data.toml`; then run a single row with `cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json` or `--row yeast_btt_err1308583_diploid_1k_pairs`.
- **Optional QUAST**: set `TREX_RUN_QUAST=1` before `phase2_illumina_benchmark_gate.sh`; artifacts land under `target/quast-phase2-synthetic/`.

Layer-specific exit codes from `phase2_illumina_benchmark_gate.sh`:

| Exit code | Layer |
|-----------|-------|
| 10 | Phase-1 `benchmark_gate.sh` |
| 20 | Phase-2 diploid reference layer |
| 30 | Phase-2 graph summaries |
| 40 | Phase-2 haplotype metrics |
| 50 | Optional QUAST |

`xtask bench` writes JSON reports with script exit codes, wall-clock time, GNU-time max RSS when `/usr/bin/time` is available, declared artifact sizes, and direct Trex row metrics when a matrix row declares `[rows.trex]`.

See [`.github/workflows/ci.yml`](../.github/workflows/ci.yml), [`tools/benchmark_matrix.toml`](../tools/benchmark_matrix.toml), and [`docs/CAPABILITIES.md`](../docs/CAPABILITIES.md).
