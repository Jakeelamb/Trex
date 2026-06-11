# Trex Capabilities

This page is the operator-facing capability matrix. It is validated by `cargo run -p xtask -- validate-capabilities` so CLI flags and benchmark scripts cannot drift out of sight when code changes.

## Product Surface

| Capability | Phase-1 default | Phase-2 Illumina --diploid | Future / deferred |
|------------|-----------------|-----------------------------|-------------------|
| Read technology | Illumina single-end or paired-end via `--r1` and optional `--r2`; FASTQ/FASTA, including gzip by suffix. | Same Illumina ingest path; no separate diploid subcommand. | Long-read, hybrid, HiFi, CLR, and ONT paths remain deferred by `CONTEXT.md` and ADRs. |
| Counting and trust | `--kmer-size` / `-k`, `--trusted-threshold` / `-T`; canonical *k*-mer counts and one global trusted threshold. | Counting and trusted membership are unchanged before the graph exists. | PE-informed trust, local thresholds, and multi-*k* ladders need a new scope decision. |
| Simplification | Tip and bounded diamond simplification with `--max-tip-bases`, `--tip-max-multiplicity`, `--max-bubble-vertices`, and `--max-bubble-internal-bases`. | Near-balanced diamond branches are retained rather than collapsed. | Complex repeats, long bubbles, tangles, and distance-aware bubble surgery remain future work. |
| Mate usage | Paired reads are validated and concatenated after parity checks; no default mate-derived graph edits. | `--insert-mean-bp` plus paired input enables conservative boosts on existing DBG edges; `--insert-stddev-bp` is stored in mode identity. | Mate-derived new edges, scaffold gaps, and distance-sensitive surgery are out of scope. |
| Outputs | `--out-dir`, `--unitigs-fasta`, `--contigs-fasta`, and `--gfa`; defaults are `unitigs.fa`, `contigs.fa`, and `graph.gfa`; `-` writes selected artifacts to stdout. | Same files, with primary FASTA collapse and GFA `XX:Z:trex-phase2-illumina`, optional `L` rows, `ctg...` `P` rows, and `p2h...` mirror paths where representable. | GFA 2 and mandatory dual haplotype FASTA are deferred. |
| Checkpoints | `--checkpoint-root`, `--resume`, `--no-resume`, `--strict-checkpoints`, and `--no-strict-checkpoints`; graph/export reload by documented identity. | Graph checkpoint identity includes diploid mode, paired input, insert priors, and mate-bridge version. | Cross-mode resume that silently reuses stale graph shape is not allowed. |
| Configuration | `--config` accepts flat assemble keys or an `[assemble]` table; CLI flags override config fields. | `[assemble.diploid]` can set `enabled`, `insert_mean_bp`, and `insert_stddev_bp`. | Additional mode surfaces need README, `CONTEXT.md`, and validator updates. |

## Benchmark And CI Matrix

| Layer | Script | CI tier | Pass/fail contract |
|-------|--------|---------|--------------------|
| Contract validation | `cargo run -p xtask -- validate` | PR, main, tags, schedule, manual | Fails when `tools/benchmark_matrix.toml` lacks required row fields, references missing fixtures/scripts, fixture digests drift from `tools/manifest.toml`, `tools/benchmark_data.toml` is malformed, or this page omits current CLI flags/scripts. |
| Rust warning gate | `cargo clippy --workspace --all-features -- -D warnings` | PR and above | Keeps core crates, CLI, and Rust automation warning-clean across `1.74.0`, `stable`, and `nightly`. |
| Rust PR gate | `cargo run -p xtask -- gate --tier pr` | PR and above | Runs validators, Clippy, the PR benchmark artifact, and `scripts/pr_smoke.sh` through the Rust automation entrypoint used by CI. |
| Matrix benchmark artifact | `cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json` | PR and above | Runs matrix scripts and direct Trex rows for a tier and writes JSON with wall time, exit codes, GNU-time max RSS when available, observed Trex counters, assembly metrics, and artifact sizes; CI uploads the JSON as a workflow artifact. |
| PR smoke | `scripts/pr_smoke.sh` | PR and above | Runs `scripts/ref_free_smoke.sh`, Phase-2 fixture checks, graph summaries, haplotype metrics, and `cargo test --workspace --all-features -q`. |
| Phase-1 reference-free golden | `scripts/ref_free_smoke.sh` | PR and above | Assembles `fixtures/tiny.fq` and byte-compares `contigs.fa`, `unitigs.fa`, and `graph.gfa` against `fixtures/expected/ref_free_smoke/`. |
| Phase-1 full benchmark gate | `scripts/benchmark_gate.sh`, `scripts/reference_smoke.sh` | main, tags, schedule, manual | Runs the reference-free golden, then checks contigs against `fixtures/tiny_ref.fa` with minimap2 when installed or substring fallback. |
| Phase-2 diploid reference layer | `scripts/phase2_illumina_diploid_reference_layer.sh` | PR and above | Verifies SHA-256 digests from `tools/manifest.toml`, parent/read consistency, and optional minimap2 PAF sanity. |
| Phase-2 graph summaries | `scripts/phase2_illumina_graph_summaries.sh` | PR and above | Runs `trex illumina assemble --diploid`, checks primary FASTA stats, GFA record counts, and the Phase-2 GFA header tag. |
| Phase-2 haplotype metrics | `scripts/phase2_illumina_haplotype_metrics.sh` | PR and above | Compares emitted `contigs.fa` against both synthetic parents using best-parent Hamming-style checks. |
| Layered Phase-2 gate | `scripts/phase2_illumina_benchmark_gate.sh` | main, tags, schedule, manual | Runs Phase-1 gate, diploid reference layer, graph summaries, haplotype metrics, and optional QUAST with layer-specific exit codes. |
| PhiX174 real-reference micro row | `cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json` | nightly, manual | Runs release `trex illumina assemble` against `fixtures/phix174/reads.fq` over pinned RefSeq `NC_001422.1`, recording reads, candidate/unique/trusted k-mers, FASTA/GFA sizes, wall time, and max RSS. |
| Biological data fetch | `cargo run -p xtask -- fetch-data` | manual/local | Prepares ignored `data/benchmarks/` subsets from ENA source files declared in `tools/benchmark_data.toml`, then verifies prepared SHA-256s where pinned. |
| Larger bacterial row | `cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json` | manual/local | Runs release Trex on the bounded E. coli MG1655 SRR001666 paired-end subset. The full source is 7,047,668 paired spots / 507,432,096 bases. |
| True diploid eukaryotic row | `cargo run -p xtask -- bench --tier manual --row yeast_btt_err1308583_diploid_1k_pairs --out target/benchmarks/yeast-btt.json` | manual/local | Runs release Trex with `--diploid` on the bounded S. cerevisiae BTT / ERR1308583 paired-end subset. The ploidy table marks BTT accession 1308583 as euploid diploid; the full source is 14,550,715 paired spots / 2,870,913,582 bases. |
| Profiling evidence | `docs/PROFILING.md` plus `target/profiles/` artifacts | manual/local | Records time/RSS/flamegraph commands, current hot symbols, and biological-row blockers without committing raw profiler output. |
| Optional QUAST row | `scripts/reference_quast.sh` | opt-in local/manual/nightly | Runs QUAST or MetaQUAST when `TREX_RUN_QUAST=1` and the tool is installed; `TREX_QUAST_MIN_CONTIG` and `TREX_QUAST_MIN_ALIGNMENT` tune smoke-scale thresholds. |

`tools/benchmark_matrix.toml` is the governed row list. Rows must name their minimum CI tier, fixtures, scripts or direct Trex invocation, optional tools, and artifact paths so adding a biological or larger row is a schema change instead of prose. External rows declare `external_data` and are backed by `tools/benchmark_data.toml`; their ignored `data/` fixtures are verified when present but are not required for default CI. Base `artifacts` are reported for every tier that runs the row; `pr_artifacts`, `main_artifacts`, `nightly_artifacts`, and `manual_artifacts` are reported only for that tier.

## Rust Automation

`xtask` owns repo automation that needs to stay portable and reviewable:

```bash
cargo run -p xtask -- validate
cargo run -p xtask -- validate-matrix
cargo run -p xtask -- validate-capabilities
cargo run -p xtask -- validate-data
cargo run -p xtask -- fetch-data
cargo run -p xtask -- gate --tier pr
cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json
cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json
cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json
cargo run -p xtask -- generate-reads --reference fixtures/phix174/reference.fa --out fixtures/phix174/reads.fq --read-len 150 --step 50 --circular
```

The benchmark artifact is intentionally separate from biological quality claims. It proves the governed row ran and records timing/resource metadata; row scripts and direct Trex row metrics decide whether assembly output passed that row's correctness contract.

## Layer Exit Codes

`scripts/phase2_illumina_benchmark_gate.sh` returns distinct non-zero codes for local and CI triage:

| Exit code | Failed layer |
|-----------|--------------|
| 10 | Phase-1 benchmark gate |
| 20 | Phase-2 diploid reference layer |
| 30 | Phase-2 graph summaries |
| 40 | Phase-2 haplotype metrics |
| 50 | Optional QUAST layer |
