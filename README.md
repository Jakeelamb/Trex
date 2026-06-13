# Trex

Rust genome assembler focused on the **Phase-2 Illumina endgame**, with Phase-1 kept as a non-relaxing compatibility and benchmark layer. Product language and policies live in [`CONTEXT.md`](CONTEXT.md). Pipeline architecture and diagrams live in [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Layout

| Path | Role |
|------|------|
| [`trex/`](trex/) | Sync **`trex`** library: FASTQ/FASTA ingest (+ gzip), preprocess, canonical *k*-mer counts (parallel sort by default), read-trust diagnostics, DBG build + tip/diamond/component simplification, unitigs/contigs, GFA 1.0 + FASTA export, checkpoints |
| [`trex-cli/`](trex-cli/) | Async **`trex-cli`** (`trex` binary): Tokio + `spawn_blocking` into the library |
| [`xtask/`](xtask/) | Rust repository automation: matrix/capability/data validators, PR gate, read/data fetchers, and benchmark artifact runner |
| [`fixtures/`](fixtures/) | Phase-1 smoke (`tiny.fq`, `tiny_ref.fa`, `expected/ref_free_smoke/`), **Phase-2 Illumina** synthetic two-parent [`phase2_synthetic/`](fixtures/phase2_synthetic/), and real-reference PhiX174 [`phix174/`](fixtures/phix174/) |
| [`scripts/`](scripts/) | Benchmark and smoke scripts; see [`scripts/README.md`](scripts/README.md) (`benchmark_gate.sh`, `phase2_illumina_benchmark_gate.sh`, …) |
| [`fuzz/`](fuzz/) | `cargo-fuzz` harness (`parse_fastq`); see [`fuzz/README.md`](fuzz/README.md) |
| [`tools/manifest.toml`](tools/manifest.toml) | Pinned external tools (e.g. **minimap2**) and prepared fixture digests |
| [`tools/benchmark_data.toml`](tools/benchmark_data.toml) | External biological benchmark catalog: ENA source files, source md5s, prepared subset SHA-256s, and ploidy provenance |
| [`docs/CAPABILITIES.md`](docs/CAPABILITIES.md) | Operator capability matrix: CLI flags, outputs, checkpoints, CI tiers, scripts, and deferred work |
| [`docs/ASSEMBLER_FRAMEWORK.md`](docs/ASSEMBLER_FRAMEWORK.md) | Literature-informed module framework for evidence, graph IR, simplification, paths/scaffolds, quality, and deferred graph adapters |
| [`docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md`](docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md) | Target Illumina assembler architecture, promotion policy, data flow, unit-test map, and wave implementation plan |
| [`docs/MODULE_MAP.md`](docs/MODULE_MAP.md) | Live assembler module inventory for tuning, optimization, memory pressure, and code-size cleanup |
| [`docs/DEVELOPMENT.md`](docs/DEVELOPMENT.md) | Orchestrator/worker protocol, claim levels, acceptance checklist, and commit discipline |
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | Active lane ledger for quality gates, evidence, graph IR, simplification, paths/scaffolds, benchmark matrix, and future adapters |
| [`docs/PROFILING.md`](docs/PROFILING.md) | Measured profiling baselines, hot symbols, and biological-row blockers |
| [`literature/`](literature/) | Assembly literature archive and review queue for OLC/DBG, long-read, diploid/T2T, metagenome, polishing, and evaluation papers |

## Build

```bash
cargo build --workspace
cargo test --workspace
```

## CLI

```bash
cargo run -p trex-cli -- illumina assemble --r1 reads.fq --kmer-size 31 --out-dir ./run1
```

Paired-end:

```bash
cargo run -p trex-cli -- illumina assemble --r1 r1.fq --r2 r2.fq --kmer-size 31 -T 2 --out-dir ./run1
```

Outputs default to `unitigs.fa`, `contigs.fa`, `graph.gfa`, typed `evidence.json`, report-only read trust `trust.json`, copy-number/repeat `annotations.json`, decision-first `simplification.json` with scheduled pass/topology metadata, evidence-backed `scaffolds.json`, optional scaffolded sequence sidecar `scaffolds.fa`, contig endpoint `fragmentation.json`, post-assembly audit reports `audit.json` / `audit.tsv`, and diploid evidence `diploid.json` under `--out-dir` (override primary FASTA/GFA paths with `--unitigs-fasta`, `--contigs-fasta`, `--gfa`). `scaffolds.json` schema v6 carries k-bimer-like mate constraint IDs, graph-context endpoints, distance bins, support histograms, and blocker reasons for bridge, endpoint-join, and scaffold-path records. Explicit multi-*k* mode also writes `multi_k.json`. Use **`-`** as the path for any of the primary FASTA/GFA outputs to write to **stdout** (Phase-1 export sentinel).

Explicit multi-*k* graph selection is opt-in with `--kmer-ladder 21,31,41` or TOML `k_ladder = [21, 31, 41]`. Trex builds and scores candidate graphs independently, selects one graph for the normal FASTA/GFA path, and records candidate completeness, contiguity, dead-end, branch/tangle, repeat-risk, graph-density, and weighted score terms in `multi_k.json`. The default remains one-*k*; with `--checkpoint-root`, multi-*k* stores counts, graph, and export checkpoints under a selected-*k* namespace so resume cannot reuse artifacts from a different chosen graph.

Simplification overrides (defaults scale with *k* unless set): `--max-tip-bases`, `--tip-max-multiplicity`, `--max-bubble-vertices`, `--max-bubble-internal-bases`, `--max-low-coverage-component-bases`, `--low-coverage-component-max-multiplicity`, or TOML `[assemble.simplify]` with the same keys. The component pass is SPAdes-inspired cleanup: it removes only short low-copy disconnected components when stronger graph components exist, and records every removal in `simplification.json`.

**Phase-2 Illumina diploid** (experimental): `--diploid`, optional `--insert-mean-bp` / `--insert-stddev-bp`, optional report-only `--parent1-reference` / `--parent2-reference`, or TOML `[assemble.diploid]` with `enabled`, insert keys, `parent1_reference`, and `parent2_reference`. When enabled, near-balanced diamond bubbles are retained and the GFA `H` line carries `XX:Z:trex-phase2-illumina`. Parent references produce parent-specific k-mer evidence in `diploid.json` and GFA path tags; Trex still emits one primary FASTA and does not claim full haplotype FASTA.

Flags: `--checkpoint-root`, `--resume`, `--no-resume`, `--strict-checkpoints`, `--no-strict-checkpoints`, `--kmer-ladder`, `--config` (TOML: optional `[assemble]` table or flat keys; CLI overrides file; `k` may come from the file alone if `--kmer-size` is omitted; a non-empty `k_ladder` can provide the requested baseline when `k` is omitted).

For fast local benchmark iteration, `TREX_XTASK_BENCH_RESUME=1 cargo run -p xtask -- bench ...` preserves the row output directory and injects a checkpoint root under `<out_dir>/checkpoints`; leave it unset for clean comparison runs.

FASTA inputs are detected from path suffixes such as `.fa`, `.fasta`, `.fna` (including `.fa.gz`).

## Smoke

```bash
bash scripts/ref_free_smoke.sh
bash scripts/pr_smoke.sh
bash scripts/benchmark_gate.sh
bash scripts/phase2_illumina_benchmark_gate.sh
cargo run -p xtask -- validate
cargo run -p xtask -- validate-data
cargo run -p xtask -- validate-framework
cargo run -p xtask -- validate-development
cargo run -p xtask -- fetch-data
cargo run -p xtask -- gate --tier pr
cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json
cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json
cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json
cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_10k_pairs --out target/benchmarks/ecoli-10k.json
cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_100k_pairs --out target/benchmarks/ecoli-100k.json
cargo run -p xtask -- bench --tier manual --row yeast_btt_err1308583_diploid_10k_pairs --out target/benchmarks/yeast-btt-10k.json
cargo run -p xtask -- bench --tier manual --row drosophila_illumina_pe_1m_benchmarks_comparison --out target/benchmarks/drosophila-illumina-1m.json
cargo run -p xtask -- bench --tier manual --row drosophila_illumina_pe_5m_benchmarks_comparison --out target/benchmarks/drosophila-illumina-5m.json
cargo run -p xtask -- bench --tier manual --row drosophila_illumina_pe_25m_benchmarks_comparison --out target/benchmarks/drosophila-illumina-25m.json
cargo run -p xtask -- bench --tier manual --row drosophila_illumina_pe_50m_benchmarks_comparison --out target/benchmarks/drosophila-illumina-50m.json
```

`ref_free_smoke.sh` writes under `target/ref-free-smoke/` and checks byte-identical `contigs.fa`, `unitigs.fa`, and `graph.gfa` against [`fixtures/expected/ref_free_smoke/`](fixtures/expected/ref_free_smoke/). See [`fixtures/README.md`](fixtures/README.md).

`phase2_illumina_benchmark_gate.sh` runs **`benchmark_gate.sh`** first, then the synthetic **two-parent** diploid reference layer, graph summaries, haplotype metrics, and optional **QUAST** when `TREX_RUN_QUAST=1` (per **Phase-2 Illumina benchmark gate** in [`CONTEXT.md`](CONTEXT.md)). CI runs the full script on **`main`/`master`**, **tags**, **schedule**, and **workflow_dispatch**; pull requests run [`pr_smoke.sh`](scripts/pr_smoke.sh) without mandatory **minimap2** (see [`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

`xtask bench --tier nightly` includes the governed PhiX174 real-reference micro row. It builds `trex-cli` in release mode, runs `trex illumina assemble` on deterministic PhiX reads, and records wall time, max RSS, observed Trex counters, FASTA/GFA artifact sizes, assembly-size metrics, typed evidence, graph annotation, simplification-decision, scaffold, multi-*k*, post-assembly audit, and diploid parent-evidence summaries when present, read-vs-assembly *k*-mer quality, reference *k*-mer quality metrics, and GIAB-style confident-region / variant truth summaries when the row declares BED/VCF truth files.

`trex illumina assemble --auto-k` derives a deterministic odd-*k* ladder from the shortest retained read, scores candidates with the same multi-*k* selector used by `--kmer-ladder`, writes `multi_k.json`, and logs the selected *k*. Use explicit `--kmer-size` for fixed-*k* reproducibility and `--kmer-ladder` when a benchmark needs a hand-curated candidate set.

`xtask fetch-data` prepares ignored `data/benchmarks/` subsets from the biological catalog. Current external rows cover **E. coli MG1655 SRR001666** at 1k, 10k, and 100k paired-read rungs, **S. cerevisiae BTT / ERR1308583** at 1k and 10k paired-read true diploid eukaryotic rungs, **HG002 / GRCh38 chr20** Tier-3 human-slice rows, and the local `/home/jake/Projects/Benchmarks` **Drosophila melanogaster DRR001444** 1M, 5M, 25M, and 50M-pair Illumina comparison rungs. The source FASTQ/BAM files stay external; the catalog records ENA/GIAB/source metadata, source md5s when published, public reference FASTA URLs, local Benchmarks source paths, derived reference-slice SHA-256s, prepared FASTQ SHA-256s, GIAB truth files, and ploidy provenance. Local Benchmarks entries are symlinked into ignored Trex `data/benchmarks/` paths and SHA-256 checked. Drosophila comparison rows record `comparison_threads = 8` to match the local SPAdes, MEGAHIT, Velvet, and Tadpole runs. Set `TREX_RUN_QUAST=1` on `xtask bench` to run QUAST / MetaQUAST for direct Trex rows with a declared reference.

MSRV is **1.74** (`rust-version` in workspace `Cargo.toml`); repo-local development defaults to nightly via [`rust-toolchain.toml`](rust-toolchain.toml), and CI runs `1.74.0`, `stable`, and `nightly`.
