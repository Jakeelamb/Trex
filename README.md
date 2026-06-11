# Trex

Rust genome assembler focused on the **Phase-2 Illumina endgame**, with Phase-1 kept as a non-relaxing compatibility and benchmark layer. Product language and policies live in [`CONTEXT.md`](CONTEXT.md). Pipeline architecture and diagrams live in [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Layout

| Path | Role |
|------|------|
| [`trex/`](trex/) | Sync **`trex`** library: FASTQ/FASTA ingest (+ gzip), preprocess, canonical *k*-mer counts (parallel sort by default), DBG build + tip/diamond bubble simplification, unitigs/contigs, GFA 1.0 + FASTA export, checkpoints |
| [`trex-cli/`](trex-cli/) | Async **`trex-cli`** (`trex` binary): Tokio + `spawn_blocking` into the library |
| [`xtask/`](xtask/) | Rust repository automation: matrix/capability/data validators, PR gate, read/data fetchers, and benchmark artifact runner |
| [`fixtures/`](fixtures/) | Phase-1 smoke (`tiny.fq`, `tiny_ref.fa`, `expected/ref_free_smoke/`), **Phase-2 Illumina** synthetic two-parent [`phase2_synthetic/`](fixtures/phase2_synthetic/), and real-reference PhiX174 [`phix174/`](fixtures/phix174/) |
| [`scripts/`](scripts/) | Benchmark and smoke scripts; see [`scripts/README.md`](scripts/README.md) (`benchmark_gate.sh`, `phase2_illumina_benchmark_gate.sh`, …) |
| [`fuzz/`](fuzz/) | `cargo-fuzz` harness (`parse_fastq`); see [`fuzz/README.md`](fuzz/README.md) |
| [`tools/manifest.toml`](tools/manifest.toml) | Pinned external tools (e.g. **minimap2**) and prepared fixture digests |
| [`tools/benchmark_data.toml`](tools/benchmark_data.toml) | External biological benchmark catalog: ENA source files, source md5s, prepared subset SHA-256s, and ploidy provenance |
| [`docs/CAPABILITIES.md`](docs/CAPABILITIES.md) | Operator capability matrix: CLI flags, outputs, checkpoints, CI tiers, scripts, and deferred work |
| [`docs/PROFILING.md`](docs/PROFILING.md) | Measured profiling baselines, hot symbols, and biological-row blockers |

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

Outputs default to `unitigs.fa`, `contigs.fa`, and `graph.gfa` under `--out-dir` (override with `--unitigs-fasta`, `--contigs-fasta`, `--gfa`). Use **`-`** as the path for any of those three to write to **stdout** (Phase-1 export sentinel).

Simplification overrides (defaults scale with *k* unless set): `--max-tip-bases`, `--tip-max-multiplicity`, `--max-bubble-vertices`, `--max-bubble-internal-bases`, or TOML `[assemble.simplify]` with the same keys.

**Phase-2 Illumina diploid** (experimental): `--diploid`, optional `--insert-mean-bp` / `--insert-stddev-bp`, or TOML `[assemble.diploid]` with `enabled` and the same insert keys. When enabled, near-balanced diamond bubbles are retained and the GFA `H` line carries `XX:Z:trex-phase2-illumina`.

Flags: `--checkpoint-root`, `--resume`, `--no-resume`, `--strict-checkpoints`, `--no-strict-checkpoints`, `--config` (TOML: optional `[assemble]` table or flat keys; CLI overrides file; `k` may come from the file alone if `--kmer-size` is omitted).

FASTA inputs are detected from path suffixes such as `.fa`, `.fasta`, `.fna` (including `.fa.gz`).

## Smoke

```bash
bash scripts/ref_free_smoke.sh
bash scripts/pr_smoke.sh
bash scripts/benchmark_gate.sh
bash scripts/phase2_illumina_benchmark_gate.sh
cargo run -p xtask -- validate
cargo run -p xtask -- validate-data
cargo run -p xtask -- fetch-data
cargo run -p xtask -- gate --tier pr
cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json
cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json
cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json
```

`ref_free_smoke.sh` writes under `target/ref-free-smoke/` and checks byte-identical `contigs.fa`, `unitigs.fa`, and `graph.gfa` against [`fixtures/expected/ref_free_smoke/`](fixtures/expected/ref_free_smoke/). See [`fixtures/README.md`](fixtures/README.md).

`phase2_illumina_benchmark_gate.sh` runs **`benchmark_gate.sh`** first, then the synthetic **two-parent** diploid reference layer, graph summaries, haplotype metrics, and optional **QUAST** when `TREX_RUN_QUAST=1` (per **Phase-2 Illumina benchmark gate** in [`CONTEXT.md`](CONTEXT.md)). CI runs the full script on **`main`/`master`**, **tags**, **schedule**, and **workflow_dispatch**; pull requests run [`pr_smoke.sh`](scripts/pr_smoke.sh) without mandatory **minimap2** (see [`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

`xtask bench --tier nightly` includes the governed PhiX174 real-reference micro row. It builds `trex-cli` in release mode, runs `trex illumina assemble` on deterministic PhiX reads, and records wall time, max RSS, observed Trex counters, FASTA/GFA artifact sizes, and assembly-size metrics in JSON.

`xtask fetch-data` prepares ignored `data/benchmarks/` subsets from the biological catalog. Current external rows cover **E. coli MG1655 SRR001666** as the larger bacterial rung and **S. cerevisiae BTT / ERR1308583** as the true diploid eukaryotic rung. The source FASTQ files stay external; the catalog records ENA md5s, source read/base counts, prepared subset SHA-256s, and ploidy provenance.

MSRV is **1.74** (`rust-version` in workspace `Cargo.toml`); repo-local development defaults to nightly via [`rust-toolchain.toml`](rust-toolchain.toml), and CI runs `1.74.0`, `stable`, and `nightly`.
