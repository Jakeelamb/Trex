# Trex

Rust genome assembler focused on the **Phase-2 Illumina endgame**, with Phase-1 kept as a non-relaxing compatibility and benchmark layer. Product language and policies live in [`CONTEXT.md`](CONTEXT.md). Pipeline architecture and diagrams live in [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Layout

| Path | Role |
|------|------|
| [`trex/`](trex/) | Sync **`trex`** library: FASTQ/FASTA ingest (+ gzip), preprocess, canonical *k*-mer counts (parallel sort by default), DBG build + tip/diamond bubble simplification, unitigs/contigs, GFA 1.0 + FASTA export, checkpoints |
| [`trex-cli/`](trex-cli/) | Async **`trex-cli`** (`trex` binary): Tokio + `spawn_blocking` into the library |
| [`xtask/`](xtask/) | Rust repository automation: matrix/capability validators and benchmark artifact runner |
| [`fixtures/`](fixtures/) | Phase-1 smoke (`tiny.fq`, `tiny_ref.fa`, `expected/ref_free_smoke/`) + **Phase-2 Illumina** synthetic two-parent [`phase2_synthetic/`](fixtures/phase2_synthetic/) |
| [`scripts/`](scripts/) | Benchmark and smoke scripts; see [`scripts/README.md`](scripts/README.md) (`benchmark_gate.sh`, `phase2_illumina_benchmark_gate.sh`, …) |
| [`fuzz/`](fuzz/) | `cargo-fuzz` harness (`parse_fastq`); see [`fuzz/README.md`](fuzz/README.md) |
| [`tools/manifest.toml`](tools/manifest.toml) | Pinned external tools (e.g. **minimap2**) for reference benchmarks |
| [`docs/CAPABILITIES.md`](docs/CAPABILITIES.md) | Operator capability matrix: CLI flags, outputs, checkpoints, CI tiers, scripts, and deferred work |

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
cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json
```

`ref_free_smoke.sh` writes under `target/ref-free-smoke/` and checks byte-identical `contigs.fa`, `unitigs.fa`, and `graph.gfa` against [`fixtures/expected/ref_free_smoke/`](fixtures/expected/ref_free_smoke/). See [`fixtures/README.md`](fixtures/README.md).

`phase2_illumina_benchmark_gate.sh` runs **`benchmark_gate.sh`** first, then the synthetic **two-parent** diploid reference layer, graph summaries, haplotype metrics, and optional **QUAST** when `TREX_RUN_QUAST=1` (per **Phase-2 Illumina benchmark gate** in [`CONTEXT.md`](CONTEXT.md)). CI runs the full script on **`main`/`master`**, **tags**, **schedule**, and **workflow_dispatch**; pull requests run [`pr_smoke.sh`](scripts/pr_smoke.sh) without mandatory **minimap2** (see [`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

MSRV is **1.74** (`rust-version` in workspace `Cargo.toml`); CI runs `1.74.0` and `stable`.
