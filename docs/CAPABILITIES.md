# Trex Capabilities

This page is the operator-facing capability matrix. It is validated by `cargo run -p xtask -- validate-capabilities` so CLI flags and benchmark scripts cannot drift out of sight when code changes. The top-level assembler architecture and promotion contract live in [`docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md`](ILLUMINA_ASSEMBLER_BLUEPRINT.md); the live tuning and optimization inventory lives in [`docs/MODULE_MAP.md`](MODULE_MAP.md).

## Product Surface

| Capability | Phase-1 default | Phase-2 Illumina --diploid | Future / deferred |
|------------|-----------------|-----------------------------|-------------------|
| Read technology | Illumina single-end or paired-end via `--r1` and optional `--r2`; FASTQ/FASTA, including gzip by suffix. | Same Illumina ingest path; no separate diploid subcommand. | Long-read, hybrid, HiFi, CLR, and ONT paths remain deferred by `CONTEXT.md` and ADRs. |
| Counting and trust | `--kmer-size` / `-k`, optional explicit `--kmer-ladder`, `--trusted-threshold` / `-T`; canonical *k*-mer counts, one global trusted threshold, and report-only `trust.json` with accepted/rejected strata, multiplicity buckets, trusted fraction, and empty correction-candidate surface. Without `--kmer-ladder`, exactly one requested graph is built. With `--kmer-ladder`, candidate graphs are scored independently and one graph is selected with `multi_k.json` evidence including completeness, contiguity, dead-end, branch/tangle, repeat-risk, graph-density, and weighted score terms. | Counting and trusted membership are unchanged before the graph exists. | PE-informed trust, local thresholds, read correction, and graph merging across *k* values need a new scope decision. |
| Simplification | Tip, bounded diamond, and short low-copy component simplification with `--max-tip-bases`, `--tip-max-multiplicity`, `--max-bubble-vertices`, `--max-bubble-internal-bases`, `--max-low-coverage-component-bases`, and `--low-coverage-component-max-multiplicity`; `simplification.json` records decision arrays plus the `spades_iterative_v1` scheduler pass order, topology deltas, and recompress/reannotation hook status; repeat-aware guardrails retain high-copy diamond branches, and component pruning requires stronger graph components before removing low-copy fragments. | Near-balanced diamond branches are retained rather than collapsed. | Complex repeats, long bubbles, tangles, distance-aware bubble surgery, and non-noop recompression passes remain future work. |
| Mate usage | Paired reads are validated and concatenated after parity checks; no default mate-derived graph edits. | `--insert-mean-bp` plus paired input enables conservative boosts on existing DBG edges and conflict-aware dead-end endpoint join promotion for absent edges; `scaffolds.json` records the active promotion policy, per-candidate promotion stage, accepted/rejected status, orientation, estimated gap, distance confidence, support, conflict counts, and competing endpoint cluster size; `--insert-stddev-bp` is stored in mode identity and contributes to confidence. | Mate-derived primary FASTA edits, scaffold gaps, and distance-sensitive surgery are out of scope. |
| Outputs | `--out-dir`, `--unitigs-fasta`, `--contigs-fasta`, and `--gfa`; defaults are `unitigs.fa`, `contigs.fa`, `graph.gfa`, JSON sidecars including `trust.json`, optional `scaffolds.fa`, contig endpoint `fragmentation.json`, post-assembly `audit.json` / `audit.tsv`, and diploid `diploid.json`; `-` writes selected primary FASTA/GFA artifacts to stdout. | Same files, with primary FASTA collapse and GFA `XX:Z:trex-phase2-illumina`, optional `L` rows, `ctg...` `P` rows, `p2h...` mirror paths where representable, tagged `scf...` scaffold sidecar `P` rows, parent-specific path tags when `--parent1-reference` / `--parent2-reference` are provided, scaffold FASTA sidecars from accepted paths only, fragmentation diagnosis that does not alter contigs, and audit reports that do not repair FASTA. | GFA 2, mandatory dual haplotype FASTA, read-corrected FASTQ output, and audit-driven primary sequence repair are deferred. |
| Checkpoints | `--checkpoint-root`, `--resume`, `--no-resume`, `--strict-checkpoints`, and `--no-strict-checkpoints`; graph/export reload by documented identity. Explicit multi-*k* and `--auto-k` modes keep preprocess checkpoints at the root and store selected graph checkpoints under `selected-k-<k>/`. | Graph checkpoint identity includes diploid mode, paired input, insert priors, and mate-bridge version. | Cross-mode resume that silently reuses stale graph shape is not allowed. |
| Configuration | `--config` accepts flat assemble keys or an `[assemble]` table; CLI flags override config fields; `k_ladder` mirrors `--kmer-ladder`. | `[assemble.diploid]` can set `enabled`, `insert_mean_bp`, `insert_stddev_bp`, `parent1_reference`, and `parent2_reference`. | Additional mode surfaces need README, `CONTEXT.md`, and validator updates. |

## Benchmark And CI Matrix

| Layer | Script | CI tier | Pass/fail contract |
|-------|--------|---------|--------------------|
| Contract validation | `cargo run -p xtask -- validate` | PR, main, tags, schedule, manual | Fails when `tools/benchmark_matrix.toml` lacks required row fields, product claims, claim levels, reference availability, artifact policy, references missing fixtures/scripts, fixture digests drift from `tools/manifest.toml`, `tools/benchmark_data.toml` is malformed, `tools/assembler_framework.toml` references missing papers/code/docs, development protocol docs drift, or this page omits current CLI flags/scripts. |
| Rust warning gate | `cargo clippy --workspace --all-features -- -D warnings` | PR and above | Keeps core crates, CLI, and Rust automation warning-clean across `1.74.0`, `stable`, and `nightly`. |
| Rust PR gate | `cargo run -p xtask -- gate --tier pr` | PR and above | Runs validators, Clippy, the PR benchmark artifact, and `scripts/pr_smoke.sh` through the Rust automation entrypoint used by CI. |
| Matrix benchmark artifact | `cargo run -p xtask -- bench --tier pr --out target/benchmarks/pr.json` | PR and above | Runs matrix scripts and direct Trex rows for a tier and writes JSON with row claim metadata, required/optional tool availability, layer status, wall time, exit codes, GNU-time max RSS when available, observed Trex counters, assembly metrics, typed `evidence.json`, `trust.json`, `annotations.json`, `simplification.json`, `scaffolds.json`, `multi_k.json`, `fragmentation.json`, `audit.json`, and `diploid.json` summaries when present, primary contig/unitig and scaffold FASTA stats, parent-specific diploid reference k-mer quality where declared, GIAB-style confident-region / variant truth summaries when BED/VCF truth is declared, and artifact sizes. Exact read-vs-assembly and full-reference k-mer quality are emitted for bounded rows and skipped with a `metric_notes` explanation above `TREX_XTASK_READ_KMER_QUALITY_MAX` candidate k-mers or `TREX_XTASK_REFERENCE_KMER_QUALITY_MAX_BASES` reference bases unless `TREX_XTASK_FULL_READ_KMER_QUALITY=1` or `TREX_XTASK_FULL_REFERENCE_KMER_QUALITY=1` is set; CI uploads the JSON as a workflow artifact. |
| PR smoke | `scripts/pr_smoke.sh` | PR and above | Runs `scripts/ref_free_smoke.sh`, Phase-2 fixture checks, graph summaries, haplotype metrics, and `cargo test --workspace --all-features -q`. |
| Phase-1 reference-free golden | `scripts/ref_free_smoke.sh` | PR and above | Assembles `fixtures/tiny.fq` and byte-compares `contigs.fa`, `unitigs.fa`, and `graph.gfa` against `fixtures/expected/ref_free_smoke/`. |
| Phase-1 full benchmark gate | `scripts/benchmark_gate.sh`, `scripts/reference_smoke.sh` | main, tags, schedule, manual | Runs the reference-free golden, then checks contigs against `fixtures/tiny_ref.fa` with minimap2 when installed or substring fallback. |
| Phase-2 diploid reference layer | `scripts/phase2_illumina_diploid_reference_layer.sh` | PR and above | Verifies SHA-256 digests from `tools/manifest.toml`, parent/read consistency, and optional minimap2 PAF sanity. |
| Phase-2 graph summaries | `scripts/phase2_illumina_graph_summaries.sh` | PR and above | Runs `trex illumina assemble --diploid`, checks primary FASTA stats, GFA record counts, and the Phase-2 GFA header tag. |
| Phase-2 haplotype metrics | `scripts/phase2_illumina_haplotype_metrics.sh` | PR and above | Compares emitted `contigs.fa` against both synthetic parents using best-parent Hamming-style checks. |
| Layered Phase-2 gate | `scripts/phase2_illumina_benchmark_gate.sh` | main, tags, schedule, manual | Runs Phase-1 gate, diploid reference layer, graph summaries, haplotype metrics, and optional QUAST with layer-specific exit codes. |
| PhiX174 real-reference micro row | `cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json` | nightly, manual | Runs release `trex illumina assemble` against `fixtures/phix174/reads.fq` over pinned RefSeq `NC_001422.1`, recording reads, candidate/unique/trusted k-mers, FASTA/GFA sizes, read-vs-assembly k-mer quality, reference *k*-mer quality, wall time, and max RSS. |
| Automatic k selection | `trex illumina assemble --auto-k --r1 reads.fq --out-dir run` | manual/local | Derives a deterministic odd-*k* ladder from retained read length, scores it through the multi-*k* selector, writes schema-v2 `multi_k.json` with graph-pressure score terms, and logs the selected *k* for auditability. |
| Biological data fetch | `cargo run -p xtask -- fetch-data` | manual/local | Prepares ignored `data/benchmarks/` FASTQ subsets, public reference FASTAs, derived reference slices, interval-paired FASTQs from indexed BAMs, and local Benchmarks repo symlinks declared in `tools/benchmark_data.toml`, then verifies prepared SHA-256s where pinned. |
| Larger bacterial row | `cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_1k_pairs --out target/benchmarks/ecoli.json` | manual/local | Runs release Trex on the bounded E. coli MG1655 SRR001666 paired-end subset. The full source is 7,047,668 paired spots / 507,432,096 bases. |
| Larger bacterial scale-up row | `cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_10k_pairs --out target/benchmarks/ecoli-10k.json` | manual/local | Runs release Trex on 10,000 E. coli MG1655 read pairs with the same pinned RefSeq reference-quality scoring. |
| Higher-depth bacterial profiling row | `cargo run -p xtask -- bench --tier manual --row ecoli_mg1655_srr001666_100k_pairs --out target/benchmarks/ecoli-100k.json` | manual/local | Runs release Trex on 100,000 E. coli MG1655 read pairs. This is the current scaling sentinel for DBG walk and memory behavior. |
| True diploid eukaryotic row | `cargo run -p xtask -- bench --tier manual --row yeast_btt_err1308583_diploid_1k_pairs --out target/benchmarks/yeast-btt.json` | manual/local | Runs release Trex with `--diploid` on the bounded S. cerevisiae BTT / ERR1308583 paired-end subset. The ploidy table marks BTT accession 1308583 as euploid diploid; the full source is 14,550,715 paired spots / 2,870,913,582 bases. |
| Higher-depth diploid eukaryotic profiling row | `cargo run -p xtask -- bench --tier manual --row yeast_btt_err1308583_diploid_10k_pairs --out target/benchmarks/yeast-btt-10k.json` | manual/local | Runs release Trex with `--diploid` on 10,000 BTT read pairs. This is the next eukaryotic scaling rung before full-run claims. |
| Drosophila Benchmarks comparison rows | `cargo run -p xtask -- bench --tier manual --row drosophila_illumina_pe_1m_benchmarks_comparison --out target/benchmarks/drosophila-illumina-1m.json` | manual/local | Runs release Trex with `--diploid` on the same local `/home/jake/Projects/Benchmarks` Drosophila melanogaster DRR001444 Illumina rungs used for SPAdes, MEGAHIT, Velvet, and Tadpole comparison work. Rows exist for 1M, 5M, 25M, and 50M read pairs and record `comparison_threads = 8` to match the local assembler bakeoff. |
| Tier-3 HG002 chr20 contract row | `cargo run -p xtask -- bench --tier manual --row hg002_giab_grch38_chr20_5k_pairs --out target/benchmarks/hg002-chr20.json` | manual/local | Runs release Trex with `--diploid` on the governed HG002 WGS-prefix row. This proves source/truth/artifact plumbing but is not a human-readiness claim. |
| Tier-3 HG002 chr20 interval row | `cargo run -p xtask -- bench --tier manual --row hg002_giab_grch38_chr20_interval_1k_pairs --out target/benchmarks/hg002-chr20-interval.json` | manual/local | Runs release Trex with `--diploid` on 1,000 paired reads extracted from the GIAB/NHGRI GRCh38 chr20 Novoalign BAM over chr20:10000000-10010000, using the matching 10,001 bp reference slice, GIAB confident-region / VCF truth summaries, and optional QUAST. |
| Profiling evidence | `docs/PROFILING.md` plus `target/profiles/` artifacts | manual/local | Records time/RSS/flamegraph commands, current hot symbols, and biological-row blockers without committing raw profiler output. |
| Optional QUAST row | `scripts/reference_quast.sh` | opt-in local/manual/nightly | Runs QUAST or MetaQUAST when `TREX_RUN_QUAST=1` and the tool is installed; `TREX_QUAST_REF`, `TREX_QUAST_ASM`, and `TREX_QUAST_OUT` target a specific assembly, while `TREX_QUAST_MIN_CONTIG` and `TREX_QUAST_MIN_ALIGNMENT` tune smoke-scale thresholds. |
| Literature review archive | `literature/sources.md` | manual/local | Tracks archived PDFs, source-only papers, and review targets for OLC/DBG, long-read, diploid/T2T, metagenome, polishing, and evaluation design work. |
| Assembler framework validation | `cargo run -p xtask -- validate-framework` | PR and above through `validate` | Fails when the literature-informed framework points at missing papers, missing implementation files, empty promotion stages, unsupported module phase/status values, missing blueprint promotion stages, or framework modules not represented in the architecture docs. |

`tools/benchmark_matrix.toml` is the governed row list. Rows must name their minimum CI tier, product claim, claim level, reference availability, artifact policy, fixtures, scripts or direct Trex invocation, required tools, optional tools, and artifact paths so adding a biological or larger row is a schema change instead of prose. External rows declare `external_data` and are backed by `tools/benchmark_data.toml`; their ignored `data/` fixtures are verified when present but are not required for default CI. Base `artifacts` are reported for every tier that runs the row; `pr_artifacts`, `main_artifacts`, `nightly_artifacts`, and `manual_artifacts` are reported only for that tier.

## Rust Automation

`xtask` owns repo automation that needs to stay portable and reviewable:

```bash
cargo run -p xtask -- validate
cargo run -p xtask -- validate-matrix
cargo run -p xtask -- validate-capabilities
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
cargo run -p xtask -- generate-reads --reference fixtures/phix174/reference.fa --out fixtures/phix174/reads.fq --read-len 150 --step 50 --circular
```

The benchmark artifact is intentionally separate from biological quality claims. It proves the governed row ran and records timing/resource metadata; row scripts and direct Trex row metrics decide whether assembly output passed that row's correctness contract.

For fast local development on expensive rows, set `TREX_XTASK_BENCH_RESUME=1` on
`cargo run -p xtask -- bench ...`. `xtask` then preserves the Trex output
directory, injects `--checkpoint-root <out_dir>/checkpoints --resume` unless the
row already declares checkpoint flags, and records `checkpoint_resume` in the
JSON artifact. This reuses Trex stage outputs under
`preprocess/reads.jsonl`, `counts/kmer_counts.json`,
`graph/simplified_dbg.json`, and `export/sequences.json`. Leave the variable
unset for clean comparison artifacts.

## Layer Exit Codes

`scripts/phase2_illumina_benchmark_gate.sh` returns distinct non-zero codes for local and CI triage:

| Exit code | Failed layer |
|-----------|--------------|
| 10 | Phase-1 benchmark gate |
| 20 | Phase-2 diploid reference layer |
| 30 | Phase-2 graph summaries |
| 40 | Phase-2 haplotype metrics |
| 50 | Optional QUAST layer |
