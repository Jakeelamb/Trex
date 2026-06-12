# Trex Roadmap Ledger

This ledger turns the assembler plan into active lanes. It is validated by
`cargo run -p xtask -- validate-development`; the machine-readable framework is
`tools/assembler_framework.toml`, and the top-level target architecture is
[`docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md`](ILLUMINA_ASSEMBLER_BLUEPRINT.md).

| Lane | Status | Current proof surface | Next narrow move |
|------|--------|-----------------------|------------------|
| Quality gates | Active | `xtask bench` JSON, reference k-mer quality, optional QUAST, benchmark scripts. | Extend read-vs-assembly k-mer quality and keep quality JSON machine-readable. |
| Evidence ledger | Active | `trex::evidence`, `evidence.json`, mate-bridge evidence records, report-only endpoint join candidates, mate orientation/distance confidence, and `xtask bench` JSON embedding. | Add conflict-cluster rejection before allowing bridge or scaffold outputs to mutate primary artifacts. |
| Graph IR | Active | `trex::dbg::graph::DbgGraph`, dense internal node-id walk adapter, unitig walks, copy-number/repeat annotations, GFA export, profiling docs. | Move the next graph hot path to accessors/id-backed views only after before/after benchmark evidence. |
| Multi-k selection | Active | Explicit `--kmer-ladder`, independent candidate graph scoring, selected graph output, and `multi_k.json`/benchmark JSON evidence. | Add selected-*k* checkpoint namespacing before allowing resume/checkpoint reuse in ladder mode. |
| Simplification policy | Active | Decision-first tip clipping and bounded diamond bubble records, Phase-2 balanced bubble retention, `simplification.json`. | Add repeat-aware guardrail policy before stronger cleanup behavior changes. |
| Assembly audit | Active | Post-assembly `audit.json` / `audit.tsv` with low read-support, mate-hint, repeat-suspicion findings, and report-only `fragmentation.json` endpoint diagnosis; no FASTA repair. | Add richer read-placement evidence before enabling any Pilon-like correction pass. |
| Diploid semantics | Active | `diploid.json`, report-only parent-specific k-mer classification, GFA parent-evidence tags, and explicit no-full-haplotype-FASTA summary fields. | Add real parental haplotype inputs before claiming haplotype-resolved biological output. |
| Path/scaffold builder | Active | Unitigs, primary contigs, GFA `S`/`L`/representable `P` rows, mate-backed `scaffolds.json` schema v3 with orientation/distance evidence, contig endpoint fragmentation diagnosis, and ranked report-only endpoint joins. | Promote high-confidence endpoint joins into scaffold/path artifacts only after conflict handling is explicit. |
| Benchmark matrix | Active | `tools/benchmark_matrix.toml`, `tools/benchmark_data.toml`, `xtask validate-matrix`, Tier-2 yeast 1k manual/release-candidate artifact row, Tier-3 HG002/GRCh38 chr20 contract row, and interval-filtered HG002/chr20 row. | Add GIAB confident-region and variant-aware scoring before claiming human-slice assembly quality. |
| Literature-derived future adapters | Deferred | `docs/ASSEMBLER_FRAMEWORK.md`, ADR 0001, ADR 0002, literature notes. | Keep ABruijn, Flye, Canu, wtdbg2, hifiasm, Verkko, MEGAHIT, and metaSPAdes ideas behind explicit adapter decisions. |

## Phase Order

1. Quality gates before graph mutation.
2. Evidence records before bridge or scaffold application.
3. Simplification audits before more aggressive cleanup.
4. Profiling before graph IR rewrites.
5. Separate path/scaffold artifacts before changing primary FASTA semantics.
6. Biological scale-out from tiny synthetic rows to bacterial, diploid synthetic, yeast/eukaryotic true diploid, then governed human-slice / GIAB-style Tier 3 rows; metagenomic and long-read or hybrid research rows remain deferred until separately governed.

## Blueprint Implementation Waves

| Wave | Coordinated scope | Output goal | Required verification |
|------|-------------------|-------------|-----------------------|
| Wave A | Commit/checkpoint current sidecar/evidence stack. | Current evidence, annotation, simplification, scaffold, multi-*k*, audit, diploid, fragmentation, endpoint join, and GFA metadata work is reviewable as one checkpoint. | `cargo fmt --all --check`; `cargo clippy --workspace --all-features -- -D warnings`; `cargo test --workspace --all-features`; `cargo run -p xtask -- validate`; `cargo run -p xtask -- gate --tier pr`. |
| Wave B | Mate orientation + distance model. | Bridge candidates carry orientation class, insert-distance estimate, support/conflict counts, and confidence. | Mate orientation unit tests, unchanged primary FASTA smoke, PR gate, and E. coli 10k bridge JSON artifact before promotion. |
| Wave C | Endpoint-join promotion into scaffold paths, still no FASTA mutation. | High-confidence endpoint joins promote into deterministic scaffold/path artifacts and representable GFA paths. | Scaffold determinism tests, GFA path tests, PR gate, E. coli 10k benchmark. |
| Wave D | Multi-*k* plus repeat-aware simplification policy. | Selected-*k* graph scoring and repeat-aware simplification guardrails interact through policy summaries. | Synthetic ladder fixtures, repeat-veto tests, one-*k* golden unchanged, E. coli 10k graph/path benchmark. |
| Wave E | Read-backed polishing/audit escalation. | Audit findings gain enough read/k-mer support detail to classify polishing candidates while remaining report-only by default. | Audit unit tests, read-vs-assembly k-mer deltas, default FASTA identity, yeast diploid audit row coverage. |
| Wave F | Diploid-aware ambiguity output. | Retained alternatives and parent-specific evidence are labeled in JSON/GFA without claiming full haplotype FASTA. | Phase-2 synthetic parent evidence, yeast diploid metrics, GFA tag tests, claim-boundary docs. |

## Claim Discipline

Every roadmap item must be either covered by a unit test, covered by `xtask`
validation, covered by a benchmark artifact, or explicitly marked
deferred/research-only. The claim levels and worker protocol live in
`docs/DEVELOPMENT.md`.

## Wave Completion Evidence

| Wave | Required output goal | Current evidence | Status |
|------|----------------------|------------------|--------|
| Evidence ledger | Typed evidence records are emitted without changing graph topology. | `evidence.json`, mate evidence unit/smoke tests, and benchmark JSON embedding. | Proven for mate-bridge evidence. |
| Graph annotation | Copy-number and repeat annotations are sidecar-only and visible in matrix artifacts. | `annotations.json`, node/unitig annotation tests, and E. coli 10k annotation artifact. | Proven as heuristic annotation. |
| Simplification policy | Tip and bubble passes expose decision records while preserving default behavior. | `simplification.json`, decision/mutation equivalence tests, and PR smoke. | Proven for current tip/diamond passes. |
| Bridge and scaffold artifact | Mate evidence produces deterministic scaffold/path sidecars without mutating primary FASTA. | `scaffolds.json` schema v3, existing-edge scaffold tests, orientation/distance evidence tests, report-only endpoint join tests, and primary-contig smoke regression. | Proven for explicit existing-edge bridges and absent-edge endpoint join ranking; conflict-cluster rejection remains the next Wave B refinement. |
| Multi-k graph selection | Multi-*k* mode is explicit, scores independent candidate graphs, and leaves one-*k* default unchanged. | `multi_k.json`, ladder selection tests, checkpoint rejection test, and matrix embedding. | Proven for select-one mode. |
| Compact graph representation | Dense ids preserve graph behavior and improve graph-walk resource use with benchmark evidence. | `CompactDbgGraph`, dense-id walk tests, `target/benchmarks/ecoli-10k-compact-id-id-stitch.json`, `target/benchmarks/ecoli-10k-gfa-path-index.json`, `target/benchmarks/ecoli-100k-stream-json-sidecars.json`. | Proven for behavior, 100k scaling, and sidecar allocation reduction; peak graph memory remains high and is tracked as the next storage-improvement target. |
| Pilon-like audit | Assembly support warnings are emitted as reports only; FASTA is unchanged. | `audit.json`, `audit.tsv`, audit unit tests, yeast 1k audit artifact. | Proven for low-support and repeat-suspicion audit classes. |
| Diploid graph semantics | Parent evidence and graph tags are report-only; no full haplotype FASTA is claimed. | `diploid.json`, parent-reference tests, synthetic/yeast benchmark artifacts, and GFA `PS` tags. | Proven for report-only parent-specific k-mer semantics. |
