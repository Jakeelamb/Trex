# Trex Roadmap Ledger

This ledger turns the assembler plan into active lanes. It is validated by
`cargo run -p xtask -- validate-development`; the machine-readable framework is
`tools/assembler_framework.toml`, and the top-level target architecture is
[`docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md`](ILLUMINA_ASSEMBLER_BLUEPRINT.md).

| Lane | Status | Current proof surface | Next narrow move |
|------|--------|-----------------------|------------------|
| Quality gates | Active | `xtask bench` JSON, reference k-mer quality, optional QUAST, benchmark scripts. | Extend read-vs-assembly k-mer quality and keep quality JSON machine-readable. |
| Read-trust diagnostics | Active | `trust.json` schema v1 reports trusted/rejected k-mer strata, multiplicity buckets, trusted fraction, and an explicit report-only correction-candidate surface. | Add read-span support and candidate correction evidence before any corrected-read output exists. |
| Evidence ledger | Active | `trex::evidence`, `evidence.json`, mate-bridge evidence records, conflict-aware endpoint join acceptance/rejection, mate orientation/distance confidence, k-bimer-like constraint IDs/support histograms/blockers, and `xtask bench` JSON embedding. | Add read-placement provenance before allowing any scaffold output to mutate primary artifacts. |
| Graph IR | Active | `trex::dbg::graph::DbgGraph`, dense internal node-id walk adapter, unitig walks, copy-number/repeat annotations, GFA export, profiling docs. | Move the next graph hot path to accessors/id-backed views only after before/after benchmark evidence. |
| Multi-k selection | Active | Explicit `--kmer-ladder`, `--auto-k`, independent candidate graph scoring, selected graph output, selected-*k* checkpoint namespacing, and schema-v2 `multi_k.json`/benchmark JSON evidence with completeness, contiguity, dead-end, branch/tangle, repeat-risk, graph-density, and weighted score terms. | Benchmark selected-*k* checkpoint reuse on biological rows before widening auto-k defaults. |
| Simplification policy | Active | Decision-first tip clipping, bounded diamond bubble records, repeat-aware high-copy diamond retention, Phase-2 balanced bubble retention, `simplification.json`. | Add SPAdes-inspired edit/recompress/reannotate/replan scheduling before stronger cleanup behavior changes. |
| Promotion policy | Active | `trex::illumina::promotion` centralizes endpoint-join promotion thresholds, conflict precedence, and evidence-to-stage decisions before scaffold/path construction. | Extend the same policy seam to audit and polishing candidates before any default sequence repair. |
| Assembly audit | Active | Post-assembly `audit.json` / `audit.tsv` with low read-support, mate-hint, repeat-suspicion findings, and report-only `fragmentation.json` endpoint diagnosis; no FASTA repair. | Add richer read-placement evidence before enabling any Pilon-like correction pass. |
| Diploid semantics | Active | `diploid.json`, report-only parent-specific k-mer classification, GFA parent-evidence tags, and explicit no-full-haplotype-FASTA summary fields. | Add real parental haplotype inputs before claiming haplotype-resolved biological output. |
| Path/scaffold builder | Active | Unitigs, primary contigs, GFA `S`/`L`/representable `P` rows, mate-backed `scaffolds.json` schema v6 with promotion-policy snapshots, k-bimer-like constraint provenance, accepted/rejected endpoint joins, orientation/distance evidence, contig endpoint fragmentation diagnosis, high-confidence endpoint sidecar path promotion, separate `scaffolds.fa` output for accepted paths, and tagged `scf...` GFA path projection. | Add read-placement provenance and larger-row scaffold artifact benchmarks before any primary FASTA mutation. |
| Benchmark matrix | Active | `tools/benchmark_matrix.toml`, `tools/benchmark_data.toml`, `xtask validate-matrix`, Tier-2 yeast 1k manual/release-candidate artifact row, Tier-3 HG002/GRCh38 chr20 contract row, interval-filtered HG002/chr20 row, and GIAB confident-region / VCF summary metrics. | Use the truth-region metrics on larger biological rows before claiming human-slice assembly quality. |
| SPAdes architecture lane | Active | [`docs/SPADES_ARCHITECTURE_INSPIRATION.md`](SPADES_ARCHITECTURE_INSPIRATION.md), read-trust diagnostics, multi-*k* evidence, k-bimer-like mate constraints, promotion policy, scaffold paths, repeat-aware simplification notes. | Build read-span trust evidence and iterative simplification as typed artifacts before mutating graph topology. |
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
| Wave D | Multi-*k* plus repeat-aware simplification policy. | Selected-*k* graph scoring, repeat-aware guardrails, and edit/recompress/reannotate/replan summaries interact through policy artifacts. | Synthetic ladder fixtures, repeat-veto tests, recompression/reannotation tests, one-*k* golden unchanged, E. coli 10k graph/path benchmark. |
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
| Simplification policy | Tip and bubble passes expose decision records while preserving default behavior. | `simplification.json` schema v2 with `spades_iterative_v1` scheduler metadata, decision/mutation equivalence tests, scheduler unit test, and PR smoke. | Proven for current tip/diamond passes and topology-delta reporting; non-noop recompression remains future work. |
| Bridge and scaffold artifact | Mate evidence produces deterministic scaffold/path sidecars without mutating primary FASTA. | `scaffolds.json` schema v6, k-bimer-like constraint tests, `scaffolds.fa`, promotion-policy snapshot tests, existing-edge scaffold tests, orientation/distance evidence tests, accepted/rejected endpoint join tests, conflict-cluster rejection tests, and primary-contig smoke regression. | Proven for explicit existing-edge bridges, high-confidence absent-edge sidecar path promotion, separate scaffold FASTA emission, conflict-cluster rejection, and constraint provenance. |
| Multi-k graph selection | Multi-*k* mode is explicit, scores independent candidate graphs, and leaves one-*k* default unchanged. | `multi_k.json` schema v2, ladder/auto selection tests, graph-pressure score term tests, selected-*k* checkpoint namespace resume test, and matrix embedding. | Proven for select-one mode with checkpoint-safe resume and richer score explanations. |
| Read-trust diagnostics | Trusted k-mer thresholds are visible as report-only evidence before correction. | `trust.json`, trust unit test, smoke sidecar assertion, and `xtask bench` embedding. | Proven for count-table strata; read-span correction candidates remain future work. |
| Compact graph representation | Dense ids preserve graph behavior and improve graph-walk resource use with benchmark evidence. | `CompactDbgGraph`, dense-id walk tests, `target/benchmarks/ecoli-10k-compact-id-id-stitch.json`, `target/benchmarks/ecoli-10k-gfa-path-index.json`, `target/benchmarks/ecoli-100k-stream-json-sidecars.json`. | Proven for behavior, 100k scaling, and sidecar allocation reduction; peak graph memory remains high and is tracked as the next storage-improvement target. |
| Pilon-like audit | Assembly support warnings are emitted as reports only; FASTA is unchanged. | `audit.json`, `audit.tsv`, audit unit tests, yeast 1k audit artifact. | Proven for low-support and repeat-suspicion audit classes. |
| Diploid graph semantics | Parent evidence and graph tags are report-only; no full haplotype FASTA is claimed. | `diploid.json`, parent-reference tests, synthetic/yeast benchmark artifacts, and GFA `PS` tags. | Proven for report-only parent-specific k-mer semantics. |
