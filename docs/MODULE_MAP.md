# Trex module map

This page is the live optimization inventory for the Illumina assembler. It lists each module by its main interface, tuning knobs, performance pressure, and current status so future speed, memory, and code-size work has a clear target.

## Pipeline modules

| Order | Module | Main files | Interface | Tuning / policy knobs | Performance pressure | Status |
|-------|--------|------------|-----------|------------------------|----------------------|--------|
| 1 | Operator/config | `trex-cli/src/main.rs` | `trex illumina assemble` flags and config-to-`AssembleParams` mapping. | `--kmer-size`, `--auto-k`, `--kmer-ladder`, `-T`, checkpoints, output paths, diploid selectors. | CLI startup is not the hot path; keep parsing simple and avoid duplicated operator surfaces. | Active |
| 2 | Read ingest | `illumina/fastq.rs`, `illumina/fasta.rs`, `illumina/io.rs`, `illumina/read.rs` | FASTQ/FASTA/gzip readers into `Read` records. | File format, gzip/plain dispatch, empty-read policy. | Allocation per read and gzip I/O throughput. | Active |
| 3 | Preprocess and pair identity | `illumina/preprocess.rs`, `illumina/paired.rs` | N-free ACGT segments, normalized reads, R1/R2 parity. | N policy, case normalization, paired layout. | Read copying and segment iteration. | Active |
| 4 | K-mer counting and trust | `illumina/counts.rs`, `illumina/trust.rs`, `kmer.rs` | Canonical k-mer enumeration, sorted merge, trusted threshold, report-only trust diagnostics. | `k`, global `T`, canonical identity. | Primary CPU and memory pressure on larger datasets; keep trust summaries streaming-friendly. | Active |
| 5 | Multi-k graph ladder | `illumina/multik.rs` | Candidate k ladder scoring and select-one report. | Explicit ladder, `--auto-k` derived ladder, selected-k checkpoint namespace. | Repeated graph builds and count reuse strategy. | Active |
| 6 | DBG construction | `dbg/graph.rs`, `dbg/orient.rs` | Trusted de Bruijn graph plus forward representatives. | Trusted vertices only; graph identity fields. | Adjacency storage, clone pressure, graph build loops. | Active |
| 7 | Mate evidence | `illumina/mate.rs` | K-bimer-like mate constraint records with stable IDs, graph-context endpoints, orientation, distance bins, support histograms, blocker reasons, and existing-edge bridge counters. | Insert mean/stddev, support and conflict counts, blocker classification. | Pair scanning and endpoint lookup. | Active |
| 8 | Promotion policy | `illumina/promotion.rs` | Evidence-to-stage decisions and rejection reasons. | Minimum support, distance confidence, conflict and cluster rejection. | Should stay tiny; optimize for locality and auditability. | Active |
| 9 | Simplification policy | `dbg/simplify.rs` | Tip, diamond, and short low-copy component decisions, graph edits, guardrails, future edit/recompress/reannotate/replan loop. | Tip length/multiplicity, bubble node/base budgets, component span/multiplicity, stronger-component guardrail, diploid/repeat retention. | Graph mutation complexity, recompression cost, and decision-log allocation. | Active |
| 10 | Unitigs and contig walks | `dbg/unitig.rs`, `dbg/walk.rs` | Unitig extraction, sequence stitching, reference contig paths. | Walk tie-break mode and simplified graph invariants. | Path traversal, sequence stitching, clone avoidance. | Active |
| 11 | Phase-2 primary collapse | `illumina/phase2_primary.rs` | Trusted-k-mer multiplicity voting over primary contig FASTA. | Diploid selector, deterministic A/C/G/T tie-break. | Window rescoring over long contigs. | Active |
| 12 | Annotation, fragmentation, audit | `dbg/annotate.rs`, `illumina/fragmentation.rs`, `illumina/audit.rs`, `illumina/diploid.rs` | Repeat labels, endpoint diagnoses, parent evidence, audit reports. | Copy-number class, endpoint reason, parent references, audit thresholds. | Sidecar generation and repeated path scans. | Active |
| 13 | Scaffold/path builder | `illumina/scaffold.rs` | `scaffolds.json` schema v6, accepted scaffold paths, scaffold FASTA records, GFA path projection, and constraint-id propagation through links/joins. | Promotion decisions, policy snapshots, orientation, distance bins, support histograms, overlap trimming, positive-gap `N` padding. | Sequence concatenation and unitig lookup. | Active |
| 14 | Export | `dbg/export.rs` | FASTA and GFA 1.0 writers. | `ctg...`, `p2h...`, `scf...`, `L`/`P` record emission. | Output streaming and path determinism. | Active |
| 15 | Checkpointing | `illumina/checkpoint.rs` | Preprocess/counts/graph/export checkpoint load/write. | Strict manifests, graph identity, selected-k namespace. | Serialization size and stale-stage invalidation. | Active |
| 16 | Benchmark and governance | `xtask/src/main.rs`, `tools/benchmark_matrix.toml` | Matrix validation, benchmark execution, artifact metrics. | CI tier, row metadata, external tools, confident-region inputs. | Benchmark runtime and artifact size accounting. | Active |

## Optimization rule

Do not optimize a module from intuition alone. For each speed or memory change, record the target module, baseline command, before/after artifact, correctness gate, and whether the change is kept or reverted. The current default gates are `cargo test --workspace --all-features -q`, `cargo clippy --workspace --all-features -- -D warnings`, `cargo run -q -p xtask -- validate`, and the relevant `xtask bench` row.
