# Trex Assembler Framework

This framework translates the literature archive into Trex modules, seams, and benchmark gates. It does not reopen long-read or hybrid scope by itself; ADR 0001 and ADR 0002 still keep active implementation on the Phase-2 Illumina endgame.

The machine-readable version is [`tools/assembler_framework.toml`](../tools/assembler_framework.toml), validated by `cargo run -p xtask -- validate-framework` and included in `cargo run -p xtask -- validate`.

## Literature-Derived Principles

| Principle | Source papers | Trex consequence |
|-----------|---------------|------------------|
| Evidence before mutation | Unicycler, Pilon, Canu, ABruijn | Mate links, bridges, corrections, and future long-read paths should be scored evidence records before they edit the graph or FASTA. |
| Preserve ambiguity structurally | Canu, Unicycler, Li/Durbin | GFA and path artifacts are first-class; a single FASTA must not hide unresolved repeats or haplotype ambiguity. |
| Coverage is contextual | metaSPAdes, Unicycler, Merqury | Simplification policies must not assume one uniform haploid coverage model. |
| K-mer set validation is core | Merqury, Li/Durbin | Trex should own read-vs-assembly k-mer quality gates, not rely only on reference alignment. |
| Graph IR follows evidence scale | ABruijn, Canu, Li/Durbin | Illumina DBG is current; long-read/sparse/repeat graph tracks need explicit adapters. |
| Correctness needs artifacts | Pilon, Merqury, QUAST practice | Every quality claim needs JSON/TSV/GFA/changes artifacts that can be diffed and bisected. |

## Deep Modules To Build

### 1. Evidence Ledger (`evidence_ledger`)

Interface:

- Accept typed evidence: read k-mers, paired-end adjacency, bridge candidates, assembly k-mer quality, future long-read path evidence.
- Return scored, auditable records.
- Never directly mutate counts, graph topology, or FASTA.

Current adapters:

- Trusted k-mer counts in `trex::illumina::counts`.
- Conservative existing-edge mate boost in `trex::illumina::mate`.
- Reference-quality metrics in `xtask`.

Next work:

- Add a bridge-candidate record type before allowing mate-derived new edges.
- Emit evidence summaries into benchmark JSON.

### 2. Graph IR (`graph_ir`)

Interface:

- Own graph vertices, edges, multiplicity, path compression, and graph invariants.
- Hide storage choices behind graph operations.
- Preserve enough metadata for unitigs, contigs, paths, and quality gates.

Current adapter:

- Canonical Illumina DBG in `trex::dbg::graph::DbgGraph`.

Next work:

- Add copy-number/repeat annotations derived from node and edge multiplicity.
- Investigate compact/sparse storage once the 100k E. coli row is fast enough to profile cleanly.

### 3. Simplification Policy (`simplification_policy`)

Interface:

- Take graph plus evidence summaries.
- Return graph edits plus an audit trail.
- Encode mode-specific thresholds explicitly.

Current adapters:

- Tip clipping.
- Bounded diamond bubble collapse.
- Diploid near-balanced diamond retention.

Next work:

- Split simplification decisions from mutation.
- Add repeat/copy-number-aware decisions before adding aggressive bridge application.

### 4. Path And Scaffold Builder (`path_scaffold_builder`)

Interface:

- Convert graph paths and evidence-backed bridges into outputs.
- Keep contigs, scaffolds, unitigs, and unresolved graph paths distinct.
- Maintain GFA/FASTA identity consistency.

Current adapters:

- Unitig extraction.
- Reference contig walks.
- GFA `S` / `L` / representable `P` export.

Next work:

- Add a scaffold artifact only after bridge evidence exists.
- Keep sequence-only postprocessing out of `contigs.fa` unless path metadata can be updated consistently.

### 5. Assembly Quality Module (`assembly_quality`)

Interface:

- Compute reference-free and reference-backed quality metrics.
- Emit machine-readable artifacts.
- Keep accuracy, completeness, phasing, and continuity separate.

Current adapters:

- `xtask bench` FASTQ/FASTA/GFA metrics.
- Reference k-mer containment.
- Optional QUAST.

Next work:

- Add Merqury-style read-vs-assembly k-mer metrics.
- Add parent-specific k-mer metrics for diploid rows.
- Add a report artifact for collapsed-repeat suspicion instead of treating short contigs as merely cosmetic.

### 6. Future Graph Adapters (`future_graph_adapters`)

Interface:

- Convert non-Illumina evidence into graph/path forms without changing the active Illumina contract.

Deferred adapters:

- A-Bruijn / landmark graph for noisy long reads.
- Overlap/string graph for accurate long reads.
- Sparse/minimizer DBG for memory-heavy large genomes.
- Long-range phasing/scaffolding adapters for Hi-C/Pore-C/trio data.

These are tracked as architecture pressure, not current implementation scope.

## Immediate Build Order

1. **Quality module:** add read-vs-assembly k-mer metrics to `xtask bench`.
2. **Evidence ledger:** represent mate/bridge candidates without mutating graph topology.
3. **Simplification audit:** return decisions and counts from tip/bubble passes.
4. **Graph metadata:** add copy-number/repeat annotations.
5. **Scaffold artifact:** emit scaffold paths only when backed by evidence records.

The order is deliberate: without quality gates and evidence records, graph edits are not falsifiable.
