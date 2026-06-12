# Trex Assembler Framework

This framework translates the literature archive into Trex modules, seams, and benchmark gates. It does not reopen long-read or hybrid scope by itself; ADR 0001 and ADR 0002 still keep active implementation on the Phase-2 Illumina endgame.

The machine-readable contract is [`tools/assembler_framework.toml`](../tools/assembler_framework.toml). The detailed architecture, promotion ladder, data flow, and wave plan live in [`docs/ILLUMINA_ASSEMBLER_BLUEPRINT.md`](ILLUMINA_ASSEMBLER_BLUEPRINT.md). Both are validated by `cargo run -p xtask -- validate-framework` and included in `cargo run -p xtask -- validate`.

## Literature-Derived Principles

| Principle | Source families | Trex consequence |
|-----------|-----------------|------------------|
| Evidence before mutation | Unicycler, Pilon, Canu, ABruijn | Mate links, bridges, corrections, and future long-read paths become scored evidence records before they edit the graph or FASTA. |
| Preserve ambiguity structurally | Canu, Unicycler, Flye-style graph thinking, Li/Durbin | GFA and path artifacts are first-class; a single FASTA must not hide unresolved repeats or haplotype ambiguity. |
| Coverage is contextual | SPAdes, metaSPAdes, Unicycler, Merqury | Simplification policies must not assume one uniform haploid coverage model. Local ratios and neighboring graph context matter before a low-coverage branch is removed. |
| Multi-*k* is evidence, not magic | SPAdes, Unicycler | Candidate graphs should be scored, rejected, selected, and checkpointed as inspectable artifacts before Trex attempts cross-*k* graph merging. |
| Cleanup is iterative | SPAdes source behavior, metaSPAdes | Graph edits must be followed by recompression, reannotation, and replanning when later decisions depend on topology. |
| K-mer set validation is core | Merqury, T2T-era work | Trex should own read-vs-assembly k-mer quality gates, not rely only on reference alignment. |
| Promotion must be explicit | Unicycler, Pilon, Merqury | Report-only candidates, scaffold artifacts, GFA paths, FASTA scaffolds, graph edits, and polishing edits are separate claims. |
| Correctness needs artifacts | Pilon, Merqury, QUAST practice | Every quality claim needs JSON/TSV/GFA/changes artifacts that can be diffed and bisected. |

## Top-Level Modules

| Module | Status | Responsibility | Promotion stages |
|--------|--------|----------------|------------------|
| `read_correction_trust` | Active | Read correction / trusted k-mer model: ingest diagnostics, canonical counts, trusted thresholds, report-only `trust.json`, and future correction candidates. | `report_only_candidate` |
| `multi_k_graph_ladder` | Active | Explicit SPAdes-style multi-*k* candidate graph build, scoring, and select-one behavior. | `report_only_candidate` |
| `repeat_annotation` | Active | Repeat-aware graph annotation: node/unitig multiplicity, endpoint class, and repeat suspicion. | `report_only_candidate` |
| `simplification_policy` | Active | Decision-first simplification for tips, bubbles, repeats, and diploid guardrails. | `report_only_candidate`, `graph_edit` |
| `mate_evidence` | Active | Mate-pair distance/orientation evidence for existing-edge boosts, endpoint joins, conflicts, and bridge ranking. | `report_only_candidate`, `scaffold_artifact` |
| `promotion_policy` | Active | Central policy that converts evidence into behavior through the promotion ladder. | all promotion stages |
| `path_scaffold_builder` | Active | Scaffold/path promotion surface for unitigs, primary contigs, GFA paths, and scaffold artifacts. | `scaffold_artifact`, `gfa_path`, `fasta_scaffold_with_gaps` |
| `polishing_audit` | Active | Pilon-like audit loop that reports low support, mate disagreement, and repeat-collapse suspicion before sequence repair. | `report_only_candidate`, `polishing_edit` |
| `diploid_ambiguity` | Active | Diploid ambiguity handling: retained alternatives, parent-specific evidence, and claim-boundary metrics. | `report_only_candidate`, `gfa_path` |
| `assembly_quality` | Active | Quality gates and benchmark claims with machine-readable artifacts and explicit claim levels. | `report_only_candidate` |
| `future_graph_adapters` | Deferred | A-Bruijn, overlap/string graph, sparse/minimizer DBG, long-range phasing, and non-Illumina adapter pressure. | `report_only_candidate` |

## Promotion Ladder

The promotion stages are:

1. `report_only_candidate`
2. `scaffold_artifact`
3. `gfa_path`
4. `fasta_scaffold_with_gaps`
5. `graph_edit`
6. `polishing_edit`

Later stages require evidence from earlier stages plus the appropriate tests, docs, changelog entry, and benchmark artifact. Primary `contigs.fa` remains conservative until a promotion stage explicitly allows FASTA changes.

## Wave Build Order

1. Wave A: commit/checkpoint the current sidecar/evidence stack.
2. Wave B: add mate orientation and distance confidence.
3. Wave C: promote high-confidence endpoint joins into scaffold/path artifacts, still without FASTA mutation.
4. Wave D: connect multi-*k* selection with repeat-aware simplification policy, including explicit edit/recompress/reannotate/replan behavior.
5. Wave E: escalate audit evidence toward read-backed polishing candidates without default repair.
6. Wave F: deepen diploid ambiguity output and parent-specific graph evidence without full haplotype FASTA claims.

## Immediate Build Bias

Work should happen in large, evidence-gated waves instead of isolated feature drips. Each wave should add the minimum durable data types, unit tests, sidecar outputs, and validator coverage needed for the next promotion stage.

The SPAdes-specific pressure test is captured in
[`docs/SPADES_ARCHITECTURE_INSPIRATION.md`](SPADES_ARCHITECTURE_INSPIRATION.md).
Use it to decide whether a Trex change is borrowing a real assembler
architecture pattern or merely adding an option.
