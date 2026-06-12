# SPAdes Architecture Inspiration For Trex

This note translates SPAdes-family design ideas into Trex module work. It is
not a plan to clone SPAdes. The contract remains Trex's evidence-first,
promotion-policy architecture on the single `trex illumina assemble` surface.

## Sources

- Bankevich et al. 2012, SPAdes: source-indexed in
  [`literature/sources.md`](../literature/sources.md).
- Nurk et al. 2017, metaSPAdes:
  [`literature/papers/2017-nurk-metaspades.pdf`](../literature/papers/2017-nurk-metaspades.pdf).
- SPAdes manual:
  <https://ablab.github.io/spades/running.html> and
  <https://ablab.github.io/spades/output.html>.

## What To Borrow

| SPAdes-family idea | Trex interpretation | First Trex module pressure |
|--------------------|---------------------|----------------------------|
| Read correction before graph construction | Separate read trust from graph mutation. Trusted k-mers, correction candidates, rejected spans, and support summaries should be visible before any corrected reads are allowed to influence the graph. | `read_correction_trust` emits report-only correction/trust diagnostics. |
| Multi-*k* assembly as a ladder | Treat *k* values as independent evidence snapshots, not a hidden auto parameter. Start with select-one behavior, then use cross-*k* signals as annotations before any graph merge. | `multi_k_graph_ladder` records completeness, branchiness, dead ends, repeat risk, and selected-*k* checkpoint provenance. |
| Paired de Bruijn / k-bimer thinking | Mate pairs are constraints between graph contexts with orientation and distance confidence, not just endpoint counters. | `mate_evidence` grows from endpoint joins into k-bimer-like path constraints and conflict clusters. |
| Iterative simplification | Graph edits should be followed by recompression, reannotation, and replanning so later decisions see the current topology. | `simplification_policy` gains an edit/recompress/reannotate/replan scheduler. |
| Bulge projection rather than blind deletion | When a branch is collapsed, retain enough provenance to explain where it went and to reuse it for future repeat/path reasoning. | `simplification.json` records projected/retained structures; GFA/path exports can cite them later. |
| Contextual coverage | Coverage interpretation must be local. metaSPAdes shows why one global coverage threshold breaks under uneven depth and related strains. | `repeat_annotation` and `simplification_policy` use local ratios and adjacent-edge context, not only global multiplicity floors. |
| Repeat resolution by path evidence | SPAdes exSPAnder-style behavior combines graph paths, paired reads, and other long-range signals to extend paths only when one choice is sufficiently supported. | `path_scaffold_builder` promotes only representable paths first; FASTA scaffolds stay separate and policy-gated. |
| Graph output as a primary artifact | SPAdes exposes assembly graph, contig paths, scaffold paths, and coverage tags. Trex should make graph/path output inspectable, not a debug afterthought. | `dbg/export.rs`, `scaffolds.json`, `scaffolds.fa`, and tagged GFA paths stay first-class outputs. |

## What Not To Copy Yet

- Do not implement full BayesHammer-style correction before Trex has
  report-only read-trust diagnostics and read-vs-assembly quality gates.
- Do not merge graphs across *k* values until select-one multi-*k* behavior is
  benchmarked on bacterial, diploid, and human-slice rows.
- Do not import metagenomic strain-consensus behavior into the default Illumina
  diploid path. metaSPAdes is useful because it proves coverage is contextual;
  metagenome claims remain a deferred adapter.
- Do not let scaffolds silently replace primary `contigs.fa`. SPAdes recommends
  scaffolds as resulting sequences, but Trex's current contract keeps primary
  contigs conservative and emits scaffold sequence in a separate sidecar.
- Do not add SPAdes-scale repeat surgery without provenance. Trex needs every
  topology edit to cite evidence, policy, and benchmark impact.

## Concrete Build Waves

1. **Read-trust diagnostics**: emit trusted/untrusted k-mer strata, suspicious
   read spans, and correction candidates without changing graph construction.
2. **Richer multi-*k* scoring**: extend `multi_k.json` with read k-mer
   completeness, dead-end rate, branch/tangle counts, repeat-risk counters,
   graph-density pressure, weighted score terms, and selected-*k* checkpoint
   provenance.
3. **K-bimer-like mate constraints**: represent mate evidence as oriented
   constraints between graph contexts with insert-distance histograms, support,
   conflict clusters, and blocker reasons.
4. **Iterative simplification scheduler**: run planned edit, record
   recompress/reannotation hook status, replan the next pass from the current
   topology, and preserve decision ids across passes.
5. **Projected-bulge provenance**: record collapsed or retained bubble
   alternatives so later path/scaffold logic can distinguish deleted noise from
   biologically ambiguous structure.
6. **Path-first repeat resolution**: use mate/path evidence to emit scaffold
   artifacts and GFA paths before any graph edit or primary FASTA mutation.

## Acceptance Rule

SPAdes-inspired work only counts when it lands as one of these durable surfaces:

- a typed Trex module or data type,
- a JSON/GFA/FASTA sidecar with stable schema intent,
- a validator rule,
- a unit or smoke test,
- a benchmark matrix artifact.

Everything else stays research-only.
