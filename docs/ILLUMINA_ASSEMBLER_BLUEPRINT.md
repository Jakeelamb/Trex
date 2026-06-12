# Ultimate Illumina Assembler Blueprint

This is the target architecture contract for making Trex a strongest-possible Illumina assembler without hiding biological uncertainty in one FASTA. It is governed with [`tools/assembler_framework.toml`](../tools/assembler_framework.toml), validated by `cargo run -p xtask -- validate-framework`, and kept in the default `cargo run -p xtask -- validate` gate.

The operating rule is conservative: evidence becomes behavior only through an explicit promotion policy. Report artifacts, GFA paths, scaffolds, graph edits, and polishing edits are separate claims.

## State-Of-The-Art Lessons

| Source family | Lesson Trex adopts | Contract consequence |
|---------------|--------------------|----------------------|
| SPAdes / metaSPAdes | Multi-*k* DBG assembly, contextual coverage, paired-end scaffolding. | Build independent candidate graphs, score them, and avoid treating coverage as globally uniform. |
| Unicycler | Bridges should be evidence-scored and promoted by confidence. | Mate/path bridges first become ranked artifacts, then paths, then sequence only under policy. |
| Pilon | Audit suspicious sequence before correcting it. | Polishing starts as `audit.json` / `audit.tsv`; sequence repair needs explicit before/after evidence. |
| Merqury | K-mer truth is a primary quality signal. | Read-vs-assembly and parent/reference k-mer metrics gate quality and promotion claims. |
| Canu / ABruijn / Flye-style graph thinking | Preserve ambiguity and expose paths instead of forcing one contig too early. | GFA/path artifacts carry unresolved ambiguity; primary FASTA remains conservative. |
| T2T-era and diploid work | Biological interpretation and phasing claims need evidence labels. | Diploid output labels ambiguity and parent-specific support without claiming full haplotype FASTA. |

## Layer 1: Assembler Architecture Blueprint

The target modules below are the top-level seams. Each module must have owned data, tests, and promotion boundaries before it is allowed to change primary output.

### read correction / trusted k-mer model (`read_correction_trust`)

Responsibilities:

- Own FASTQ/FASTA ingestion quality checks, canonical k-mer counts, trusted k-mer thresholds, and future correction candidates.
- Explain why a k-mer, read segment, or corrected base is trusted, rejected, or ambiguous.
- Emit diagnostics before modifying reads.

Core data types:

- `TrustedKmerSet`: canonical k-mer plus count, quality class, and source read support.
- `CorrectionCandidate`: read id/span, original sequence, proposed sequence, support, and rejection reason.
- `TrustSummary`: candidate count, accepted trusted count, rejected count, and low-support strata.

Test obligations:

- Unit tests for canonical count stability, threshold edges, ambiguous bases, and paired input parity.
- Property tests for reverse-complement canonicalization.
- Integration tests proving report-only correction diagnostics do not change `contigs.fa`.

### multi-k graph ladder (`multi_k_graph_ladder`)

Responsibilities:

- Build candidate DBG graphs independently for an explicit k ladder.
- Score candidates using read k-mer completeness, branchiness, contiguity, repeat risk, and benchmark metadata.
- Select one graph first; graph merging is a later contract.

Core data types:

- `KGraphCandidate`: k, graph summary, trusted-read completeness, branch/tangle counters, repeat-risk counters.
- `KSelectionDecision`: selected k, rejected candidates, score terms, and reason.
- `KGraphArtifact`: `multi_k.json` plus selected-k namespace for future checkpoints.

Test obligations:

- Unit tests for score ordering, tie-breaks, empty ladders, duplicate k values, and invalid checkpoint reuse.
- Synthetic tests where low k over-branches and high k fragments.
- E. coli 10k benchmark artifact whenever graph/path behavior changes.

### repeat-aware graph annotation (`repeat_annotation`)

Responsibilities:

- Annotate nodes, edges, unitigs, contig endpoints, and paths with multiplicity and repeat suspicion.
- Separate single-copy, high-copy, low-confidence, and mixed unitig evidence.
- Feed simplification and promotion policies without directly mutating the graph.

Core data types:

- `NodeAnnotation`: multiplicity, depth class, repeat suspicion, and confidence.
- `UnitigAnnotation`: min/mean/max multiplicity, mixed-copy flag, endpoint class, and repeat-risk reason.
- `RepeatRiskSummary`: counts by class and benchmark-exportable totals.

Test obligations:

- Unit tests for single-copy, high-copy, mixed unitig, and low-coverage false-positive cases.
- Smoke tests proving annotation output does not change FASTA/GFA by itself.
- Regression tests for repeat-suspected endpoint join blocking.

### decision-first simplification (`simplification_policy`)

Responsibilities:

- Convert graph plus evidence into planned tip, bubble, repeat, and diploid decisions before mutation.
- Preserve default behavior unless a policy mode is explicitly enabled and benchmarked.
- Refuse edits on repeat-like or diploid-ambiguous structures unless policy allows them.

Core data types:

- `SimplificationDecision`: action, reason, affected vertices, support score, and policy mode.
- `PlannedGraphEdit`: delete, collapse, retain, tag, or defer operation with preconditions.
- `SimplificationSummary`: decisions by action/reason plus mutation-equivalence counters.

Test obligations:

- Unit tests proving decisions match current mutations for default tip and diamond cases.
- Guardrail tests for repeat-like unitigs and retained balanced diploid bubbles.
- Golden smoke tests for unchanged default contigs unless intentionally promoted.

### mate-pair distance/orientation evidence (`mate_evidence`)

Responsibilities:

- Interpret paired-end orientation, insert-distance priors, endpoint support, conflicts, and mate clusters.
- Rank bridge candidates without adding DBG edges by default.
- Provide evidence that can be promoted into scaffold artifacts or paths.

Core data types:

- `MatePairObservation`: read pair, orientation, observed span class, and mapped endpoint/unitig evidence.
- `BridgeCandidate`: source, target, support count, distance estimate, orientation, conflict count, and confidence.
- `MateEvidenceSummary`: existing-edge boosts, absent-edge endpoint joins, rejected/conflicting clusters.

Test obligations:

- Unit tests for FR/RF/FF orientation classes, insert-window acceptance, and conflict rejection.
- Regression test proving existing-edge boost behavior remains unchanged.
- E. coli 10k bridge artifact whenever ranking or promotion changes.

### scaffold/path promotion (`promotion_policy`, `path_scaffold_builder`)

Responsibilities:

- Convert evidence to artifacts through a single explicit promotion ladder.
- Keep unitigs, contigs, GFA paths, scaffold JSON, FASTA scaffolds, graph edits, and polishing edits distinct.
- Require deterministic output and artifact-backed explanations for every promoted claim.

Core data types:

- `PromotionDecision`: candidate id, target stage, evidence ids, confidence, blockers, and required gates.
- `ScaffoldPath`: ordered segments, gap estimates, support, source evidence, and export status.
- `PathExportSummary`: emitted GFA paths, scaffold artifacts, FASTA scaffold candidates, and skipped candidates.

Test obligations:

- Unit tests for promotion threshold edges and conflict precedence.
- Smoke tests proving `contigs.fa` is unchanged when promotion stops before FASTA.
- GFA/path determinism tests for stable ordering and stable evidence ids.

### polishing/audit loop (`polishing_audit`)

Responsibilities:

- Audit low read support, abnormal mate support, repeat-collapse suspicion, and future base-level corrections.
- Emit report artifacts before editing sequence.
- Promote polishing only when k-mer quality, local support, and regression gates agree.

Core data types:

- `AssemblyAuditFinding`: contig/span, class, evidence, severity, and suggested next action.
- `PolishingCandidate`: base/span proposal, read support, k-mer delta, blockers, and expected artifact change.
- `PolishingSummary`: findings by severity and edits accepted/rejected by policy.

Test obligations:

- Unit tests for low-support, mate-disagreement, and repeat-collapse classes.
- Integration tests proving audit-only mode is FASTA-identical.
- Before/after k-mer quality tests before any polishing edit is enabled.

### diploid ambiguity handling (`diploid_ambiguity`)

Responsibilities:

- Label retained heterozygous-like bubbles, parent-specific k-mer support, and ambiguous paths.
- Report balance and support without overclaiming phasing.
- Emit GFA/path annotations that preserve alternatives.

Core data types:

- `DiploidBubbleAnnotation`: branch ids, support balance, parent-specific support, and ambiguity class.
- `ParentKmerEvidence`: parent label, shared/specific k-mer support, and confidence.
- `DiploidSummary`: retained ambiguity, collapsed candidates, parent balance, and claim boundary.

Test obligations:

- Synthetic tests for balanced and imbalanced bubbles.
- Parent-reference tests for parent-specific k-mer labels.
- Yeast diploid benchmark rows for report coverage, not full haplotype claims.

### quality gates and benchmark claims (`assembly_quality`)

Responsibilities:

- Own PR, main, nightly, and manual benchmark contracts.
- Emit JSON artifacts that can substantiate or reject claims.
- Separate correctness, continuity, completeness, phasing, and performance claims.

Core data types:

- `QualityClaim`: claim level, evidence artifact, benchmark row, and allowed wording.
- `BenchArtifactSummary`: row id, tools, metrics, sidecars, hashes, wall time, and max RSS.
- `RegressionGate`: command, expected artifact class, and fail-fast reason.

Test obligations:

- Validator tests for missing artifacts, missing row metadata, and missing claim levels.
- PR smoke for every behavior-changing wave.
- E. coli 100k only for performance/storage claims; yeast diploid only for diploid claims.

## Layer 2: Promotion Policies

Promotion stages are intentionally ordered. A later stage must cite evidence from earlier stages and a gate that is appropriate for the output it changes.

| Stage | Token | Allowed output | Required evidence | Disallowed by default |
|-------|-------|----------------|-------------------|-----------------------|
| Report-only candidate | `report_only_candidate` | JSON/TSV summaries, benchmark metadata, tracing. | Typed evidence id, support count, confidence, and blocker list. | Any FASTA, GFA path, or graph topology change. |
| Scaffold artifact only | `scaffold_artifact` | `scaffolds.json` or equivalent path artifact. | Non-conflicting bridge/path evidence, deterministic ordering, insert/orientation model when mate-derived. | Primary `contigs.fa` mutation and DBG edge insertion. |
| GFA path | `gfa_path` | GFA `P`/path-compatible metadata and tags. | Representable path through existing graph segments, evidence ids, and stable path naming. | Creating sequence not backed by graph/path representation. |
| FASTA scaffold with gaps | `fasta_scaffold_with_gaps` | Separate scaffold FASTA with explicit gaps. | Promotion decision, gap estimate/confidence, path artifact, and benchmark row proving no primary FASTA regression. | Replacing `contigs.fa` silently. |
| Graph edit | `graph_edit` | DBG topology changes such as clipping, collapse, or bridge application. | Planned edit, reason, precondition check, repeat/diploid guardrail result, and golden/regression coverage. | Applying mate-derived absent edges without explicit bridge policy. |
| Polishing edit | `polishing_edit` | Sequence correction artifact or future corrected FASTA. | Audit finding, local read support, k-mer quality improvement, before/after artifact diff, and rollback path. | Base edits from alignment-only evidence or one-row benchmark wins. |

Promotion rules:

1. Default behavior may stop at `report_only_candidate` even when evidence is strong.
2. Any stage that changes FASTA or graph topology requires docs, changelog, unit tests, and a benchmark artifact in the same wave.
3. Repeat-suspected, parent-specific, or conflicting evidence lowers promotion unless the policy names an explicit override.
4. A performance claim cannot be made from tests; it needs stored benchmark artifacts.
5. A biological claim must name whether it is reference-free, reference-backed, parent-backed, diploid-aware, or deferred.

## Layer 3: Wave-Based Implementation

Wave A: commit/checkpoint current sidecar/evidence stack.

- Output goal: clean reviewable diff for evidence, annotation, simplification, scaffold, multi-*k*, audit, diploid, fragmentation, endpoint join, and GFA metadata work.
- Verification: `cargo fmt --all --check`, `cargo clippy --workspace --all-features -- -D warnings`, `cargo test --workspace --all-features`, `cargo run -p xtask -- validate`, `cargo run -p xtask -- gate --tier pr`.

Wave B: mate orientation + distance model.

- Output goal: `BridgeCandidate` records include orientation class, insert-distance estimate, support/conflict counts, and confidence.
- Verification: mate orientation unit tests, endpoint conflict tests, unchanged primary FASTA smoke, E. coli 10k bridge JSON artifact.

Wave C: endpoint-join promotion into scaffold paths, still no FASTA mutation.

- Output goal: high-confidence endpoint joins promote from report-only candidates to deterministic scaffold/path artifacts.
- Verification: scaffold artifact determinism tests, GFA/path metadata tests where representable, PR gate, E. coli 10k benchmark.

Wave D: multi-k plus repeat-aware simplification policy.

- Output goal: multi-*k* selected graph and repeat-aware simplification guardrails interact through policy summaries rather than hidden mutation.
- Verification: synthetic ladder fixtures, repeat-veto unit tests, default one-*k* golden unchanged, E. coli 10k graph/path benchmark.

Wave E: read-backed polishing/audit escalation.

- Output goal: audit findings gain enough read/k-mer support detail to classify future polishing candidates while remaining report-only by default.
- Verification: audit unit tests, read-vs-assembly k-mer quality deltas, no FASTA repair in default mode, yeast diploid audit row coverage.

Wave F: diploid-aware ambiguity output.

- Output goal: retained heterozygous-like alternatives and parent-specific evidence are labeled in JSON/GFA without claiming full haplotype FASTA.
- Verification: Phase-2 synthetic parent evidence, yeast diploid metrics, GFA tag tests, claim-boundary docs.

## Data Flow

1. Ingest reads and optional references.
2. Build trusted k-mer evidence and read-quality diagnostics.
3. Build one or more candidate DBG graphs.
4. Annotate graph nodes, unitigs, repeats, endpoints, and ambiguity.
5. Plan simplification decisions and apply only policy-approved graph edits.
6. Interpret mate orientation/distance support and emit bridge candidates.
7. Promote evidence into scaffold artifacts, GFA paths, FASTA scaffolds, graph edits, or polishing candidates only through policy.
8. Export primary unitigs, primary contigs, GFA, scaffold sidecars, audit reports, diploid evidence, and benchmark JSON.
9. Use quality gates to decide whether a wave can claim observed, tested, benchmarked, or production-gated status.

## Unit Test Map

| Module | Minimum unit tests | Integration or benchmark tests |
|--------|--------------------|--------------------------------|
| `read_correction_trust` | canonical counts, trust thresholds, ambiguous-base handling, correction candidate rejection. | Report-only correction leaves smoke FASTA unchanged. |
| `multi_k_graph_ladder` | score terms, tie-breaks, invalid ladder input, checkpoint rejection. | Synthetic low/high-*k* ladder and E. coli 10k artifact. |
| `repeat_annotation` | single-copy, high-copy, mixed unitig, endpoint repeat class. | Repeat-suspected endpoint join remains unpromoted. |
| `simplification_policy` | decision/action equivalence, tip/bubble reasons, repeat/diploid vetoes. | Golden smoke unchanged unless policy changes are intentional. |
| `mate_evidence` | orientation, insert window, conflict clusters, bridge confidence. | E. coli 10k bridge artifact and unchanged primary FASTA. |
| `promotion_policy` | stage ordering, conflict precedence, missing-evidence rejection. | Scaffold/GFA artifact determinism and no accidental FASTA mutation. |
| `polishing_audit` | low support, mate disagreement, repeat collapse, candidate blockers. | Default audit-only FASTA identity and yeast diploid report coverage. |
| `diploid_ambiguity` | balanced bubble, parent-specific k-mer labels, ambiguous path class. | Phase-2 synthetic parent metrics and yeast diploid claim boundary. |
| `assembly_quality` | validator failure modes, JSON summary shapes, claim-level parsing. | PR gate, E. coli 10k, E. coli 100k, yeast diploid rows by wave. |

## Claim Boundaries

- Primary `contigs.fa` is conservative and must not silently absorb scaffolding or polishing behavior.
- `scaffolds.json`, `fragmentation.json`, `audit.json`, `diploid.json`, and benchmark JSON are evidence artifacts, not sequence claims.
- GFA path annotations can expose alternatives and evidence before FASTA changes.
- Full haplotype FASTA, hybrid assembly, long-read assembly, and metagenomic interpretation remain deferred until their adapters and claims are explicitly promoted.
