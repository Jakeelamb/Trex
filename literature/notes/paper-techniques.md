# Paper Technique Breakdowns

This note breaks down the archived PDFs into the main functions, techniques, and concrete Trex implications. It is intentionally implementation-oriented: every paper should either shape a module, a benchmark gate, or a deferred track.

## Lin et al. 2016 — ABruijn

Paper: [`2016-lin-abruijn-long-error-prone-dbg.pdf`](../papers/2016-lin-abruijn-long-error-prone-dbg.pdf)

Main functions:

- Build an A-Bruijn graph by gluing vertices with identical labels across read paths instead of treating the classical short-read DBG as the only graph form.
- Select reliable landmarks / solid sequence features from noisy long reads.
- Compress long noisy reads into graph paths over those landmarks.
- Derive draft paths, split them into local mini-alignments, and correct local segments by consensus.
- Iterate graph construction and correction so noisy reads become usable assembly evidence.

Core techniques:

- DBG/OLC hybridization: keep graph leverage while using long-read overlap-like signal.
- Landmark-driven sparsification: avoid storing every base-level event from noisy reads.
- Local consensus over small graph regions rather than global all-read alignment.
- Repeat sensitivity: treat long reads as repeat-spanning evidence but avoid masking all frequent k-mers blindly.

Trex implications:

- Keep Illumina DBG as the active Phase-2 path, but reserve a future graph adapter for sparse landmark graphs.
- Do not hard-code “DBG means short reads only” into architecture names.
- Build an evidence layer that can score paths from more than one evidence type without immediately mutating the trusted k-mer count table.

## Koren et al. 2017 — Canu

Paper: [`2017-koren-canu.pdf`](../papers/2017-koren-canu.pdf)

Main functions:

- Correct reads.
- Trim reads.
- Overlap reads with adaptive k-mer weighting.
- Build a sparse assembly graph.
- Separate repeats and diverged haplotypes rather than collapsing them too early.
- Emit graph outputs for structures that cannot be represented linearly.

Core techniques:

- Adaptive k-mer weighting: downweight uninformative/repetitive seeds without deleting all repeat signal.
- tf-idf-style MinHash overlap filtering.
- Sparse graph construction.
- Repeat/haplotype separation through graph structure.
- Graph-first output when linear contigs would hide ambiguity.

Trex implications:

- Add copy-number and repeat evidence to Trex graph metadata before adding aggressive simplification.
- Treat GFA as the carrier for unresolved structure, not only as a debug sidecar.
- Future long-read tracks should not reuse the Illumina DBG interface directly; they need an evidence-to-graph adapter.

## Bankevich et al. 2012 — SPAdes

Paper: source-indexed in [`literature/sources.md`](../sources.md). The original
article is not archived locally; the Trex architecture reading is cross-checked
against the SPAdes manual and the archived metaSPAdes paper.

Main functions:

- Correct reads before assembly.
- Build de Bruijn graphs over a ladder of k-mer sizes.
- Use paired-read information as graph-context constraints rather than only as
  final contig ordering.
- Simplify the graph through repeated cleanup of errors, tips, bulges, and
  repeat-like structures.
- Emit sequence outputs together with graph/path artifacts.

Core techniques:

- Multi-*k* graph construction and checkpointable per-*k* iterations.
- Paired de Bruijn graph / k-bimer-style reasoning for distance-aware graph
  constraints.
- Bulge projection/provenance instead of losing all information about collapsed
  alternatives.
- Iterative graph cleanup where topology changes are followed by graph
  normalization before later decisions.
- Assembly graph and path outputs as inspectable deliverables.

Trex implications:

- Keep `--auto-k` and explicit ladders as inspectable select-one evidence until
  cross-*k* graph merging has its own tests and benchmark gates.
- Add read-trust diagnostics before any correction path is allowed to change
  graph construction.
- Grow mate evidence toward oriented graph-context constraints with distance
  histograms and conflict clusters.
- Add a simplification scheduler that explicitly alternates edit,
  recompress, reannotate, and replan.
- Preserve collapsed/retained branch provenance for future repeat and scaffold
  reasoning.

## Nurk et al. 2017 — metaSPAdes

Paper: [`2017-nurk-metaspades.pdf`](../papers/2017-nurk-metaspades.pdf)

Main functions:

- Assemble metagenomic short reads under uneven species coverage.
- Preserve strain variation where isolate-style simplification would collapse it.
- Adapt SPAdes graph ideas to complex communities.
- Resolve graph structures while accounting for microdiversity and nonuniform coverage.

Core techniques:

- Uneven-coverage graph interpretation.
- Strain-aware repeat/bubble handling.
- Conservative graph simplification when coverage cannot be explained by one isolate genome.
- Metagenome-specific benchmark expectations.

Trex implications:

- Do not bake a uniform-coverage assumption into graph simplification interfaces.
- Keep coverage interpretation behind a simplification policy module.
- Add metagenome rows only after Trex has explicit uneven-coverage gates.

## Wick et al. 2017 — Unicycler

Paper: [`2017-wick-unicycler.pdf`](../papers/2017-wick-unicycler.pdf)

Main functions:

- Build an accurate short-read assembly graph.
- Score and select the best k-mer graph.
- Estimate contig multiplicity from depth and graph connectivity.
- Create bridges from short-read repeat-resolution paths and long-read graph alignments.
- Apply bridges in quality order.
- Finalize graph by removing redundant bridged contigs, merging simple paths, rotating circular replicons, and polishing.

Core techniques:

- Multi-k graph selection rather than one fixed k for all cases.
- Multiplicity inference from depth plus graph degree.
- Bridge candidates as scored evidence records before graph mutation.
- Conservative/normal/bold modes that change bridge acceptance thresholds.
- Graph-derived bridge sequence preferred over noisy long-read sequence when possible.

Trex implications:

- Introduce a bridge evidence ledger before allowing mate-derived new edges or scaffolds.
- Keep contigs and scaffolds separate: a bridge is not automatically a contig unless the graph path is representable.
- Add simplification modes only when they are tied to different thresholds and benchmark rows.

## Walker et al. 2014 — Pilon

Paper: [`2014-walker-pilon.pdf`](../papers/2014-walker-pilon.pdf)

Main functions:

- Align reads back to a draft assembly.
- Detect base errors, small indels, gaps, misassemblies, and suspicious collapsed repeats.
- Collect local evidence around suspicious regions.
- Reassemble or correct local regions.
- Emit corrected FASTA plus VCF/changes/tracks as audit artifacts.

Core techniques:

- Evidence-driven local repair: mark suspicious loci, collect anchored reads, then repair.
- Misassembly signals from invalid pairs, soft clips, abnormal coverage, and local coverage drops.
- Collapsed-repeat reporting instead of pretending every suspicious repeat can be safely fixed.
- Rich audit outputs for changes.

Trex implications:

- Add a post-assembly evaluation/repair track later, but keep it separate from primary graph construction.
- For now, add metrics that identify suspicious collapsed/repeated regions rather than silently editing them.
- Every future correction needs a changes artifact, not only a rewritten FASTA.

## Rhie et al. 2020 — Merqury

Paper: [`2020-rhie-merqury.pdf`](../papers/2020-rhie-merqury.pdf)

Main functions:

- Count reliable k-mers from high-accuracy reads.
- Count k-mers in assembly outputs.
- Compare read and assembly k-mer sets.
- Estimate consensus QV and k-mer completeness.
- For trio/parental data, estimate haplotype completeness, phase block continuity, and switch errors.
- Produce spectra plots for copy-number interpretation.

Core techniques:

- Reference-free validation from k-mer set operations.
- Error k-mers: assembly k-mers absent from reliable read k-mers.
- Completeness: reliable read k-mers present in the assembly.
- Copy-number spectra: compare read multiplicity to assembly multiplicity.
- Parent-specific k-mers for haplotype correctness.

Trex implications:

- The next quality framework should be Merqury-like before BUSCO-like: Trex already owns k-mer counting.
- Add read-vs-assembly k-mer metrics to `xtask bench`, not just reference-containment metrics.
- For diploid rows, add parent/haplotype-specific k-mer containment as a first phasing gate.

## Li and Durbin — Genome Assembly In The T2T Era

Paper: [`2023-li-durbin-t2t-era.pdf`](../papers/2023-li-durbin-t2t-era.pdf)

Main functions:

- Classify sequencing evidence by length, accuracy, and long-range information.
- Explain the current near-T2T recipe: accurate long reads, ultra-long reads, and long-range phasing/scaffolding data.
- Compare overlap graphs, string graphs, DBGs, multiplex DBGs, and minimizer/sparse DBGs.
- Summarize polishing and evaluation practices.
- Identify unsolved challenges in repeats, phasing, polyploidy, metagenomes, and structurally variable samples.

Core techniques:

- Match graph IR to evidence: overlap/string graphs for accurate long reads; DBG/multiplex DBG and minimizer-space graphs where they preserve the right information.
- Use sparse/minimizer graph representations to reduce memory.
- Use trio or Hi-C evidence for chromosome-scale phasing.
- Evaluate with assembly size, N50, gene completeness, k-mer quality, alignment-based checks, and manual review for hard genomes.

Trex implications:

- Keep long-read and T2T tracks as explicit future adapters, not hidden changes to Illumina.
- Add minimizer/sparse graph investigation to the performance roadmap.
- Keep evaluation layered: reference-free k-mer quality, reference alignment, gene completeness, and artifact review each answer different questions.

## Cross-Paper Framework

The shared framework is:

1. **Evidence first:** reads, pairs, long reads, parent-specific k-mers, Hi-C, and alignments should become typed evidence before they mutate graph structure.
2. **Graph IR second:** choose the graph representation that preserves the evidence scale: Illumina DBG now; sparse/minimizer or repeat graph later.
3. **Policy-driven simplification:** tips, bubbles, repeat bridges, strain/haplotype branches, and scaffolds need explicit policies, thresholds, and modes.
4. **Structure-preserving output:** unresolved ambiguity belongs in GFA/path/scaffold artifacts rather than being hidden in one FASTA.
5. **Evaluation as a module:** reference-free k-mer quality, reference alignment, phasing, and gene completeness are separate gates.
6. **Auditability:** every graph mutation or post-assembly correction needs evidence and an artifact trail.
