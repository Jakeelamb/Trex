# Trex

Rust genome assembler: correctness and resource use are first-class; biological scope expands in deliberate phases. **Active engineering** is now the **Phase-2 Illumina endgame**: Illumina-only eukaryotic diploid assembly on the single **`trex illumina assemble`** path, while long-read and hybrid work remain governed by **Phase-2 deferral**. Further Phase-2 detail ships via grilling, **`ARCHITECTURE.md`**, or **ADRs**, not unbounded glossary growth.

## Language

**Phase-1 target**:
Small, effectively haploid genomes (e.g. viral or bacterial) from a single sample; the baseline for correctness claims, benchmarks, and core data structures before scaling complexity.
_Avoid_: MVP, toy genome (unless explicitly synthetic controls)

**Phase-2 target**:
Eukaryotic diploid assembly, including human-relevant repeat load and heterozygosity; phasing and diploid-specific logic belong here, not in Phase-1.
_Avoid_: “Later”, generic scaling (be specific: diploid eukaryote)

**Phase-2 Illumina endgame acceptance**:
Trex is accepted as a governed **Illumina diploid assembler** when **Phase-1** compatibility remains green, **Phase-2 Illumina** has a tiered benchmark matrix from tiny synthetic cases through real eukaryotic or human-slice datasets, outputs include primary FASTA plus inspectable haplotype graph structure, and reproducibility, checkpointing, CI cadence, tool pins, fixture provenance, and operator docs are stable enough that a new governed matrix row can be added without redesigning the pipeline.
_Avoid_: Treating one synthetic diploid fixture as completion, accepting aggregate contig stats without haplotype graph inspection, adding benchmark rows through ad hoc scripts or undocumented downloads

**Phase-2 Illumina primary export**:
The operator-facing **primary** contig sequence stream is a **single** FASTA representing a documented **collapse** at unresolved heterozygous sites; **haplotype-resolved** or **dual-haplotype** structure is expressed primarily in **GFA 1.0** (segments, links, paths, and Trex-documented tags), not as a mandatory pair of full-length haplotype FASTA files in the first Illumina diploid milestone. **Collapse** in that FASTA chooses **exactly one** **A/C/G/T** per column by **deterministic** rules grounded in **k-mer multiplicity** or equivalent count-derived signals, in the same spirit as **Phase-1 bubble resolution**; **IUPAC ambiguity** and **`N`** in the primary stream for SNP-level heterozygous collapse are **out of scope** unless explicitly rescoped.
When parent references are supplied, Trex may emit parent-specific k-mer evidence in sidecars and GFA tags. That evidence is a support annotation, not a claim that emitted FASTA is fully phased or haplotype-resolved.
_Avoid_: Treating two independent haplotype FASTAs as the default Phase-2 v1 contract, hiding diploid structure only in sidecar formats not covered by Trex GFA documentation, stochastic or undocumented tie-breaks at het collapse, default IUPAC or **N** in the primary FASTA for simple hets

**Phase-2 mate usage (Illumina)**:
Paired-end reads **may** inform **graph simplification** and **diploid-aware bubble handling** on the **Phase-1 contigging IR** after the **trusted de Bruijn graph** for the run is fixed from **Phase-1 k-mer count representation**; they **do not** change which canonical *k*-mers are counted, merged, or admitted as trusted for graph construction. When mate-derived signals depend on **numeric insert distance** or an insert-size **distribution**, **Phase-1 insert prior policy** applies in full, the same as for scaffold **`U`** records; **adjacency-only** hints that do not invent numeric gap lengths remain valid **without** that prior only where documented for Trex. This **relaxes** **Phase-1 mate usage** only for **Phase-2 Illumina** work, not for **Phase-1 target** claims.
_Avoid_: PE-aware *k*-mer enumeration or trust rules that rewrite the counted multiset, using mate hints to silently change *k* or *T*, claiming **Phase-1 mate usage** unchanged while wiring PE into simplification without a **Phase-2** anchor, silent insert models for distance-sensitive simplification while **Phase-1 insert prior policy** still governs scaffolding

**Phase-2 Illumina benchmark gate**:
**Phase-2 Illumina** acceptance runs the **Phase-1 benchmark gate** unchanged as a **first layer**, then runs **separate** diploid-specific checks documented for Trex; operator and CI reporting must be able to attribute failure to **Phase-1** vs **Phase-2 Illumina** layers without conflating vocabulary. The diploid-specific **reference-aligned** layer requires at least one fixture with **two parental haplotypes** or equivalent **diploid ground truth** so metrics can target **haplotype** correctness, not only **collapse** against a single haploid reference. The **default CI** diploid fixture is **synthetic** (reads and truth fully pinned in-repo for reproducibility); **biological** read sets are **out of scope** for that default gate unless explicitly rescoped. On the non–reference-aligned side, **Phase-1 reference-free metrics** apply **only** to the **primary** FASTA stream; any **GFA-derived** checklist items are **Phase-2 Illumina graph summaries**, not **Phase-1 reference-free metrics**.
_Avoid_: A single undifferentiated gate script that mixes haploid and diploid metrics without layered failure semantics, dropping or weakening **Phase-1 benchmark gate** requirements when enabling diploid CI, labeling haplotype-derived scores as **Phase-1 reference-free** without documentation, claiming **Phase-2 Illumina** diploid benchmark completeness with only a single haploid truth sequence, pinning nondeterministic or externally hosted biological data in the default CI gate without documentation, labeling **GFA-only** counts as **Phase-1 reference-free metrics**

**Phase-2 Illumina benchmark tiers**:
Benchmark progress is named by evidence class: **Tier 0** is in-repo synthetic PR smoke, **Tier 1** is governed small synthetic or real-reference nightly coverage, **Tier 2** is governed real eukaryotic diploid evidence that may remain manual or release-candidate when external data is too heavy for CI, and **Tier 3** is governed human-slice or GIAB-style evidence. Passing a lower tier does not imply claims from a higher tier.
_Avoid_: Calling synthetic diploid smoke a biological diploid gate, treating a manual biological row as default CI, starting Tier 3 without a governed human-slice row and explicit artifact contract

**Phase-2 Illumina graph summaries**:
Counts or summaries computed from emitted **GFA 1.0** (and Trex-documented tags) **without** reference sequence or alignments; they may appear in **Phase-2 Illumina benchmark gate** checklists alongside **Phase-1 reference-free metrics** but **must not** reuse the **Phase-1 reference-free metrics** name or definition, which stays restricted to **N50**, **total assembled sequence length**, **contig count**, and documented siblings on the **primary** FASTA stream.
_Avoid_: Treating graph topology statistics as interchangeable with **Phase-1 reference-free metrics**, omitting **Phase-1 reference-free metrics** from the Phase-2 Illumina layer, smuggling alignment-derived quantities into **Phase-2 Illumina graph summaries**

**Phase-2 Illumina CLI surface**:
**Phase-2 Illumina** behaviour is selected through the same **`trex illumina assemble`** surface as **Phase-1**, using **explicit flags** and/or **documented configuration**; **Phase-1 target** semantics remain the **default** when those selectors are off. A **separate subcommand** is **not required** for the first Illumina diploid milestone. With selectors **off**, emitted **FASTA**, **GFA**, and related streams contain **no diploid-only** records, tags, or paths; output shape stays **Phase-1-shaped** so **Phase-1 target** claims and tooling parity hold unless **Phase-1** itself is deliberately replanned.
_Avoid_: Silent diploid behaviour with no explicit operator selector, implying a second binary for diploid assembly, duplicating the entire flag matrix under a parallel subcommand without documentation, emitting **Phase-2 Illumina**-only **GFA** tags or paths in default mode, undocumented optional fields that force parsers to branch on build flavor

**Phase-2 Illumina versioning**:
Public **semver** treats **Phase-2 Illumina** features as **additive** (minor or patch bumps) unless a change deliberately breaks **Phase-1 target** contracts; **`CHANGELOG.md`** and operator-facing docs label **Phase-2 Illumina** selectors and outputs **experimental** until the **Phase-2 Illumina benchmark gate** is declared **stable** in project documentation. Prerelease tags such as **`-beta`** are **not** required solely because diploid code ships behind explicit selectors.
_Avoid_: Major semver solely for off-by-default diploid surfaces, removing experimental labelling while diploid gates are still in flux, conflating marketing “Phase-2” with **semver** bumps

**Phase-2 Illumina counting and trust**:
Through construction of the **trusted de Bruijn graph**, **Phase-1 k-mer count representation**, **Phase-1 trusted k-mer rule**, **Phase-1 local threshold policy**, **Phase-1 k-mer identity**, **Phase-1 N policy**, and related **Phase-1** counting vocabulary apply **unchanged** in **Phase-2 Illumina** mode unless **Phase-1** is deliberately replanned; **Phase-2 mate usage (Illumina)** applies only **after** that counting boundary.
_Avoid_: PE-informed *k*-mer trust, spatially varying *T*, duplicate forward and reverse keys in counting without a glossary replan

**Phase-2 Illumina simplification**:
**Phase-2 Illumina** extends **Phase-1 graph simplification** with **documented diploid-aware** motifs (including bubble classes) on the fixed **trusted** graph, respecting **Phase-1 bubble bounds** in spirit with **documented** diploid-specific caps at least as strict unless the operator widens them; **Phase-1 simplified graph invariants**, including **loop-free** self-edge prohibition after automatic simplification, remain in force unless **Phase-1** invariants are explicitly revised.
_Avoid_: Unbounded diploid bubble surgery, silent widening of topology or length budgets, tolerating self-loops as normal output

**Phase-2 Illumina checkpointing**:
Checkpoint trees **may** include **Phase-2 Illumina**-only artifacts **only** when **Phase-2 Illumina** selectors are **on**; with selectors **off**, checkpoint payloads remain **Phase-1-shaped**. **Phase-1 checkpoint integrity** and **strict** resume semantics apply to shared stages without silent cross-mode resume.
_Avoid_: Diploid-only blobs in default **Phase-1** checkpoints, ambiguous resume across mismatched mode manifests

**Phase-2 Illumina observability**:
With **Phase-2 Illumina** selectors **on**, **`tracing`** must expose the **active diploid profile** and whether an **explicit insert-size prior** was supplied wherever **Phase-2 mate usage (Illumina)** touches **distance-sensitive** rules; with selectors **off**, **Phase-1 observability** shape is unchanged.
_Avoid_: Mode-blind logs, reconstructing insert priors from telemetry instead of operator configuration

**Phase-2 Illumina phasing ladder**:
Early milestones keep **haplotype-resolved** structure in **GFA 1.0** per **Phase-2 Illumina primary export**; **mandatory additional haplotype FASTA** products beyond the **primary** collapsed stream wait for an explicit glossary and **Phase-2 Illumina benchmark gate** amendment.
_Avoid_: v1 promises of full-length phased pair FASTAs, diploid claims without inspectable **GFA** structure

**Phase-2 Illumina reproducibility**:
**Phase-1** multiset, merged-count, and **unitig** determinism policies apply unchanged through counting; **Phase-2 Illumina** **GFA** and **primary FASTA** collapse semantics document their **stability tier** in **`CHANGELOG.md`** until explicitly frozen; **Phase-1 randomness policy** extends to diploid simplification choices.
_Avoid_: RNG in diploid tie-breaks, silent promises of byte-identical diploid **GFA** across versions before freeze

**Phase-2 Illumina preprocess**:
Illumina feeds in **Phase-2 Illumina** mode inherit **Phase-1 preprocess scope**, **Phase-1 N policy**, **Phase-1 quality usage**, **Phase-1 default qual filters**, **Phase-1 IUPAC ambiguity policy**, **Phase-1 case normalization**, **Phase-1 empty read policy**, **Phase-1 paired read layout**, and **Phase-1 pair parity** unless **Phase-1** preprocess is explicitly replanned.
_Avoid_: Long-read preprocess rules in the Illumina-only slice, silent loosening of **Phase-1 pair parity** for diploid

**Phase-2 Illumina Tier-1 platform**:
**Phase-1 Tier-1 platform** remains the **only** mandatory CI triple for **Phase-2 Illumina** until this glossary is amended.
_Avoid_: Declaring additional Tier-1 triples without dedicated runners and policy updates

**Phase-2 Illumina contig construction**:
The **primary** contig stream obeys **Phase-1 contig construction**, **Phase-1 contig walk rule**, **Phase-1 disconnected graph policy**, **Phase-1 contig walk score**, and **Phase-1 contig tie-break** unless revised with documentation; additional documented walkers or exports for diploid structure must not contradict those rules for the **primary** stream without a glossary amendment.
_Avoid_: Undocumented vertex revisits on **primary** contigs, random walk tie-breaks

**Phase-2 Illumina API stability**:
The **`trex`** crate remains **0.x** with **Phase-1 API stability** expectations during **Phase-2 Illumina** evolution; public Rust surfaces for diploid selectors may change between releases until explicitly stabilized.
_Avoid_: Stable semver promises for diploid APIs while the **GFA** contract is still moving

**Phase-2 Illumina scaffold inheritance**:
Whenever **Phase-2 Illumina** emits scaffold or gap structure, **Phase-1 scaffold encoding**, **Phase-1 GFA gap policy**, and **Phase-1 insert prior policy** apply unchanged.
_Avoid_: Numeric **`U`** gaps without prior, diploid-specific gap records outside **GFA 1.0** without documentation

**Phase-2 Illumina FASTA and GFA surface**:
With selectors **on**, **FASTA** headers obey **Phase-1 FASTA header policy** for **utg**/**ctg** namespaces; any diploid-specific segment or path naming uses **documented disjoint** prefixes or auxiliary tags that **cannot** collide with those reserved namespaces. **GFA** interchange stays **GFA 1.0** per **Phase-1 GFA export**; diploid annotations use a **documented** tag and record subset, not **GFA 2** as default.
_Avoid_: Header collisions, undocumented magic tags, default **GFA 2** interchange

**Phase-2 Illumina benchmark automation**:
Automation implementing the **Phase-2 Illumina benchmark gate** is **separate** from **`scripts/benchmark_gate.sh`**, **invokes** the **Phase-1 benchmark gate** first, then runs the diploid layer so failures stay attributable. External programs for the diploid reference-aligned leg follow **Phase-1 tool pin manifest** and **Phase-1 fixture integrity** expectations.
_Avoid_: Monolithic scripts that blur layer boundaries, unpinned aligners or simulators on governed rows

**Phase-2 Illumina benchmark fixtures**:
Allow-listed **Phase-2 Illumina** matrix rows follow **Phase-1 benchmark data policy** and **Phase-1 fixture integrity**, with digests and provenance for **reads**, **both parental haplotypes**, and any simulators; **default CI** rows remain **synthetic** per **Phase-2 Illumina benchmark gate**.
_Avoid_: Undocumented biological defaults in CI, missing parental truth digests, conflating matrix rows with **Phase-1 benchmark fixtures** without separate labelling

**Phase-2 Illumina error and CLI inheritance**:
**Phase-1 error typing**, **Phase-1 CLI exit policy**, **Phase-1 CLI configuration** precedence, **Phase-1 CLI I/O model**, **Phase-1 CPU/async boundary**, and **Phase-1 default log verbosity** apply to **`trex-cli`** and **`trex`** during **Phase-2 Illumina** work unless narrowly rescoped; diploid-specific fatal outcomes remain **typed** and **documented**, not severity-driven exit hacks.
_Avoid_: `anyhow::Error` as the primary **`trex`** error surface, tying exit status to **`ERROR`** logs without a documented rule

**Phase-2 Illumina fuzzing**:
**Phase-1 fuzzing policy** applies; optional **Phase-2 Illumina** fuzz targets (for example **GFA** round-trips) are **developer or scheduled** tools, not mandatory **per-pull-request** merge gates, unless this entry is amended.
_Avoid_: PR-blocking diploid fuzz without runner budget

**Phase-2 Illumina scope and CI cadence**:
A **Phase-2 Illumina** run consumes **Illumina short reads only**; **same-run hybrid** with **long-read** data waits for the **long-read** sub-track described in **Phase-2 deferral**. **Phase-1 PR smoke** stays the default **pull-request** burden; the **full Phase-2 Illumina benchmark gate** runs on **`main`**, **scheduled nightly**, and **release tags** (or an equivalent documented schedule), not on every pull request by default, until this entry is amended. **`main`** must remain green for the **Phase-1 benchmark gate**; once **Phase-2 Illumina benchmark automation** exists, that gate also runs on **`main`** without relaxing **Phase-1** requirements.
_Avoid_: HiFi or ONT in the first Illumina diploid slice, mandatory full diploid matrix on every PR without runner justification, omitting the diploid gate on **`main`** once it exists

**Phase-2 deferral**:
When **Phase-2** reopens after Illumina **Phase-1** is stable: **Phase-2 target** remains the biological contract. **Engineering may begin with Illumina short-read extensions** (diploid-aware behaviour on the **Phase-1 contigging IR** and exports) before any long-read sub-track ships. When the **long-read** sub-track opens: **PacBio HiFi** ingest and error modeling first, then an **overlap/string-graph** HiFi contigging IR (not the Illumina DBG as default for HiFi); **CLR** waits until HiFi is stable; **ONT** follows HiFi; **diploid-aware long-read graph** work waits until a **haploid HiFi** path is gated; operators still use a **single `trex` binary** with modes. The **Phase-2** glossary is bounded by **Phase-2 target**, the **Phase-2 Illumina** entries from **Phase-2 Illumina primary export** through **Phase-2 Illumina scope and CI cadence**, and this **Phase-2 deferral** entry until focus shifts—finer decisions return via grilling or ADRs then.
_Avoid_: Sprawling ad hoc Phase-2 entries outside that bounded block without grilling or ADRs, treating this sketch as a day-to-day implementation spec, mandating long-read before Illumina diploid work when Phase-2 opens

**Current Phase-2 Illumina implementation boundary**:
The current code path is an **experimental overlay** on **`trex illumina assemble`**: **Phase-1** preprocess, canonical counting, trusted *k*-mer filtering, and DBG construction run first and are unchanged; only after the trusted DBG exists may **Phase-2 Illumina** code alter edge weights, simplification behavior, primary contig collapse, checkpoint identity, observability, and **GFA 1.0** annotations. Today that overlay is limited to **existing-edge** mate boosts gated by **`--diploid`**, paired input, and an explicit insert mean; conflict-aware endpoint join promotion in sidecar scaffold/path artifacts, including separate **`scaffolds.fa`** and tagged **`scf...`** GFA **`P`** rows when accepted paths exist, without replacing or mutating primary **`contigs.fa`**; repeat-aware and near-balanced diamond retention; deterministic trusted-*k*-mer primary FASTA collapse; **`XX:Z:trex-phase2-illumina`**, optional **`L`** rows, primary **`ctg...`** **`P`** rows where full-unitig partitioning is possible, and **`p2h...`** unphased mirror **`P`** rows. This entry documents the shipped boundary; product vocabulary remains the surrounding **Phase-2 Illumina** block and implementation details live in **`ARCHITECTURE.md`** and ADRs.
_Avoid_: Treating current **`p2h...`** mirror paths as fully phased haplotypes, assuming mate boosts can create new graph edges, moving diploid logic before the trusted graph without a replan, reading this boundary as permission for long-read or hybrid work

**Active engineering focus**:
Scheduling and implementation attention is promoted to the **Phase-2 Illumina endgame**: stable Illumina-only eukaryotic diploid assembly, richer GFA/haplotype structure, governed benchmark matrices, and production-grade infrastructure on the existing **`trex illumina assemble`** surface. **Phase-1 target** remains a non-relaxing compatibility layer, not the ceiling; long-read, hybrid, HiFi, CLR, and ONT work remain deferred.
_Avoid_: Treating Phase-1 smoke success as the final assembler target, reopening long-read or hybrid work without amending **Phase-2 deferral**, adding a second Illumina assembler surface for the endgame

**Phase-1 read technology**:
Illumina-style short reads only (single-end and/or paired-end as supported); long-read platforms and hybrid pipelines are explicitly out of scope for the Phase-1 milestone.
_Avoid_: HiFi, ONT, PacBio, hybrid (for Phase-1 claims)

**Phase-1 contigging IR**:
The de Bruijn graph on *k*-mers; Phase-1 contigging is defined as well-defined transformations on this graph until contig sequences are emitted.
_Avoid_: Overlap graph, string graph, dual IR (for Phase-1 scope)

**Phase-1 k policy**:
By default, exactly one user-chosen *k* per assembly run; one de Bruijn graph is built at that *k* without an internal multi-*k* ladder or automatic *k* switching. Explicit `--kmer-ladder` / `k_ladder` and `--auto-k` modes may build and score candidate graphs, but they must be opt-in, emit `multi_k.json`, select only one graph for the normal assembly path, and isolate selected-*k* checkpoints so counts, graph, and export artifacts from one chosen *k* cannot be reused for another.
_Avoid_: Hidden *k* selection, graph merging across *k* values, or treating explicit multi-*k* / auto-*k* mode as the Phase-1 default

**Phase-1 k feasibility**:
If *k* exceeds the **shortest post-preprocess read length** among all reads feeding the run, **`trex-cli`** fails with a **hard error** before *k*-mer enumeration rather than clamping *k* or completing a meaningless count pass in Phase-1.
_Avoid_: Silent *k* clamping, warning-only behaviour that still emits empty or misleading assemblies (for Phase-1)

**Phase-1 mate usage**:
Paired-end reads affect linking, orientation, or scaffolding only after unitigs or contigs exist; mates do not change de Bruijn graph construction or core graph simplification in Phase-1.
_Avoid_: PE-aware DBG, paired-end graph traversal (for Phase-1)

**Phase-1 insert prior policy**:
Distance-sensitive scaffold joins or **`U`** gap records require an **explicit insert-size prior** supplied by the operator (documented parameterization); without a prior, mate-based steps may only assert **relative orientation or adjacency hints** that do not invent numeric gap lengths in Phase-1.
_Avoid_: Silent default insert distributions, mandatory priors even for orientation-only workflows, undocumented automatic prior fitting (for Phase-1)

**Phase-1 scaffold encoding**:
Scaffold or mate-derived order, orientation, and gap structure is expressed using **GFA 1.0** path and gap conventions only, as documented for Trex; Phase-1 does not require a separate mandatory scaffold interchange file (e.g. AGP) alongside FASTA/GFA.
_Avoid_: Scaffold information only in undocumented TSV, requiring AGP as a parallel Phase-1 contract (unless explicitly rescoped)

**Phase-1 GFA gap policy**:
Unknown-distance gaps use **GFA 1.0 `U` lines** with **explicit integer gap-size fields** per spec and Trex docs; **`N`** padding inside `S` sequences appears only where documentation ties it to those **`U`** records, not as a silent replacement for gap metadata in Phase-1.
_Avoid_: Long `N` runs without accompanying `U` records, non-integer or undocumented gap sizes in emitted GFA (for Phase-1)

**Phase-1 quality usage**:
Phred scores may affect trimming or read dropping only in preprocess stages before graph construction; the de Bruijn layer consumes sequence plus a fixed **N** policy, not per-base qualities.
_Avoid_: Quality-weighted k-mers, probabilistic DBG (for Phase-1)

**Phase-1 default qual filters**:
Phred-based **end trimming** and **aggregate whole-read quality dropping** are both **off by default** and must be explicitly enabled in preprocessing configuration; any enabled filter is declared in benchmark metadata.
_Avoid_: Default trimming, default whole-read drops, undocumented implicit quality filtering (for Phase-1)

**Phase-1 N policy**:
No *k*-mer spans an **N**; reads split into **N**-free segments and *k*-mers are emitted only inside segments (no base invention across ambiguity).
_Avoid_: Mapping **N** to a fixed concrete base for counting, whole-read discard as the default for a single **N**

**Phase-1 benchmark gate**:
Phase-1 acceptance requires both reference-aligned metrics against a known truth genome and reference-free summary assembly statistics; one class alone is not enough.
_Avoid_: Reference-only gates, summary-only gates (for Phase-1)

**Phase-1 reference-free metrics**:
The reference-free side of **Phase-1 benchmark gate** is defined on **N50**, **total assembled sequence length**, and **contig count** (plus any explicitly documented sibling summaries), without reference-derived quantities or gene-completeness scores in Phase-1.
_Avoid_: Labeling NG50 or alignment-based stats as reference-free, requiring BUSCO in the Phase-1 reference-free gate (unless deliberately rescoped)

**Phase-1 k-mer count representation**:
*k*-mers are enumerated, sorted, and equal neighbors merged into counts; this sorted run-length form is the canonical frequency table backing the Phase-1 pipeline.
_Avoid_: Hash-map-first as the long-term design center, sketch-only counts (for Phase-1)

**Phase-1 k-mer sort key**:
Enumeration and sorting use **only** the canonical *k*-mer string as the ordering key; read identifiers, coordinates, or other metadata are not mixed into the sort or merge that builds multiplicity tables in Phase-1.
_Avoid_: Stable sort tie-breakers on read name inside the count pipeline, asserting a canonical permutation of duplicate identical *k*-mers before merge (for Phase-1)

**Phase-1 parallel sort policy**:
Parallel *k*-mer sorts **need not be stable** with respect to duplicate identical *k*-mer strings; implementations must preserve the **multiset of *k*-mers** and resulting **merged counts**, not a canonical permutation of duplicates before merging.
_Avoid_: Tests that assert a specific ordering among equal keys, assuming single-threaded sort side effects (for Phase-1)

**Phase-1 default k-mer sort parallelism**:
*k*-mer sorting for **Phase-1 k-mer count representation** may use **parallel algorithms by default**; correctness is defined by **Phase-1 parallel sort policy** and merged counts, not by bitwise-identical intermediate sort streams across threads or runs.
_Avoid_: Mandating single-threaded default sorts, claiming reproducible parallel intermediate dumps (for Phase-1)

**Phase-1 k-mer identity**:
Each counted *k*-mer is identified by the lesser of its forward sequence and reverse complement under a fixed total order; graph vertices use this same canonical identity.
_Avoid_: Treating forward and reverse complement as unrelated keys, unspecified ordering (for Phase-1)

**Phase-1 counting orientation**:
Frequency tables count **canonical *k*-mers only**, matching **Phase-1 k-mer identity** end-to-end, so strand collapse happens **before** graph construction rather than by merging duplicate forward and reverse keys later in Phase-1.
_Avoid_: Separate forward and reverse tallies that must be reconciled at graph build, forward-only counts against canonicalized graph nodes (for Phase-1)

**Phase-1 canonical alphabet**:
Canonical *k*-mer comparison uses the total order **A < C < G < T** only; counted *k*-mers never contain **N**, consistent with **Phase-1 N policy**, so **N** does not participate in lexicographic orientation decisions.
_Avoid_: IUPAC symbols inside counted *k*-mers without an explicit separate policy, silently including **N** in *k*-mers while claiming **Phase-1 N policy** (for Phase-1)

**Phase-1 graph simplification**:
Contigging removes dead-end tips and resolves only small, topology- and length-bounded bubble motifs; arbitrary complex repeat bouquets are out of scope for Phase-1.
_Avoid_: Unbounded bubble merging, SPAdes-scale repeat surgery (for Phase-1)

**Phase-1 simplified graph invariants**:
After **Phase-1 graph simplification**, the working graph is **loop-free** in the sense that **self-adjacent edges on a single vertex are forbidden**; their presence is treated as a **hard error** or another documented **fatal abort** of contigging, not a silent repair, in Phase-1.
_Avoid_: Exporting self-loop `L` records as normal, silently stripping loops without failing the run (for Phase-1)

**Phase-1 bubble resolution**:
Within bounded bubbles, simplification uses **coverage from existing k-mer multiplicity information** (per branch or edge weights derived from sorted counts) to choose retained sequence; ties are broken by a fixed deterministic rule, not by randomness.
_Avoid_: Lexicographic-only bubble resolution ignoring counts, stochastic bubble collapse (for Phase-1)

**Phase-1 bubble bounds**:
Automatic bubble simplification applies **both** a **maximum internal sequence-length budget** and a **maximum graph node or edge budget**, each with documented defaults and operator overrides; motifs exceeding either bound are left unresolved by automatic rules in Phase-1.
_Avoid_: Relying on only length or only topology caps, silently enlarging bounds without documentation (for Phase-1)

**Phase-1 tip clipping**:
Dead-end **tips** are removed only when they are shorter than a configured length bound **and** carry aggregated multiplicity below a configured floor, both derived from **Phase-1 k-mer count representation** signals; **documented global defaults** apply unless operators override them via CLI or **Phase-1 CLI configuration**.
_Avoid_: Length-only clipping ignoring coverage, multiplicity-only clipping with no length guard, requiring operators to invent thresholds with no defaults (for Phase-1)

**Phase-1 trusted k-mer rule**:
A canonical *k*-mer is used in graph construction iff its total count is at least a single global threshold *T* chosen per run; no per-position trust model in Phase-1.
_Avoid_: Threshold-free inclusion of all observed *k*-mers, per-read Bayesian trust (for Phase-1)

**Phase-1 local threshold policy**:
The count floor *T* in **Phase-1 trusted k-mer rule** is **spatially uniform** for a run; Phase-1 does **not** apply **position- or region-dependent** thresholds or silent local rescues that change *T* inside subgraphs.
_Avoid_: Automatic local *T* relaxation in tips or bubbles, undocumented geography-sensitive counting rules (for Phase-1)

**Phase-1 sequence artifacts**:
Phase-1 emits **unitigs** as the graph-faithful sequences after simplification and separately emits **contigs** produced by an explicit walk or compaction policy; structural tests anchor on unitigs while release gates may score contigs.
_Avoid_: Contigs without unitigs as the structural contract, collapsing the two artifacts into one unnamed product (for Phase-1)

**Phase-1 FASTA header policy**:
FASTA headers for **unitigs** and **contigs** use **disjoint identifier namespaces** (distinct reserved prefixes such as `utg` and `ctg`) so IDs cannot collide between artifact kinds in a combined export or log stream.
_Avoid_: Shared prefixes across artifacts, ambiguous numeric-only headers with no type tag (for Phase-1)

**Phase-1 FASTA sequence payload**:
FASTA records for **unitigs** and **contigs** carry the **full assembled sequence strings** for those artifacts; headers are not stand-ins that omit sequence in favor of hashes or external blobs in Phase-1.
_Avoid_: Sequence-less FASTA stubs for primary outputs, any mismatch between FASTA bodies and matching **GFA 1.0 `S`** sequences (for Phase-1)

**Phase-1 read ingest**:
Inputs are FASTQ with optional gzip and FASTA records; FASTA has no per-base qualities and still passes through preprocess and **Phase-1 N policy**; paired-end layout follows **Phase-1 paired read layout**.
_Avoid_: SAM/BAM as a Phase-1 primary read source, interleaved paired FASTQ (for Phase-1)

**Phase-1 empty read policy**:
Any read whose sequence length is **zero** after preprocess is a **hard error** in Phase-1; empty reads are not silently dropped or treated as ignorable padding.
_Avoid_: Silent truncation to empty, continuing assembly with zero-length records (for Phase-1)

**Phase-1 case normalization**:
ASCII nucleotide letters in reads are **uppercased deterministically** during preprocess so *k*-mer machinery and **Phase-1 canonical alphabet** logic see a uniform case; lowercase alone is not a hard error in Phase-1.
_Avoid_: Case-sensitive DNA comparisons in the graph layer, rejecting otherwise valid reads solely for lowercase letters (for Phase-1)

**Phase-1 IUPAC ambiguity policy**:
Non-**ACGT** IUPAC codes (other than **`N`**) are **mapped to `N`** using a **documented translation table** during preprocess so downstream counting sees at most **A/C/G/T/N** before **Phase-1 N policy** applies; ambiguity symbols are not an automatic hard error in Phase-1.
_Avoid_: Silent ad hoc mappings, expanding the counted alphabet beyond **Phase-1 canonical alphabet** plus **N** handling (for Phase-1)

**Phase-1 paired read layout**:
Paired-end Illumina is specified only as **two explicit read files** per documented CLI flags or positional conventions; **interleaved** paired FASTQ and other single-stream pair encodings are not supported in Phase-1.
_Avoid_: Interleaved-first UX, silently treating one file as interleaved pairs (for Phase-1)

**Phase-1 pair parity**:
After preprocessing, paired **R1** and **R2** inputs must contain the **same number of reads**, and reads at each index must share the **same pair identifier** parsed from headers per documented rules; **count mismatches** or **name mismatches** are **hard errors** that abort assembly rather than silently truncating, padding, or reordering in Phase-1.
_Avoid_: Pairing to the shorter file without error, padding with synthetic reads, ignoring header identity while checking only counts (for Phase-1)

**Phase-1 contig construction**:
Contigs come from pluggable walk strategies over the simplified graph; Phase-1 ships exactly one reference walker that **Phase-1 benchmark gate** exercises, without claiming that walker is the final algorithm.
_Avoid_: Only unitigs forever, a single monolithic contigger with no strategy boundary (for Phase-1 architecture)

**Phase-1 contig walk rule**:
The shipped reference walker emits walks that visit each graph vertex **at most once**; vertex revisits and cyclic traversals are not valid Phase-1 contigs under this rule.
_Avoid_: Allowing arbitrary vertex repeats in contigs without documenting a different mode, promising full Eulerian tours as contigs (for Phase-1)

**Phase-1 disconnected graph policy**:
If the simplified graph has **multiple connected components**, the reference walker produces **independent contigs per component**, each obeying **Phase-1 contig walk score** and **Phase-1 contig tie-break** within that component, rather than discarding secondary components or failing solely for disconnectedness in Phase-1.
_Avoid_: Emitting only the highest-scoring component, hard errors for benign disconnectivity during early simplification (for Phase-1)

**Phase-1 contig walk score**:
Candidate **vertex-simple** walks are ranked primarily by **total supporting multiplicity** along their arcs, derived from **Phase-1 k-mer count representation**; **Phase-1 contig tie-break** lex rules apply only when this aggregate score ties.
_Avoid_: Ranking walks by length alone before coverage, using lexicographic ordering before any declared numeric score (for Phase-1)

**Phase-1 contig tie-break**:
When several **vertex-simple** walks tie on **Phase-1 contig walk score**, the reference walker keeps the walk whose contig sequence has the **lexicographically smallest prefix** under **Phase-1 canonical alphabet** ordering at the first differing position; deeper parity follows the same deterministic lex extension rule in documentation.
_Avoid_: Random tie breaks, emitting every tied maximal walk as its own contig by default (for Phase-1)

**Phase-1 unsafe policy**:
`unsafe` is allowed only in workspace crates whose names begin with **`trex-sys-`** or **`trex-simd-`**, each with required safety documentation and review expectations; all other crates remain free of `unsafe` unless a narrowly documented exemption exists.
_Avoid_: `unsafe` scattered in unprefixed crates, project-wide `#![forbid(unsafe_code)]` that blocks planned SIMD work without documented exceptions, unconstrained `unsafe` outside the prefixed crates (for Phase-1)

**Phase-1 Tier-1 platform**:
CI and release binaries treat **x86_64-unknown-linux-gnu** as the only mandatory target; other triples are not Phase-1 gates.
_Avoid_: Declaring aarch64 or macOS as Tier-1 without dedicated runners (for Phase-1)

**Phase-1 reference metrics toolchain**:
The reference-aligned portion of **Phase-1 benchmark gate** is computed using **pinned external executables** (aligner plus a QUAST-class or equivalent reporter), invoked from bench or CI with explicit version pins.
_Avoid_: Unpinned “whatever minimap2 is on PATH”, requiring in-Rust aligners for Phase-1 gates (unless deliberately rescoped)

**Phase-1 tool pin manifest**:
Pinned external tool **names, versions, and install or verification metadata** for **Phase-1 reference metrics toolchain** and for any **read simulator or other generator** whose outputs feed **Phase-1 benchmark matrix** entries live in **one machine-readable manifest** consumed by scripts and CI, rather than scattered prose-only notes.
_Avoid_: README-only pinning, divergent duplicate pin tables per dataset without documented exceptions, unpinned simulators that change read bytes on matrix rows (for Phase-1)

**Phase-1 export formats**:
The supported CLI emits FASTA for unitigs and contigs and emits GFA for the simplified graph or unitig encoding per documented dialect and version.
_Avoid_: FASTA-only exports, unspecified GFA version (for Phase-1)

**Phase-1 export layout**:
Output **paths and grouping** of FASTA and GFA are **operator-configurable**, but the **default** is **separate files per artifact** (distinct default filenames for unitigs, contigs, and graph GFA) rather than a single merged FASTA unless explicitly overridden per documented CLI rules. **Standard output** is used **only** when an artifact’s output path is **explicitly `-` (or an equivalent documented sentinel)**, not by omitting paths in Phase-1.
_Avoid_: A single combined FASTA as the only option, undocumented default filenames, implicit stdout when no `-o` is given (for Phase-1)

**Phase-1 GFA export**:
Normative graph interchange is **GFA 1.0**, using a documented subset of segment, link, and path records; **GFA 2** is out of scope for Phase-1 normative exports.
_Avoid_: Defaulting on GFA 2, mixing multiple GFA versions without explicit flags (for Phase-1)

**Phase-1 GFA segment naming**:
`S` line segment identifiers match the **FASTA headers** for the same **unitigs** and **contigs** per **Phase-1 FASTA header policy**, so GFA and FASTA share one namespace for those artifacts without an extra mapping table in Phase-1.
_Avoid_: Internal-only GFA segment names that diverge from FASTA, silent renaming between outputs (for Phase-1)

**Phase-1 GFA link emission**:
`L` link lines are written **only when required** to express documented **path**, **scaffold**, or other relationships that are not already captured by **`P`** lines or segment sequences alone; purely redundant `L` records for edges implied elsewhere are omitted in Phase-1.
_Avoid_: Dumping a fully redundant `L` mesh, emitting scaffolds without the `L`/`P` structure they need (for Phase-1)

**Trex CLI topology**:
End users invoke a **single `trex` binary** with explicit **subcommands or modes** for distinct read regimes (Illumina versus HiFi and later additions), instead of a family of separately named binaries per regime.
_Avoid_: One-off binaries per modality as the default UX, hiding modes inside undocumented flags (for Phase-2+ product shape)

**Phase-1 stage checkpointing**:
Each **major assembly phase** exposes **operator-visible** **checkpoint export** and **resume from checkpoint** through **`trex-cli`** so runs can be **interrupted and continued** at stage boundaries and each phase can be **benchmarked or swapped** in isolation without always restarting from raw reads in Phase-1.
_Avoid_: Stage boundaries that exist only for tests or the library with **no** user-visible checkpoint path (for Phase-1)

**Phase-1 checkpoint integrity**:
Checkpoint artifacts **may** carry **manifest** metadata (format version and **SHA-256** or documented equivalent of payloads); **default resume** tolerates skipping full digest verification for fast local iteration, while a **documented strict mode** (e.g. **`--strict-checkpoints`**) requires manifest match before continue—on par with **Phase-1 fixture integrity** when operators opt in.
_Avoid_: No strict path for cross-machine or CI resume, strict mode that ignores mismatch or truncation (for Phase-1)

**Phase-1 checkpoint strictness in benchmark gate**:
**Phase-1 benchmark gate** and the **full-matrix** leg of **Phase-1 CI matrix cadence** are **not** required to enable **strict checkpoint verification** by default when they exercise checkpoint resume; strict mode stays an **explicit** flag or job configuration unless this glossary is amended.
_Avoid_: Mandating `--strict-checkpoints` (or equivalent) on every gate or scheduled full-matrix run that touches resume (for Phase-1)

**Phase-1 checkpoint layout anchor**:
Exported **stage checkpoints** for a run are grouped under a **documented per-run root** (path or equivalent operator-visible root) with **stable, documented identifiers** for **artifact kinds** by stage and role; exact byte layouts, on-disk tree shape beyond that anchor, and manifest file naming are specified in **`ARCHITECTURE.md`** once the merge gate in **Phase-1 architecture documentation** applies, or in **ADRs** until then—not as fine-grained glossary entries in Phase-1.
_Avoid_: Undocumented scatter of checkpoint files, encoding binary layout minutiae in `CONTEXT.md` before formats exist (for Phase-1)

**Phase promotion policy**:
**Phase-2 target** work may proceed on parallel branches without an automatic switch that demotes **Phase-1 target** gating on `main`; changing what `main` represents is an explicit decision, not a calendar default. When **`main`** adopts the **Phase-2 Illumina benchmark gate**, **Phase-1 benchmark gate** remains a **non-relaxing** subset per **Phase-2 Illumina scope and CI cadence**.
_Avoid_: Calendar-driven Phase-2 takeover, silently weakening Phase-1 gates when side branches merge unrelated work (for governance), turning off **Phase-1 benchmark gate** on **`main`** solely because diploid checks exist

**Phase-1 CI matrix cadence**:
Pull requests run a **fast smoke** subset of **Phase-1 benchmark gate** checks; the **full benchmark matrix** runs on a fixed schedule such as nightly and on **release tags**.
_Avoid_: Full matrix on every PR without runner budget, only-on-tag full matrix with no scheduled baseline (for Phase-1 drift)

**Phase-1 PR smoke scope**:
Per-pull-request smoke includes **Rust automated tests** and **Phase-1 reference-free metrics** evaluated on a **small named set** of **tiny in-repository synthetic** fixtures (typically a handful of deliberately distinct cases), not the full **Phase-1 benchmark matrix**; it does **not** invoke the external **Phase-1 reference metrics toolchain** by default on every pull request.
_Avoid_: Demanding full reference-based gates on each PR, smoke that omits reference-free checks entirely, relying on a single undifferentiated toy genome as the only smoke input (for Phase-1)

**Phase-1 workspace layout**:
The Cargo workspace centers on a **`trex` library crate** and a separate **`trex-cli`** binary crate from the first merged end-to-end pipeline onward; any **`trex-sys-*`** or **`trex-simd-*`** crates mandated by **Phase-1 unsafe policy** live in the **same repository workspace**, not as separate submodules or unpublished side repos in Phase-1.
_Avoid_: Indefinite single-crate layout, mandatory third benchmark-only crate in Phase-1 (unless later rescoped), hiding `unsafe` crates outside the main workspace without documentation (for Phase-1)

**Phase-1 MSRV policy**:
The workspace declares a **minimum supported Rust version** in Cargo metadata, and CI runs at least one job on that exact toolchain, not only on latest stable.
_Avoid_: Implicit “whatever stable is today”, distro-only MSRV claims without CI enforcement (for Phase-1)

**Phase-1 MSRV extensions**:
`trex-sys-*` and assembly-heavy **`trex-simd-*`** integration ships behind **default-off Cargo features** so plain `trex` / `trex-cli` builds keep the **workspace MSRV** from **Phase-1 MSRV policy**; opting into those features may advertise a **higher rustc floor** in the relevant crate metadata without silently changing the base MSRV in Phase-1.
_Avoid_: Newer compiler requirements on default builds, undocumented nightly-only hooks (for Phase-1)

**Phase-1 CLI I/O model**:
The default `trex-cli` build is **async-first**, using an async runtime for streaming reads and writes rather than a purely synchronous I/O stack in Phase-1.
_Avoid_: Sync-only Phase-1 CLI, async relegated to non-default features while still advertising streaming ingestion (for Phase-1)

**Phase-1 CLI configuration**:
`trex-cli` supports an **optional TOML configuration file** at a documented default path or an explicit path flag; its fields override defaults when present, **command-line flags override the file**, and the file overrides built-in defaults when fields are unset on the CLI.
_Avoid_: Requiring a config file for basic invocations, supporting multiple ambiguous structured formats without precedence rules, unspecified conflict rules between flags and files (for Phase-1)

**Phase-1 CLI async runtime**:
The default `trex-cli` build uses **Tokio** with the **multi-thread** scheduler as the canonical async runtime in Phase-1.
_Avoid_: `async-std` or `smol` as the default runtime, current-thread Tokio as the default without an explicit rescoping decision (for Phase-1)

**Phase-1 CPU/async boundary**:
The **`trex` library** is **synchronous** at its public API; **`trex-cli`** uses Tokio for async streaming I/O and invokes blocking calls into the sync core rather than async-coloring the library in Phase-1.
_Avoid_: `async fn` graph-wide in `trex`, doing heavy CPU work directly on Tokio async worker threads without an explicit off-load boundary (for Phase-1)

**Phase-1 observability**:
Both **`trex`** and **`trex-cli`** use **`tracing`** spans and events for assembly stages; ad hoc stderr logging is not the primary operator interface in Phase-1.
_Avoid_: `println!`-driven UX as the norm, tracing only in the CLI while the library stays silent (for Phase-1)

**Phase-1 default log verbosity**:
At default settings, **`tracing`** emits **INFO**-level events for **major assembly stage boundaries** only; finer-grained diagnostics require explicit higher verbosity flags.
_Avoid_: Default **DEBUG** noise, default **WARN**-only logs that hide stage transitions (for Phase-1)

**Phase-1 CLI exit policy**:
`trex-cli` uses **non-zero exit codes** only for **documented fatal outcomes** represented by **`trex` error returns** (or their wrapped equivalents); **`ERROR`-level `tracing` events** do not by themselves change exit status unless explicitly specified for a narrow diagnostic mode in Phase-1.
_Avoid_: Tying exit code to log severity heuristics, zero exit on returned hard errors (for Phase-1)

**Phase-1 reproducibility**:
Identical inputs, parameters, and pinned tools reproduce **deterministic k-mer frequency semantics and unitig outputs**; **contig** sequences from **Phase-1 contig construction** may change across Trex versions until walker tie-break rules are explicitly stabilized.
_Avoid_: Promising byte-identical contigs in Phase-1, disclaiming all structural reproducibility and relying only on aggregate benchmarks (for Phase-1)

**Phase-1 randomness policy**:
The **`trex`** library uses **no stochastic algorithms** in Phase-1: no random number generators, no randomized default hashers, and tie-breaks are pure functions of inputs and parameters.
_Avoid_: Hidden RNG in graph simplification, statistical algorithms without explicit deterministic modes (for Phase-1)

**Phase-1 preprocess scope**:
Adapter and barcode removal, beyond quality-based trimming or filtering, exists only as an **explicitly opt-in** preprocess mode; default runs rely on upstream cleanliness except where **Phase-1 quality usage** and **Phase-1 N policy** apply.
_Avoid_: Always-on silent adapter trimming, pretending adapters never matter for benchmarks (for Phase-1)

**Phase-1 fuzzing policy**:
FASTQ, FASTA, and gzip ingest paths include **documented fuzz targets** intended for developer or scheduled execution; fuzzing is **not** a mandatory per-pull-request merge gate in Phase-1. **Deterministic fuzz seeds** live in-repo under a documented path such as **`fuzz/corpus/`**, kept **small** and version-controlled rather than gitignored or downloaded ad hoc in Phase-1.
_Avoid_: Skipping fuzz harnesses entirely, per-PR fuzz gates without runner support, huge binary-only corpora committed without justification (for Phase-1)

**Phase-1 API stability**:
The **`trex` crate** stays on **0.x** with an explicitly **unstable** public Rust API during Phase-1; releases may break library APIs without a major bump story, and consumers pin exact versions or revisions.
_Avoid_: Strict semver promises while the graph IR is still churning, undocumented “anything goes” breakage (for Phase-1)

**Phase-1 error typing**:
The **`trex` library** surfaces failures through **`thiserror`-based enums scoped to modules**, with **`#[non_exhaustive]`** on public error enums at crate boundaries as documented; `trex-cli` may add human-oriented reporting without collapsing everything to untyped strings in Phase-1.
_Avoid_: `anyhow::Error` as the primary library error type, a single undifferentiated crate-wide error enum for all subsystems, exhaustively matched public error enums that block evolution (for Phase-1)

**Phase-1 changelog policy**:
User-visible **release notes** for **`trex`** and **`trex-cli`** accumulate in **`CHANGELOG.md` (or a documented equivalent) at tagged releases**, summarizing API and behaviour changes since the prior tag, not necessarily on every intermediate merge in Phase-1.
_Avoid_: Absent changelog entries on tags, requiring changelog edits for every pull request regardless of user impact (for Phase-1)

**Phase-1 architecture documentation**:
The **merge gate** requiring a root **`ARCHITECTURE.md`** plus **at least one versioned diagram** (ASCII or Mermaid) of the Illumina de Bruijn **Phase-1 contigging IR** pipeline starts with the **first `main` pull request that implements graph / DBG** (or equivalent contigging-IR) code—**not** checkpoint-only changes, **not** empty-repo bootstrap (optional stub until then). Once required, the document tracks workspace layout, dataflow, and IR boundaries and the diagram stays aligned with stages; **`rustdoc`** augments but does not replace it.
_Avoid_: Forcing full architecture before DBG work exists, treating checkpoint I/O alone as the DBG trigger, prose-only overview with no pipeline diagram after the gate applies (for Phase-1)

**Phase-1 architecture diagram scope**:
The merge-gated **versioned diagram** must cover the **de Bruijn / contigging** pipeline; **checkpoint export and resume** may live in a **separate** versioned diagram or major section with **cross-references** to the contigging map instead of being crammed into one figure in Phase-1.
_Avoid_: Mandating one combined diagram when split views are clearer, leaving checkpoint I/O undocumented in `ARCHITECTURE.md` once formats exist (for Phase-1)

**Trex license**:
Original Trex sources are released under **`MIT OR Apache-2.0`** at the recipient’s choice unless an individual file carries a different SPDX notice.
_Avoid_: Ambiguous “open source” without SPDX, mixing third-party or benchmark bundle licensing into this statement without separate attribution files (for the project as a whole)

**Phase-1 benchmark data policy**:
Datasets used in the automated **Phase-1 benchmark matrix** appear on a **maintained allow-list** that records **license**, **provenance URL**, **digest** expectations, and a **required coverage or depth class label** (for example an expected Illumina depth band) used when reporting matrix results; ad hoc URLs are not part of the governed matrix.
_Avoid_: Benchmarking against arbitrary download links in CI without license review, conflating dataset licenses with **Trex license**, omitting coverage metadata on allow-listed rows (for Phase-1)

**Phase-1 benchmark fixtures**:
The repository ships **only compact synthetic fixtures** by default; full reference genomes and read sets for the **Phase-1 benchmark matrix** are fetched or built via **pinned scripts** or CI steps, not committed wholesale to git.
_Avoid_: Giant default clones, ad hoc “put E. coli here” instructions without automation (for Phase-1)

**Phase-1 fixture integrity**:
Download and preparation scripts for **Phase-1 benchmark fixtures** must verify each fetched artifact against a **recorded SHA-256 digest** (or stronger documented hash) before those files participate in **Phase-1 benchmark gate** runs.
_Avoid_: Trust-on-HTTPS without local digest verification, verifying only uncompressed sizes (for Phase-1)

**Phase-1 default trusted threshold**:
The CLI’s default global *k*-mer count floor for **Phase-1 trusted k-mer rule** is **T = 2** unless the user overrides it; documentation calls out lower values as non-default expert behaviour.
_Avoid_: No default threshold, undocumented default of **T = 1** (for Phase-1)
