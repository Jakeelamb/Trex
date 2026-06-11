# ADR 0001 — Phase-2 deferrals, primary collapse, and GFA scope

## Status

Accepted (documentation anchor for glossary alignment).

## Context

`CONTEXT.md` bounds **Phase-2 Illumina** vs **Phase-2 deferral** (long-read / hybrid / default **GFA 2**).

## Decisions

1. **Deferred (explicit)**  
   - Long-read / hybrid assemblies and **GFA 2** as default interchange remain **out of scope** until the **Phase-2 deferral** block is amended.  
   - **IUPAC** and **`N`** in the **primary** collapsed FASTA for het sites stay **out of scope** unless rescoped (**Phase-2 Illumina primary export**).

2. **Primary FASTA het collapse (in scope)**  
   - Implemented as **overlapping trusted *k*-mer multiplicity** voting per base with **A < C < G < T** tie-break, skipping positions with zero supporting trusted mass (see `trex::illumina::phase2_primary`).

3. **GFA `P` coverage**  
   - **`P`** lines list **full** unitig segment traversals only (including multi-unitig greedy partitions).  
   - Contig vertex paths that are **strict subpaths** of a single unitig still emit **no** `P` row (spec-safe; documented product gap vs “always tie ctg to utg”).

4. **Mate usage**  
   - **Phase-2** mate bridge **only** increments weights on **existing** DBG edges between **R1** last forward *k*-mer and **R2** first forward *k*-mer when an **insert prior** is configured; requires `preprocess/pair_layout.json` on resume for the same boost identity.

## Consequences

- Operators comparing Trex **GFA** to tools that assume richer `P` / CIGAR-trimmed subpaths should treat **`L`** + dual **`p2h`** mirror paths as the current diploid carrier, not full phased walks.
