# Trex Literature Review

This directory is the project reading queue for assembler architecture, benchmark design, and product positioning.

## Layout

| Path | Purpose |
|------|---------|
| [`papers/`](papers/) | Public PDFs successfully archived from publisher, PLOS, Europe PMC, BMC, or genome.gov endpoints. |
| [`sources.md`](sources.md) | Source ledger for every requested paper or paper family, including DOI/source links and archive status. |
| [`notes/`](notes/) | Review notes that should become Trex design decisions, benchmarks, or ADRs. |

## Reading Order

1. **DBG baseline:** Velvet, SPAdes, MEGAHIT, and metaSPAdes. Extract graph simplification, paired-read handling, and memory representation ideas.
2. **Long-read graph evolution:** ABruijn, Flye, wtdbg2, and Canu. Separate ideas that require long reads from ideas that improve graph representation generally.
3. **Diploid and T2T frontier:** hifiasm, Verkko, T2T-CHM13, and Li/Durbin. Convert these into Trex phase gates for haplotype graphs, repeat resolution, and validation artifacts.
4. **Hybrid and finishing:** Unicycler, Pilon, Merqury, QUAST/BUSCO. Use these to define what Trex measures before claiming assembly quality.
5. **Reviews and recommendations:** Yang 2025, EBP recommendations, plant-genome review. Use these for benchmark matrix expansion and product-roadmap pressure tests.

## Trex Extraction Targets

- **Graph representation:** succinct or compact DBG storage, repeat graphs, fuzzy/approximate graph nodes, and path compression.
- **Paired and long-range evidence:** when to add new edges, when to only score existing edges, and when to output scaffolds separately from contigs.
- **Diploid semantics:** phased graph outputs, bubble preservation/collapse rules, and validation against haplotype-aware metrics.
- **Quality gates:** reference-free k-mer quality, QUAST/BUSCO/Merqury-style checks, and artifact retention for bisectable regressions.
- **Performance posture:** memory ceilings, single-node scalability, and deterministic bounded benchmark ladders.

Do not turn review notes into product claims until there is a runnable Trex gate for the claim.
