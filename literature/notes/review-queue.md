# Review Queue

Use this file to turn papers into concrete Trex work. Each review should end in one of:

- benchmark row
- ADR
- Rust implementation issue
- rejected idea with measured reason

## Immediate Reviews

1. **Merqury:** define a Trex reference-free k-mer quality row beyond simple reference containment.
2. **Canu + ABruijn:** identify which long-read ideas translate to compact graph storage without changing Trex's Illumina contract.
3. **Unicycler + SPAdes:** decide when Trex should emit scaffolds separately from contigs, especially for paired reads.
4. **Li/Durbin T2T era:** convert recommended data/evaluation layers into future benchmark matrix tiers.
5. **metaSPAdes:** extract uneven-coverage and strain/bubble handling ideas before adding metagenome rows.

## Missing Specific Citations

The user requested several tool families without exact papers: ABySS, SOAPdenovo, ALLPATHS-LG, Miniasm/minimap, Shasta, NECAT, NextDenovo, Raven, SmartDenovo, MetaMDBG, 3D-DNA, ALLHiC, HapHiC, SALSA, GPhase, Racon, Medaka, Nanopolish, BUSCO, WebQUAST, HERRO, and recent 2024-2025 benchmarks.

Before archiving those, choose the canonical citation for each family and add it to [`../sources.md`](../sources.md).
