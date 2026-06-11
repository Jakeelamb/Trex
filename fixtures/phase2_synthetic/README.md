# Phase-2 Illumina synthetic benchmark fixture

Two **32 bp** parental haplotypes (`parent1.fa`, `parent2.fa`) differing at one SNP (`p2` has **`T`** at column 10, 1-based). **`reads.fq`** holds four short Illumina-style reads used by **`scripts/phase2_illumina_diploid_reference_layer.sh`**: each read must occur as a substring of **`p1`** or **`p2`** (reference-aligned sanity without invoking Trex diploid assembly).

**SHA-256** digests for these three files are recorded under **`[fixtures.phase2_synthetic]`** in **`tools/manifest.toml`**; the diploid reference layer verifies them before other checks.

**`scripts/phase2_illumina_graph_summaries.sh`** runs **`trex illumina assemble --diploid`** on `reads.fq`, prints **Phase-1-style** summary statistics on **`contigs.fa`** only, and prints **Phase-2 Illumina graph summaries** (GFA record counts and the `trex-phase2-illumina` header tag).

This satisfies **Phase-2 Illumina benchmark gate** requirements for **synthetic** two-parent truth in the **default CI** diploid layer until **`trex illumina assemble`** gains explicit **Phase-2 Illumina** selectors.
