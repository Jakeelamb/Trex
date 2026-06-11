# Assembly Literature Source Ledger

Status key:

- `archived`: verified PDF exists under [`papers/`](papers/).
- `source-indexed`: source/DOI is recorded, but direct public PDF archiving was blocked, non-OA, or requires manual access.
- `selection-needed`: category was named, but no single paper was specified in the request.

## Foundational Concepts

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| source-indexed | Li Z. et al. 2012, OLC vs DBG comparison | <https://academic.oup.com/bfg/article/11/1/25/191455> | - | Terminology and tradeoff baseline for OLC/DBG positioning. |
| archived | Lin Y. et al. 2016, ABruijn long error-prone reads with DBG | <https://www.pnas.org/doi/10.1073/pnas.1604560113>, Europe PMC `PMC5206522` | [`2016-lin-abruijn-long-error-prone-dbg.pdf`](papers/2016-lin-abruijn-long-error-prone-dbg.pdf) | Long-read DBG generalization and OLC/DBG hybrid thinking. |
| source-indexed | Kolmogorov M. et al. 2019, Flye repeat graphs | <https://doi.org/10.1038/s41587-019-0072-8> | - | Repeat graph model for long noisy reads and future repeat-resolution design. |

## Short-Read De Novo Assembly

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| source-indexed | Zerbino and Birney 2008, Velvet | <https://pmc.ncbi.nlm.nih.gov/articles/PMC2336801/>, <https://genome.cshlp.org/content/18/5/821> | - | Classic DBG error cleanup, simplification, and paired-end use. |
| source-indexed | Bankevich A. et al. 2012, SPAdes | <https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/> | - | Multi-k, paired DBG, and single-cell/standard assembly design pressure. |
| source-indexed | Prjibelski A. et al. 2020, SPAdes protocol/update | <https://doi.org/10.1002/cpbi.102> | - | Operator workflow and reproducible SPAdes usage. |
| source-indexed | Li D. et al. 2015, MEGAHIT SdBG | <https://doi.org/10.1093/bioinformatics/btv033> | - | Succinct DBG and large-metagenome single-node scalability. |
| selection-needed | ABySS, SOAPdenovo, ALLPATHS-LG | User requested as “other notables” | - | Add after choosing exact canonical citations. |

## Long-Read Assembly

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| archived | Koren S. et al. 2017, Canu | <https://pmc.ncbi.nlm.nih.gov/articles/PMC5411767/> | [`2017-koren-canu.pdf`](papers/2017-koren-canu.pdf) | OLC repeat separation, adaptive k-mer weighting, and scalable correction/trim/assembly pipeline. |
| source-indexed | Kolmogorov M. et al. 2019, Flye | <https://doi.org/10.1038/s41587-019-0072-8> | - | Same repeat-graph core as above. |
| source-indexed | Ruan J. and Li H. 2020, wtdbg2 / Redbean | <https://doi.org/10.1038/s41592-019-0669-3>, `PMC7004874` | - | Fuzzy Bruijn graph and speed/accuracy tradeoffs. |
| selection-needed | Miniasm/minimap, Shasta, NECAT, NextDenovo, Raven, SmartDenovo | User requested as frequently benchmarked tools | - | Add exact papers as we decide which competitors define Trex benchmarks. |

## HiFi, Diploid, And Phased Assembly

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| source-indexed | Cheng H. et al. 2021, hifiasm | <https://doi.org/10.1038/s41592-020-01056-5>, `PMC7961889` | - | Phased assembly graphs and diploid/haplotype-resolved output contracts. |
| source-indexed | Cheng H. et al. 2025, hifiasm ONT | <https://www.biorxiv.org/content/10.1101/2025.04.14.648685v1> | - | Emerging standard ONT simplex route to near-T2T assemblies. |

## T2T And Near-Complete Genomes

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| source-indexed | Rautiainen M. et al. 2023, Verkko | <https://doi.org/10.1038/s41587-023-01662-6>, `PMC10427740` | - | Iterative graph pipeline for diploid T2T assembly. |
| source-indexed | Verkko2 extensions | User mentioned recent 2025 work | - | Add exact papers once selected. |
| source-indexed | Nurk S. et al. 2022, T2T-CHM13 | <https://www.science.org/doi/10.1126/science.abj6987>, `PMC9186530` | - | Finished-genome validation and artifact standards. |
| source-indexed | Liu J. et al. 2024, T2T mouse | <https://www.science.org/doi/10.1126/science.adq8191> | - | Mammalian benchmark target for near-complete assembly workflows. |
| archived | Li H. and Durbin R., genome assembly in the T2T era | <https://pmc.ncbi.nlm.nih.gov/articles/PMC10462168/>, <https://www.genome.gov/sites/default/files/media/files/2024-10/Genome-assembly-in-the-telomere-to-telomere-era.pdf> | [`2023-li-durbin-t2t-era.pdf`](papers/2023-li-durbin-t2t-era.pdf) | Review-level map for data types, algorithms, evaluation, and T2T workflows. |

## Hybrid Assembly

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| archived | Wick R. et al. 2017, Unicycler | <https://doi.org/10.1371/journal.pcbi.1005595>, <https://pmc.ncbi.nlm.nih.gov/articles/PMC5481147/> | [`2017-wick-unicycler.pdf`](papers/2017-wick-unicycler.pdf) | Hybrid bacterial completion and SPAdes-backed bridge/scaffold behavior. |

## Metagenome Assembly

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| archived | Nurk S. et al. 2017, metaSPAdes | <https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/> | [`2017-nurk-metaspades.pdf`](papers/2017-nurk-metaspades.pdf) | Strain-aware/metagenomic graph handling and uneven coverage pressure. |
| source-indexed | metaFlye / long-read metagenome assembly | <https://doi.org/10.1038/s41592-020-00971-x> | - | Repeat/strain graph behavior under metagenomic variation. |
| source-indexed | MEGAHIT | <https://doi.org/10.1093/bioinformatics/btv033> | - | Same SdBG metagenome scalability target as above. |
| selection-needed | MetaMDBG and related long-read metagenome assemblers | User mentioned as recent resources | - | Add exact papers before archiving. |

## Scaffolding, Phasing, Hi-C, And Finishing

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| selection-needed | 3D-DNA, ALLHiC, HapHiC, SALSA, GPhase, optical mapping, Pore-C | User listed tool families, not specific papers | - | Add exact papers when we scope long-range scaffolding. |

## Polishing, Error Correction, And Evaluation

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| archived | Walker B. et al. 2014, Pilon | <https://doi.org/10.1371/journal.pone.0112963>, <https://pmc.ncbi.nlm.nih.gov/articles/PMC4237348/> | [`2014-walker-pilon.pdf`](papers/2014-walker-pilon.pdf) | Short-read polishing and draft-improvement contracts. |
| selection-needed | Racon, Medaka, Nanopolish | User listed tool families, not specific papers | - | Add exact papers when polishing enters Trex scope. |
| archived | Rhie A. et al. 2020, Merqury | <https://doi.org/10.1186/s13059-020-02134-9> | [`2020-rhie-merqury.pdf`](papers/2020-rhie-merqury.pdf) | K-mer assembly quality, completeness, and phasing assessment. |
| source-indexed | QUAST | <https://doi.org/10.1093/bioinformatics/btt086> | - | Reference assembly evaluation used by Trex benchmark gates. |
| selection-needed | BUSCO, WebQUAST, HERRO, recent benchmark papers | User listed categories and examples | - | Add exact papers as gates are defined. |

## Reviews And Recommendations

| Status | Work | Source | Local artifact | Why Trex cares |
|--------|------|--------|----------------|----------------|
| source-indexed | Yang Y. et al. 2025, recent advances and challenges in de novo genome assembly | <https://www.maxapress.com/article/doi/10.48130/gcomm-0025-0015> | - | Current cross-tool review and roadmap pressure test. |
| source-indexed | Earth BioGenome Project assembly recommendations | <https://www.earthbiogenome.org/report-on-assembly-recommendations> | - | Practical data/pipeline requirements for reference-quality genomes. |
| source-indexed | Kong W. et al. 2023, complex plant genome assembly review | <https://www.sciencedirect.com/science/article/pii/S1672022923000700> | - | Large/repetitive/polyploid genome context. |

## Archived PDF Checksums

```text
005706c0e50cf1b97d771b1666b11589273dbd38c9db8f5b61e32a5b95cacf55  papers/2014-walker-pilon.pdf
8f941223a12d4df3e44fd813851bc12c9047227df2008f14cf969f122621f80a  papers/2016-lin-abruijn-long-error-prone-dbg.pdf
c7758fa994b66a5cc8a8ef58e2dfbdb55711d96a3eaa93ad649e7d5e1493ed53  papers/2017-koren-canu.pdf
b6e82e3f17e2b294488eaf2fbdce9da9c18ea05acee66e52bdece24b0b9966fe  papers/2017-nurk-metaspades.pdf
eb6adb9ea4a79ffde666a0246d01be38f2acea0935e36677cf0979e243d89a03  papers/2017-wick-unicycler.pdf
3b0b49041d72854d77a916408631a67b1ef12a772c9478b049a46e4d923c08f5  papers/2020-rhie-merqury.pdf
b74fbbe6a4cf40e3cdab8932ddc1ee854d878afd42659311d6db90358c6e9724  papers/2023-li-durbin-t2t-era.pdf
```
