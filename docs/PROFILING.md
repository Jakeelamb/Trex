# Profiling

This page records measured Trex performance evidence. Keep raw profiler outputs under
`target/profiles/`; commit only commands, summaries, and decisions.

## 2026-06-10 PhiX174 Baseline

Command:

```bash
mkdir -p target/profiles
/usr/bin/time -v target/release/trex illumina assemble \
  --r1 fixtures/phix174/reads.fq \
  --kmer-size 31 \
  --trusted-threshold 1 \
  --out-dir target/profiles/phix174-baseline \
  > target/profiles/phix174.stdout \
  2> target/profiles/phix174.time.txt

cargo flamegraph -p trex-cli --bin trex \
  -o target/profiles/phix174-flamegraph.svg -- \
  illumina assemble \
  --r1 fixtures/phix174/reads.fq \
  --kmer-size 31 \
  --trusted-threshold 1 \
  --out-dir target/profiles/phix174-flamegraph-run
```

Observed baseline:

| Input | Reads | Candidate k-mers | Unique k-mers | Unitigs | Contigs | Wall | Max RSS |
|-------|-------|------------------|---------------|---------|---------|------|---------|
| PhiX174 deterministic 150 bp reads | 108 | 12,960 | 5,386 | 1 | 1 | 27.25 s | 12,668 KiB |
| PhiX174 after unorientable-walk filtering | 108 | 12,960 | 5,386 | 1 | 1 | 26.65 s | 12,596 KiB |
| PhiX174 after `pick_best_neighbor` lookup/set optimization | 108 | 12,960 | 5,386 | 1 | 1 | 17.42 s | 12,408 KiB |

Top `perf report` symbols from `target/profiles/phix174.perf.data`:

| Symbol | Self |
|--------|------|
| `trex::dbg::walk::pick_best_neighbor` | 21.27% |
| `trex::kmer::reverse_complement` | 15.70% |
| `trex::dbg::walk::reference_contig_paths` | 7.24% |
| `BTreeMap<Vec<u8>, ...>::insert` | 4.28% |
| `trex::dbg::unitig::stitch_sequence` | 2.53% |
| `malloc` / `cfree` | 2.97% combined |

Immediate read:

- The current passing micro benchmark is dominated by graph walking and repeated allocation-heavy
  sequence orientation work, not FASTQ parsing.
- `pick_best_neighbor` clones neighbor k-mers and consults `BTreeSet<Vec<u8>>` forbidden sets inside
  greedy walks from every component seed.
- `reverse_complement` returns a fresh `Vec<u8>` and appears both in canonicalization and stitch/orientation paths.
- The next optimization slice should target the DBG walk representation before broad parser work.

## 2026-06-10 Biological Rows

Prepared with:

```bash
cargo run -p xtask -- fetch-data
```

Rows:

| Row | Source | Prepared subset | Current result |
|-----|--------|-----------------|----------------|
| `ecoli_mg1655_srr001666_1k_pairs` | ENA `SRR001666` + RefSeq `GCF_000005845.2` / `NC_000913.3`; 7,047,668 paired spots / 507,432,096 bases | 1,000 R1 + 1,000 R2 reads | Passes: 2,000 reads, 11,880 unique/trusted k-mers, 1,979 contigs, N50 36 bp, 0.16 s, 19,172 KiB RSS. Reference-quality: 10,557 / 11,880 contig k-mers in reference (88.86%), 1,818 / 1,979 contigs with a reference k-mer hit. QUAST basic stats only: all contigs reported unaligned. |
| `ecoli_mg1655_srr001666_10k_pairs` | ENA `SRR001666` + RefSeq `GCF_000005845.2` / `NC_000913.3` | 10,000 R1 + 10,000 R2 reads | Passes: 20,000 reads, 117,829 unique/trusted k-mers, 19,433 contigs, N50 36 bp, 6.61 s, 155,348 KiB RSS. Reference-quality: 97,729 / 117,546 contig k-mers in reference (83.14%), 16,902 / 19,433 contigs with a reference k-mer hit. QUAST basic stats only: all contigs reported unaligned. |
| `yeast_btt_err1308583_diploid_1k_pairs` | ENA `ERR1308583` + RefSeq `GCF_000146045.2_R64`; 14,550,715 paired spots / 2,870,913,582 bases; BTT ploidy table = euploid diploid | 1,000 R1 + 1,000 R2 reads | Passes: 2,000 reads, 149,866 unique/trusted k-mers, 1,884 contigs, N50 101 bp, 15.00 s, 176,680 KiB RSS. Reference-quality: 129,452 / 146,994 contig k-mers in reference (88.07%), 1,853 / 1,884 contigs with a reference k-mer hit. QUAST: genome fraction 2.000%, NA50 101, 1 misassembly, 210 unaligned contigs. |

Immediate read:

- Biological data ingestion and fixture governance are now present.
- The yeast row exposed that candidate walks through a canonical undirected graph may not have one
  valid forward orientation. Trex now filters those candidate paths instead of aborting the whole run.
- The E. coli row exposed that tip clipping was treating isolated sparse linear components as
  removable tips. Trex now only clips a low-coverage tip when it reaches a higher-degree junction.
- Reference FASTAs and SHA-256-pinned larger bacterial data are now governed data, not local notes.
- `xtask bench` now emits fast reference k-mer quality metrics when a Trex row declares a reference,
  and `TREX_RUN_QUAST=1` records optional QUAST / MetaQUAST timing/status per direct Trex row.
- QUAST currently does not align the short E. coli fragments even though reference k-mer containment
  is high; yeast aligns against the public S288C reference but remains a sparse-subset assembly with
  low genome fraction.
- The remaining biological benchmark work is to improve contig continuity and reduce fragment-only
  outputs, then add longer eukaryotic/diploid ladders.
- The next performance blocker on the passing path is `dbg::walk`; profile-guided work should start
  there while keeping biological rows as manual regressions.
- First walk optimization removed redundant edge lookups in `pick_best_neighbor` and replaced ordered
  forbidden sets with hash membership. PhiX wall time dropped from 26.65 s to 17.42 s with unchanged
  assembly metrics.
- SPAdes-style exact overlap merging was tested as a throwaway prototype on the E. coli 10k contigs.
  At `k - 1` overlap it found only 37 joins, kept N50 at 36 bp, and moved max contig length from
  44 bp to 47 bp with unchanged reference k-mer containment. That is not the current production
  bottleneck.
- The next data rungs are governed but manual: E. coli 100k pairs and yeast BTT diploid 10k pairs.
  A pre-optimization E. coli 100k run was stopped after roughly 3 minutes while CPU-bound at about
  1.3 GiB RSS, making component-walk scaling the next performance target rather than post-FASTA
  cleanup.
- Linear/cyclic component walking now uses one deterministic traversal plus reverse-orientation
  scoring instead of launching a greedy walk from every vertex. The E. coli 10k regression row kept
  identical output/reference-quality metrics and improved wall time from 6.61 s to 5.67 s.
