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
| `yeast_btt_err1308583_diploid_1k_pairs` | ENA `ERR1308583` + RefSeq `GCF_000146045.2_R64`; 14,550,715 paired spots / 2,870,913,582 bases; BTT ploidy table = euploid diploid | 1,000 R1 + 1,000 R2 reads | Passes: 2,000 reads, 149,866 unique/trusted k-mers, 1,885 contigs, N50 101 bp, 9.10 s, 214,292 KiB RSS. Reference-quality: 129,490 / 147,049 contig k-mers in reference (88.06%), 1,854 / 1,885 contigs with a reference k-mer hit. QUAST remains optional for this manual row. |

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
- `pick_best_stitchable_path` now computes candidate edge score before stitching and skips candidates
  that cannot beat the best already-stitched score. On `ecoli_mg1655_srr001666_10k_pairs`, wall time
  moved from 6.20 s (`target/benchmarks/ecoli-10k-before-walk-prune.json`) to 5.40 s
  (`target/benchmarks/ecoli-10k-after-walk-prune-final.json`) with identical contig/unitig counts,
  reference k-mer fraction, read containment, and assembly-only k-mer count.
- `stitch_from` now preallocates the stitched sequence and picks between direct/reverse-complement
  orientations without allocating and sorting a temporary option vector at every step. On the same
  row, repeated measurement landed at 5.37 s and 154,768 KiB RSS
  (`target/benchmarks/ecoli-10k-after-stitch-opts-2.json`) with unchanged assembly metrics.
- `pick_best_stitchable_path` now accepts any candidate iterator, so branching components can stream
  greedy paths directly instead of first collecting every candidate path. The same E. coli 10k row
  measured 4.98 s and 153,440 KiB RSS (`target/benchmarks/ecoli-10k-after-candidate-stream.json`),
  again with identical contig/unitig counts and k-mer quality metrics.
- Wave 6 introduced a dense internal `NodeId` view for `DbgGraph` and moved reference contig walking
  onto that view for connected components, linear-component walks, greedy neighbor selection, edge
  weights, and diploid multiplicity tie-breaks. This is a representation step only; do not attach a
  performance claim until a fresh biological before/after artifact is captured.
- Post-change E. coli 10k completed through the governed manual row at
  `target/benchmarks/ecoli-10k-compact-id.json`: 20,000 reads, 117,829 unique/trusted k-mers,
  19,433 unitigs/contigs, N50 36 bp, 6.19 s GNU time, 206,716 KiB max RSS, 99.74% reliable-read
  k-mer containment, and 0 assembly-only k-mers. The E. coli 100k sentinel was attempted at
  `target/benchmarks/ecoli-100k-compact-id.json` and stopped after 606.58 s / 1,499,588 KiB RSS
  before export completion (exit 143), so 100k remains the next graph-storage proof target.
- Keeping the compact walk in `NodeId` form until the selected path is known improved the E. coli 10k
  compact row to 5.73 s at `target/benchmarks/ecoli-10k-compact-id-id-stitch.json` with unchanged
  contig/unitig counts, reference quality, read containment, and 0 assembly-only k-mers. Max RSS was
  206,904 KiB, so compact graph memory reduction is still not proven by this slice.
- A fresh E. coli 100k attempt on the current compact-walk code wrote
  `target/benchmarks/ecoli-100k-current-compact-walk.json` and was stopped after 900.25 s /
  1,499,620 KiB RSS (exit 143). It reached graph simplification logging
  (`tips_removed=13784`) but emitted no assembly metrics, so the remaining 100k bottleneck is
  downstream of simplification and before completed exports/metrics.
- Stage timing showed `primary_contig_paths_for_gfa` was scanning all unitigs for each contig offset:
  on E. coli 10k it took 5.08 s inside `target/benchmarks/ecoli-10k-stage-timing-probe.json`.
  Indexing unitigs by first/last vertex reduced that stage to 51 ms and the full row to 3.29 s in
  `target/benchmarks/ecoli-10k-gfa-path-index.json`, with unchanged contigs, read containment, and
  0 assembly-only k-mers.
- The same indexed GFA path code cleared the E. coli 100k sentinel in
  `target/benchmarks/ecoli-100k-gfa-path-index.json`: 200,000 reads, 1,075,156 unique/trusted k-mers,
  163,132 unitigs, 163,124 contigs, N50 36 bp, 55.54 s GNU time, 1,752,444 KiB max RSS, 98.16%
  reliable-read k-mer containment, and 0 assembly-only k-mers. This proves the 100k scaling gate for
  the current representation; peak memory remains high and is not a memory-reduction claim.
- Streaming JSON sidecar writes removed the large pretty-JSON string buffers for evidence,
  annotations, simplification, scaffold, audit, and diploid reports. The E. coli 100k row improved to
  40.49 s / 1,586,644 KiB in `target/benchmarks/ecoli-100k-stream-json-sidecars.json` with identical
  contig count, total bases, reliable-read containment, and 0 assembly-only k-mers. This is a
  sidecar-allocation reduction only; the graph itself still needs a deeper memory pass.
- Post-assembly audit sidecars are now emitted without repair. The yeast BTT diploid 1k row wrote
  `target/benchmarks/yeast-btt-audit.json` and `target/benchmarks/yeast_btt_err1308583_diploid_1k_pairs/trex/audit.tsv`:
  1,885 contigs, 147,049 assembly k-mers, 147,031 trusted-supported k-mers, 18 low-support k-mers
  in one region, and one collapsed-repeat suspicion finding.
- Report-only mate endpoint joins exposed a Phase-2 diploid GFA metadata hotspot on E. coli 10k:
  `prepare_gfa_metadata` took 28,506 ms in
  `target/benchmarks/ecoli_mg1655_srr001666_10k_pairs/trex-mate-endpoint-joins/`. Indexing unitigs
  by first vertex and checking graph neighbors from each tail reduced the same stage to 24 ms in
  `target/benchmarks/ecoli_mg1655_srr001666_10k_pairs/trex-mate-endpoint-joins-fast-gfa/`; `graph.gfa`,
  `contigs.fa`, and `scaffolds.json` were byte-identical between the two runs.
