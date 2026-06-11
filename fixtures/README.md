# Fixtures (Phase-1)

| Path | Role |
|------|------|
| [`tiny.fq`](tiny.fq) | Single-end Illumina-style read (48 bp) for smoke tests and `scripts/ref_free_smoke.sh`. |
| [`tiny_ref.fa`](tiny_ref.fa) | Same sequence as the assembled contig; used by `scripts/reference_smoke.sh` (minimap2 or substring check). |
| [`expected/ref_free_smoke/`](expected/ref_free_smoke/) | Golden FASTA/GFA from `trex illumina assemble` on `tiny.fq` with `--kmer-size 4` and `--trusted-threshold 1`. The benchmark gate compares live `target/ref-free-smoke/` output to these files. `graph.gfa` includes **`P`** lines when the contig path matches the unitig path. |
| [`phase2_synthetic/`](phase2_synthetic/) | Two parental **32 bp** haplotypes + `reads.fq` for the Phase-2 Illumina gate; digests in **`tools/manifest.toml`** under `[fixtures.phase2_synthetic]` (see [`phase2_synthetic/README.md`](phase2_synthetic/README.md)). |
| [`phix174/`](phix174/) | Real public RefSeq `NC_001422.1` reference plus deterministic circular 150 bp synthetic reads for the nightly/manual direct Trex benchmark row; digests in **`tools/manifest.toml`** under `[fixtures.phix174]`. |

## Refresh expected outputs

After an intentional assembly/export change:

```bash
cargo run -p trex-cli --release -- illumina assemble \
  --r1 fixtures/tiny.fq --kmer-size 4 --trusted-threshold 1 \
  --out-dir fixtures/expected/ref_free_smoke
```

Then commit the updated files under `fixtures/expected/ref_free_smoke/`.
