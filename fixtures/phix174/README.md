# PhiX174 real-reference benchmark fixture

This fixture is the first Trex governed row using a real public reference.

- `reference.fa`: NCBI RefSeq `NC_001422.1` (`Escherichia phage phiX174, complete genome`), fetched through NCBI E-utilities.
- `reads.fq`: deterministic 150 bp single-end synthetic reads generated from `reference.fa` every 50 bp with circular wraparound.
- Digests are pinned in `tools/manifest.toml` under `[fixtures.phix174]`.

This is a micro benchmark for reproducible direct Trex metrics, not a biological readset benchmark: the reference is real, while the reads are deterministic synthetic reads generated from that reference.

Regenerate reads:

```bash
cargo run -p xtask -- generate-reads \
  --reference fixtures/phix174/reference.fa \
  --out fixtures/phix174/reads.fq \
  --read-len 150 \
  --step 50 \
  --circular
```

Run the row:

```bash
cargo run -p xtask -- bench --tier nightly --out target/benchmarks/nightly.json
```
