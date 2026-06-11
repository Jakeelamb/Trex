# Trex fuzzing (Phase-1 + optional Phase-2 hooks)

Trex FASTQ ingest is exercised with **`cargo-fuzz`** (install: `cargo install cargo-fuzz`).

## Targets

- `parse_fastq` — `trex::illumina::fastq::parse_fastq` on arbitrary bytes.

## Phase-2 Illumina (optional / scheduled)

Additional targets (e.g. **GFA** path round-trips, diploid export parsers) are **not** PR-blocking per **CONTEXT: Phase-2 Illumina fuzzing**; add bins under `fuzz/fuzz_targets/` when a runner budget exists.

## Run

```bash
cd fuzz
cargo fuzz run parse_fastq corpus/parse_fastq -- -runs=1000
```

Seeds live under `fuzz/corpus/` (small, version-controlled).
