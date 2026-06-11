# Trex scripts

| Script | Purpose |
|--------|---------|
| [`ref_free_smoke.sh`](ref_free_smoke.sh) | **Phase-1 PR smoke**: assemble `fixtures/tiny.fq`, print reference-free stats on FASTA, compare `contigs.fa` / `unitigs.fa` / `graph.gfa` to `fixtures/expected/ref_free_smoke/`. |
| [`reference_smoke.sh`](reference_smoke.sh) | **Phase-1 benchmark gate** (reference leg): check contigs from `ref_free_smoke` against `fixtures/tiny_ref.fa` via **minimap2** PAF or substring fallback. Expects `target/ref-free-smoke/contigs.fa` (or `REF_FREE_SMOKE_OUT`). |
| [`benchmark_gate.sh`](benchmark_gate.sh) | Full **Phase-1 benchmark gate**: `ref_free_smoke.sh` then `reference_smoke.sh`. |
| [`phase2_illumina_diploid_reference_layer.sh`](phase2_illumina_diploid_reference_layer.sh) | **Phase-2 Illumina** diploid fixture leg: SHA-256 vs [`tools/manifest.toml`](../tools/manifest.toml), parental/read consistency checks, optional minimap2. |
| [`phase2_illumina_graph_summaries.sh`](phase2_illumina_graph_summaries.sh) | **Phase-2 Illumina graph summaries**: `trex illumina assemble --diploid` on `fixtures/phase2_synthetic/reads.fq`, stats on `contigs.fa`, GFA record counts + `trex-phase2-illumina` header tag. |
| [`phase2_illumina_benchmark_gate.sh`](phase2_illumina_benchmark_gate.sh) | Layered **Phase-2 Illumina benchmark gate**: Phase-1 gate → diploid reference layer → graph summaries. |

CI: pull requests run `ref_free_smoke`, both Phase-2 scripts, and `cargo test`; **minimap2** + `reference_smoke` run on `main`/`master`, tags, `schedule`, and `workflow_dispatch` (see [`.github/workflows/ci.yml`](../.github/workflows/ci.yml)).
