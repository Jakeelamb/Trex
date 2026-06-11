# ADR 0002 — Promote active scope to the Phase-2 Illumina endgame

## Status

Accepted.

## Context

Trex shipped an Illumina assembler with Phase-1 gates and an experimental Phase-2 diploid overlay. The code path already supports `trex illumina assemble --diploid`, Phase-2 GFA annotations, primary FASTA collapse, mate-pair edge boosts, and synthetic two-parent benchmark gates.

## Decision

Trex active engineering is promoted from Phase-1-only Illumina stability to the **Phase-2 Illumina endgame**: Illumina-only eukaryotic diploid assembly on the single `trex illumina assemble` path, with Phase-1 retained as a non-relaxing compatibility and benchmark layer.

## Consequences

Long-read and hybrid assembly remain explicitly deferred so the promotion deepens the Illumina assembler instead of opening every Phase-2 track at once.
