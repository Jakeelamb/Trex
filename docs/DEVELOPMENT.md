# Trex Development Protocol

This document is the durable handoff contract for Trex development agents. It is
validated by `cargo run -p xtask -- validate-development`.

## Scope

Trex is a Rust-first, Illumina-first assembler. The active product contract is
the Phase-2 Illumina endgame on top of the existing Phase-1 compatibility layer.
Long-read, hybrid, HiFi, ONT, Hi-C, and T2T-style work may shape interfaces, but
they remain deferred until an ADR changes the active scope.

## Orchestrator

The orchestrator owns repository proof and integration. Only the orchestrator may:

- Run `git`, inspect branch state, create commits, push branches, or manage remotes.
- Run `cargo`, `xtask`, benchmark scripts, formatters, clippy, and CI-equivalent gates.
- Edit the main repo worktree.
- Create, assign, merge, or delete temporary worker worktrees.
- Decide whether a worker packet is accepted, rejected, or needs another packet.
- Maintain the docs and ledgers that say what is proven, planned, or deferred.

Every behavioral development wave must end with these gates unless the final
handoff names the exact blocker:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-features -- -D warnings
cargo test --workspace --all-features
cargo run -p xtask -- validate
```

Run a subsystem benchmark gate when the change touches a benchmarked surface.

## Worker Rules

Workers operate in one of two modes assigned by the orchestrator.

Read-only workers provide constrained packets to the orchestrator and do not
mutate any worktree.

Implementation workers may edit only the orchestrator-assigned temporary
worktree, only inside their declared file ownership, and only for the assigned
task. They do not integrate their own changes into the main repo.

Allowed:

- Read files with `rg`, `sed`, `nl`, `ls`, `wc`, or equivalent read-only commands.
- Inspect docs, source files, tests, manifests, literature notes, and benchmark artifacts.
- Return concise findings, interface proposals, test proposals, and patch sketches.
- For implementation workers only: edit assigned files in the assigned temporary worktree.

Forbidden:

- No `git`.
- No `cargo`.
- No edits outside the assigned file ownership.
- No formatting.
- No benchmark execution.
- No broad repo rewrites.
- No generated dependency, cache, or benchmark artifact churn.
- No claims without file references or paper-derived rationale.

## Worker Packet Template

```markdown
# Worker Packet Result

Objective:
Files Read:
Paper/Technique Basis:
Findings:
Proposed Change:
Patch Sketch:
Tests Orchestrator Should Run:
Risks / Unknowns:
```

Worker packets should be narrow: one objective, enough exact paths and line
references to verify the claim, and no attempt to make policy decisions for the
orchestrator.

## Claim Levels

| Level | Meaning | Required evidence |
|-------|---------|-------------------|
| `observed` | A fact was seen in source, docs, logs, or an artifact. | File path, command output, artifact path, or literature note. |
| `tested` | A behavior is covered by a focused automated test. | Test name plus the command that ran it. |
| `benchmarked` | A behavior or resource claim has a stored benchmark artifact. | Matrix row, command, JSON/script artifact, and tier. |
| `production-gated` | A claim is enforced by the normal Trex validation or CI path. | `xtask validate`, `xtask gate`, CI job, or documented release gate. |

Do not promote a claim to a stronger level because it is plausible. Promote it
only when the matching evidence exists.

## Acceptance Checklist

Before accepting a worker packet or a local implementation, the orchestrator
checks:

- The change stays inside the active Illumina contract or is clearly marked deferred.
- The change has a test, an `xtask` validation check, a benchmark artifact, or an explicit deferred label.
- Public behavior changes are reflected in `README.md`, `CONTEXT.md`, `ARCHITECTURE.md`, `docs/CAPABILITIES.md`, an ADR, or this document as appropriate.
- Benchmark or quality claims name the exact row, command, and artifact.
- Missing optional tools are reported as unavailable instead of hidden in logs.

## Commit Discipline

Use one coherent behavioral change per commit. Documentation must move with
public behavior. Performance commits need before/after artifacts. Research or
future adapter commits must not change active Illumina behavior unless an ADR
has already changed the scope.
