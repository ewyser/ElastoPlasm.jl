# Ongoing task state — typed strain/stress tensor work & follow-ups

Local, gitignored (as of this session — `.claude/` is now in the repo's own
`.gitignore`, not just the user's global git ignore, which previously only covered
`.claude/settings.local.json` specifically) trace of in-progress Claude Code work on
this repo. Check this file, and `.claude/plans.md` (queued-but-not-started plans), at
the start of a session here.

**Note on this file's history**: an earlier, more detailed version of this file
existed earlier in this same session (tracking the typed-tensor port round by round)
but went missing partway through — found gone, cause unknown, when checked again
later in the session. Given `.claude/` was never actually gitignored until this
rewrite (see above), it's possible something touched it; no direct evidence found.
Rebuilt from conversation context below rather than left blank.

## Done, merged to `dev`

The full typed strain/stress tensor port and its follow-up cleanups, on branch
`feat-typed-strain-stress-tensors` (5 commits: `a026b4f` tensor port,
`fbc64ec` Δλ fixes, `dfc2811` `AbstractElasticity` removal, `67e5297` `get_voigt`
convention fix + `mutate`/`_mutate` retirement, `941e8d4` constructor/comment
cleanup) — merged into `mimic-implementation-logic-of-ample`, then into `dev`
(merge commits `ff6134a`, `78838bf`), and `mimic` fast-forwarded to match `dev`. One
follow-up doc commit landed directly on `mimic`/`dev` after the merge: `2d474c5`
(`docs: Confirm test_workflow.jl baseline holds after tensor-port merge to dev`).

All of this is independently verified (not just self-reported) at multiple points
this session: bit-exact/machine-epsilon numerical equivalence checks against pre-change
math, zero-allocation checks, the full `test_workflow.jl` 96-case sweep unchanged at
71 passed/25 failed throughout, `test_collapse.jl` 4/4, `test_basis.jl` passing,
`test_performance.jl` showing real improvement (-17.2% memory/-5.4% allocs/0.0% time
vs. the `dev`-generated baseline) — see CLAUDE.md's "Typed strain/stress tensor
storage" and related entries for the full detail, since those are committed and won't
go missing the way this file did.

`fancy-implementation` (the source branch these ideas were ported from) still exists
locally and on `origin` — user was weighing deleting it since everything useful from
it is now ported; no decision recorded on whether/when to actually delete.

## Queued, not started

**See `.claude/plans.md`** for a fully-detailed, approved-but-not-implemented plan:
dispatch transfer scheme (`std`/`tpic`/`apic`) via a new `Basis` type parameter
(mirroring `Basis`'s existing `kind::K<:AbstractBasis` pattern), replacing the
`std_p2n`/`tpic_p2n`/`apic_p2n` separately-named-function + `Cairn`-closure-table
dispatch with ordinary Julia multiple dispatch. Includes an addendum: `Point`'s
`Δnp` field is dead code (delete it), and `Bᵢⱼ`/`Dᵢⱼ` (APIC-only, currently always
allocated on every `Point` regardless of transfer scheme) move onto a new
`ApicTransfer` type. Start a new branch off current `dev` when picking this up.

## Standing session notes

- User prefers `function ... end` block syntax over compact `f(...) = expr` Julia
  syntax (saved as its own memory: `julia-function-syntax-style`).
- User prefers saving ongoing-task-state traces in this repo's own `.claude/` folder
  in addition to the cross-session global memory (saved as its own memory, on the
  `track-ongoing-task-state` file).
- Watch for subagents spawning further subagents instead of doing delegated work
  directly and reporting real numbers — happened repeatedly earlier this session;
  redirect explicitly and check `ListAgents` before trusting a "done" report.
