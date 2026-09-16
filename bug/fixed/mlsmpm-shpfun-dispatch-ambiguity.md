# `test_workflow.jl`'s 96-case sweep reported 47/49 instead of the documented 71/25 baseline

**Status: fixed.**

Root cause: **every one of `mlsmpm`'s 24 cases was throwing `MethodError: ...
is ambiguous`** — a genuine `shpfun!` dispatch ambiguity between the generic
method (`src/home/core/common/shpfun.jl`) and `mlsmpm.jl`'s MLS-specific
method, introduced by a commit that narrowed the generic method's `Point`
pattern from `Point{T1,T2,D,CM}` (4 explicit type params) to
`Point{T1,T2,D}` (3 explicit). Both patterns match the exact same set of
concrete `Point` types (`CM` unconstrained either way) — but Julia's
method-specificity algorithm empirically treats an explicitly-named-but-free
type parameter differently from an implicit "rest" pattern when compared
against a sibling method (`mlsmpm.jl`'s) constraining a *different* parameter
(`Basis{...,K<:MLSBasis}`) — a real, obscure dispatch quirk, not reviewable
by eye (confirmed via `Test.detect_ambiguities` on an isolated reproduction).
Fixed originally by naming `CM` explicitly again; later made structurally
robust by constraining the generic method's `K` to
`Union{BSplineBasis,GimpBasis,LinearBasis}` (explicitly excluding
`MLSBasis`) instead of relying on the naming quirk — the two `shpfun!`
methods' `Basis` type patterns are now disjoint by construction. Verified:
sweep exactly matches the documented 71/25 baseline, identical per-basis-kind
breakdown (smpm 12 fail, gimpm 13 fail, bsmpm 0 fail, mlsmpm 0 fail — see
`bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md` for that
baseline's own history); zero `shpfun!` ambiguities via
`Test.detect_ambiguities`.

## Methodology lesson from chasing this

An earlier pass at this investigation wrongly concluded
`elastodynamic!`/`elastoplastic!` produce exactly zero velocity/stress for
any `slump_problem` run, treating it as a separate severe bug. That was
**not a real bug** — it was caused by calling `elastoplasm(...)` (no `!`)
and then re-loading `mpts` from the JLD2 file afterward. `elastoplasm` opens
the file read-only and never writes the post-workflow state back — only
`elastoplasm!` does — so the re-load always returned the untouched,
zero-initialized pre-simulation state. **Always use `elastoplasm!`, never
`elastoplasm`, when the goal is to inspect post-run `mpts`/`mesh` state from
the JLD2 file afterward** — `elastoplasm`'s return value doesn't carry the
mutated objects either (just a bare `(;simulation,success)`), so there is no
way to recover real results from it short of re-reading the file.
