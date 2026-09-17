# `solver.plast.status` is dead config in the explicit solver path — it gates nothing

**Status: open** (design decision needed, not a crash).

Found while investigating 3D conformity: `elastoplastic!` (one of the two
built-in explicit workflows) calls `elastoplast()`
(`src/home/core/solver/explicit/update/update.jl`), which calls
`solver.cairn.update.retmap!(...)` (the plastic return-mapping dispatcher)
**unconditionally** — there is no `if solver.plast.status` branch at that
call site, nor anywhere else in `src/home/`; the only other reference to
`plast.status` in the entire directory is a docstring comment
(`collapse.jl`). Whether plasticity actually runs is decided entirely by
*which workflow function you pass to `elastoplasm!`* — `elasto`/
`elastodynamic!` never touches `retmap!`, `elastoplast`/`elastoplastic!`
always does — the config flag plays no role. Concretely: `test_workflow.jl`'s
sweep (and the 3D sweep run against it, see
`.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`) calls
`elastoplasm!(jld2; workflows=[elastodynamic!, elastoplastic!])` with
`plast.status` left at its `false` default, and plasticity genuinely runs
anyway — so that sweep *does* exercise the plastic corrector, not just the
elastic path.

Contrast with `nonloc.status`, a same-shaped, similarly-named flag that
**is** live (gates the nonlocal regularization block inside `elastoplast`
itself) — the two flags look parallel but behave completely differently,
which is exactly the kind of thing worth checking directly rather than
assuming from the name.

Not fixed: either `plast.status` should actually gate `retmap!`'s call
(matching what its name implies), or it should be removed/documented as
vestigial if `elasto`/`elastoplast` being separate workflow functions is
considered sufficient by design.
