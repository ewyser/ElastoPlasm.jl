## Persistence (JLD2)

Simulation setup is saved/loaded as `ic["problem"]` (a `MechanicalProblem`, bundling
`mesh`/`mpts`/`time`) and `ic["basis"]` (its own top-level key, deliberately
independent of `Problem`); solver config as `cfg["solver"]` (NOT `cfg["instr"]` — that
key name is stale/wrong despite appearing in some scripts). `elastoplasm`/
`elastoplasm!` unpack `problem = file["ic/problem"]; mesh,mpts,time =
problem.mesh, problem.mpts, problem.time; basis = file["ic/basis"]`; `elastoplasm!`'s
write-back after each workflow rebuilds `file["ic/problem"]`/`file["ic/basis"]` in
place. Both take a `workflows::Vector{Function}` kwarg (plural — not `workflow`).
`mpts.s.cmp` (constitutive constants) persists as part of `ic["problem"].mpts`; there
is no separate `ic["cmpr"]` key. `test_column.jl` (formerly `test_collapse.jl`)
predated the current three-key layout for a long stretch (undefined `kwargser`, stale
`cfg["instr"]`, matrix-style `mpts.x` indexing) but was fully fixed — see
`.claude/bug/fixed/test-collapse-renamed-and-fixed.md` if picking that file up again ever
surfaces a similar staleness.
