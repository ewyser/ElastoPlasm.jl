## Persistence (JLD2)

A simulation's `.jld2` is a saved `Simulation`, one `"simulation"` key: the `Problem` (`mpts`,
`mesh`, `basis`, `time`, `solver`), the workflow name and load steps, and the folder. It is written
to `<folder>/<name>_simulation.jld2` by `save_simulation(sim)` (or to a given file,
`save_simulation(sim, file)`, the one writer of the key) and read by `load_simulation(file)`
(`src/home/api/io/simulation.jl`), whose folder is then the file's (a run folder can be moved or
copied); `setup_simulation` saves it, and `execute!(jld2)` loads it and writes it back after each
workflow. `mpts.s.cmp` (constitutive constants) and APIC's `Bᵢⱼ`/`Dᵢⱼ` (on
`basis.transfer`) persist inside it. **Format change:** before the `Simulation` type the file held
`"problem"`, `"paths"` and `"misc"` (`save_problem`/`load_problem`), and before
`refactor-init-problem-spec` five keys (`ic/problem`, `ic/basis`, `cfg/solver`, `cfg/paths`,
`cfg/misc`); such files don't load (`load_simulation` says to set them up again from the case TOML). (`cfg["instr"]`, seen in old
scripts, was never the key.) `test_column.jl` once lagged a layout change for a long stretch, see
`.claude/bug/fixed/test-collapse-renamed-and-fixed.md`.

**Format change:** stress/strain storage switched to a scalar plus a 3D Voigt deviator
(`feat-refactor-tensor-storage`), and the 2D `mpts.s.ϵzz` field was removed. `.jld2` files written before that don't load;
regenerate them with `setup_simulation(<example>_case(...))`.

**The problem as TOML beside the `.jld2`:** `setup_simulation` also writes `<name>_case.toml`
(`save_case`) next to `<name>_simulation.jld2`, so a run's setup can be read, edited, reloaded
(`load_case`) and rebuilt without the `.jld2` (provenance; skipped with a warning when a zone is a
`Region` predicate, which TOML cannot hold). The `.jld2` layout is unchanged.

**Format change:** `MeshSolidPhase` lost its unused `oobp` and `p` fields, `MeshProperties` and
`Point` their `nn`, and the file moved to the one `"problem"` key (`refactor-init-problem-spec`), so
`.jld2` files written before don't load; regenerate them with
`setup_simulation(<example>_case(...))` or `setup_simulation(load_case(<name>_case.toml))`.
