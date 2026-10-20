## Persistence (JLD2)

A simulation's `.jld2` holds three keys: `"problem"`, the whole `Problem` (`mpts`, `mesh`,
`basis`, `time`, `solver`); `"paths"`, the output folders (`mkpaths`); `"misc"`. Write and read it
with `save_problem(file, problem; paths)`/`load_problem(file)` (`src/home/api/io/problem.jl`);
`setup_simulation` saves it, `elastoplasm`/`elastoplasm!` load it, and `elastoplasm!` writes
`"problem"` back after each workflow. `mpts.s.cmp` (constitutive constants) and APIC's `Bᵢⱼ`/`Dᵢⱼ`
(on `basis.transfer`) persist inside it. Before `refactor-init-problem-spec` the file had five keys
(`ic/problem` a `MechanicalProblem` of `mesh`/`mpts`/`time`, `ic/basis`, `cfg/solver`, `cfg/paths`,
`cfg/misc`); readers that used them (`elastoplasm`, `test_setup`, `test_column`, `test_basis`,
`test_performance`, `test_case`, `metanalysis.jl`) read `"problem"` now. (`cfg["instr"]`, seen in old
scripts, was never the key.) `test_column.jl` once lagged a layout change for a long stretch, see
`.claude/bug/fixed/test-collapse-renamed-and-fixed.md`.

**Format change:** stress/strain storage switched to a scalar plus a 3D Voigt deviator
(`feat-refactor-tensor-storage`), and the 2D `mpts.s.ϵzz` field was removed. `.jld2` files written before that don't load;
regenerate them with the example's `*_problem(...)` call.

**The problem as TOML beside the `.jld2`:** `setup_simulation` also writes `<name>_case.toml`
(`save_case`) next to `<name>_simulation.jld2`, so a run's setup can be read, edited, reloaded
(`load_case`) and rebuilt without the `.jld2` (provenance; skipped with a warning when a zone is a
`Region` predicate, which TOML cannot hold). The `.jld2` layout is unchanged.

**Format change:** `MeshSolidPhase` lost its unused `oobp` and `p` fields, `MeshProperties` and
`Point` their `nn`, and the file moved to the one `"problem"` key (`refactor-init-problem-spec`), so
`.jld2` files written before don't load; regenerate them with the
example's `*_problem(...)` call or `setup_simulation(load_case(<name>_case.toml))`.
