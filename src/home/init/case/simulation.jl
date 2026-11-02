export setup_simulation, case_workflows

"""
    Solver(case; run=(;), overrides…) -> ExplicitSolver or ImplicitSolver

How `case` is solved: the numerics read from its settings (the option sections that are not physics: `dtype`,
`basis`, `stab`, `nonloc`, `solution`, `plot`, `output`, `backend`) and its run settings (`workflow`, `te`, `tg`,
`tep`, `lstps`, `folder`), with the kernels built on them. The physics comes from the case (`case.options`, the
boundary conditions); `solution` follows the workflow (`"implicit"` for the quasi-static ones) unless set.

`overrides` replace numerical options field by field (`Solver(case; basis=(; which="gimpm"))`) and `run` run settings
(`run=(; te=5.0)`), so one case is solved several ways without rebuilding its physics. The solver keeps the case it
solves, its settings with these overrides (`solver.case`), which `Problem(solver)` builds from and `setup_simulation`
writes beside the run.

# Throws
- `ArgumentError`: when an override is physics (a `PHYSICS_SECTIONS` section, or `bcs`: set it in the case), or an
  explicit workflow is asked of an implicit solver (`solution = "implicit"`); and for an unknown or invalid option or
  run setting.

# Examples
```julia
case = slump_case([64.16, 16.04], [40, 10])
a = Solver(case)                                   # as its settings say
b = Solver(case; basis=(; which="gimpm"), run=(; te=5.0))
execute!(setup_simulation(b))
```
"""
function Solver(case::Case; run=(;), overrides...)
    o = named_tuple((; overrides...))
    physics = [k for k ∈ keys(o) if k ∈ (PHYSICS_SECTIONS..., :bcs)]
    isempty(physics) || throw(ArgumentError("$(join(physics, ", ")): the case's physics; set it in the case (its settings)"))
    c = with_numerics(case, o, run)
    r = case_run(c.layout, c.settings)
    config = merge(case_options(c.settings.options), (; bcs=(; dirichlet=c.bcs)))
    haskey(config, :solution) || (config = merge(config, (; solution=quasistatic(r.workflow) ? "implicit" : "explicit")))
    quasistatic(r.workflow) || config.solution == "explicit" ||
        throw(ArgumentError("workflow \"$(r.workflow)\" needs an explicit solver (solution = \"explicit\")"))
    D = length(c.layout.L)
    return build_solver(merge(solver_instr(; dim=D, config...), (; run=r)), D; case=c)
end
# the case as a solver solves it: its settings with the numerical overrides (in their file form) and run settings
function with_numerics(case::Case, overrides::NamedTuple, run)
    isempty(overrides) && isempty(run) && return case
    s = case.settings
    options = merge_dicts(s.options, option_labels(config_dict(overrides)))
    settings = EasyDraw.Settings(; zones=s.zones, defaults=s.defaults, sides=s.sides, options,
                                   run=merge(s.run, Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(run))))
    return Case(case.layout, settings)
end
# `b`'s entries over `a`'s, nested tables merged entry by entry
function merge_dicts(a::AbstractDict, b::AbstractDict)
    out = Dict{String,Any}(a)
    for (k, v) ∈ b
        out[k] = haskey(out, k) && out[k] isa AbstractDict && v isa AbstractDict ? merge_dicts(out[k], v) : v
    end
    return out
end

"""
    setup_simulation(solver) -> String
    setup_simulation(case) -> String

Build a solver's simulation, `Simulation(Problem(solver))` in its run folder, and save it (`save_simulation`) with
its case's TOML (`save_case`, the settings it was solved with) beside it, the initial fields plotted when
`solver.plot.status`; removes an output (`output_file(sim)`) left there by an earlier simulation, and returns the
`.jld2` path, for `execute!`. `setup_simulation(case)` is `setup_simulation(Solver(case))`.
"""
function setup_simulation(solver::AbstractSolver)
    s = solver.case
    isnothing(s) && throw(ArgumentError("this solver solves no case (get_solver): build one with Solver(case)"))
    @info "Setting up mesh & material point system for a $(length(s.layout.L))d problem ($(length(s.layout.zones)) zones)"
    sim = Simulation(Problem(solver))
    p   = sim.problem
    @info ic_log(p.mesh, p.mpts, p.time, p.solver)
    mkpath(sim.folder)
    p.solver.plot.status && plot_initial_fields(p, sim.folder)
    try
        save_case(joinpath(sim.folder, "$(basename(sim.folder))_case.toml"), s)
    catch e
        @warn "Case not written as TOML: $(sprint(showerror, e))"
    end
    rm(output_file(sim); force=true)                         # a new simulation starts a new output
    return save_simulation(sim)
end
function setup_simulation(case::Case)
    return setup_simulation(Solver(case))
end

"""
    case_workflows(sim) -> Vector{Function}
    case_workflows(workflow, lstps) -> Vector{Function}

The workflow functions of a workflow name (one of `WORKFLOWS`), the quasi-static ones in `lstps` load
steps: what `execute!` runs by default, for a simulation its solver's (`solver.run`).
"""
function case_workflows(sim::Simulation)
    run = sim.problem.solver.run
    return case_workflows(run.workflow, run.lstps)
end
function case_workflows(w::AbstractString, lstps::Integer)
    w == "elastodynamic! + elastoplastic!" && return Function[elastodynamic!, elastoplastic!]
    w == "elastodynamic!"                  && return Function[elastodynamic!]
    w == "elastoplastic!"                  && return Function[elastoplastic!]
    plastic = w == "elastoquasistatic! (plastic)"
    return Function[(a...) -> elastoquasistatic!(a...; lstps, plastic)]
end
