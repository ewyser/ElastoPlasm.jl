export setup_simulation, case_workflows

"""
    Solver(case; run=(;), overrides…) -> ExplicitSolver or ImplicitSolver

How a case is solved, as its settings say: the numerical options (the option sections that are not physics: `dtype`,
`basis`, `stab`, `nonloc`, `solution`, `plot`, `output`, `backend`) and the run settings (`workflow`, `te`, `tg`,
`tep`, `lstps`), with the kernels built on them; `solution` follows the workflow (`"implicit"` for the quasi-static
ones) unless set. The solver keeps nothing of the case: `Problem(case, solver)` builds the two together, and the same
solver solves other cases of its dimension (`Problem(other, solver)`).

`overrides` replace numerical options field by field (`Solver(case; basis=(; which="gimpm"))`) and `run` run settings
(`run=(; te=5.0)`), so one case is solved several ways.

# Throws
- `ArgumentError`: when an override is physics (a `get_physics_sections()` section, or `bcs`: set it in the case), or
  the run folder (the case's, or `setup_simulation`'s `folder`), or an explicit workflow is asked of an implicit solver
  (`solution = "implicit"`); and for an unknown or invalid option or run setting.

# Examples
```julia
case = slump_case([64.16, 16.04], [40, 10])
a = Solver(case)                                   # as its settings say
b = Solver(case; basis=(; which="gimpm"), run=(; te=5.0))
execute!(setup_simulation(case, b))
```
"""
function Solver(case::Case; run=(;), overrides...)
    o = named_tuple((; overrides...))
    refuse_physics(o)
    haskey(config_dict(run), "folder") &&
        throw(ArgumentError("run.folder: a simulation's folder is its case's (its run settings) or setup_simulation's `folder`"))
    s = case.settings
    options = isempty(o) ? s.options : merge_dicts(s.options, option_labels(config_dict(o)))
    settings = merge(Dict{String,Any}(s.run), config_dict(run))
    delete!(settings, "thermal")                             # the case's
    r = solver_run(settings)
    config = (; (k => v for (k, v) ∈ pairs(case_options(options)) if k ∉ get_physics_sections())...)
    haskey(config, :solution) || (config = merge(config, (; solution=quasistatic(r.workflow) ? "implicit" : "explicit")))
    quasistatic(r.workflow) || config.solution == "explicit" ||
        throw(ArgumentError("workflow \"$(r.workflow)\" needs an explicit solver (solution = \"explicit\")"))
    D = length(case.layout.L)
    return build_solver(merge(solver_instr(; dim=D, config...), (; run=r)), D)
end
# solver options that are a case's physics are refused: set them in the case
function refuse_physics(options)
    physics = [k for k ∈ keys(options) if k ∈ (get_physics_sections()..., :bcs)]
    isempty(physics) || throw(ArgumentError("$(join(physics, ", ")): the case's physics; set it in the case (its settings)"))
    return nothing
end
# `b`'s entries over `a`'s, nested tables merged entry by entry
function merge_dicts(a::AbstractDict, b::AbstractDict)
    out = Dict{String,Any}(a)
    for (k, v) ∈ b
        out[k] = haskey(out, k) && out[k] isa AbstractDict && v isa AbstractDict ? merge_dicts(out[k], v) : v
    end
    return out
end
# the fields of a configuration that differ from `ref`'s, sections recursively (unchanged sections left out)
function config_changes(c::NamedTuple, ref::NamedTuple)
    out = Pair{Symbol,Any}[]
    for k ∈ keys(c)
        v, r = getfield(c, k), get(ref, k, nothing)
        if v isa NamedTuple && r isa NamedTuple
            d = config_changes(v, r)
            isempty(d) || push!(out, k => d)
        elseif !isequal(v, r)
            push!(out, k => v)
        end
    end
    return (; out...)
end

"""
    setup_simulation(case, solver; folder=case_folder(case)) -> String
    setup_simulation(case) -> String

Build a case's simulation, `Simulation(Problem(case, solver); folder)`, and save it (`save_simulation`) with the case
as it was solved beside it (`save_case(file, case, solver)`, `<name>_case.toml`), the initial fields plotted when
`solver.plot.status`; removes an output (`output_file(sim)`) left there by an earlier simulation, and returns the
`.jld2` path, for `execute!`. The folder is the case's (its run settings' `folder`, else its file's folder), a name
under `dump/` or a path. `setup_simulation(case)` is `setup_simulation(case, Solver(case))`.
"""
function setup_simulation(case::Case, solver::Solver; folder::AbstractString=case_folder(case))
    @info "Setting up mesh & material point system for a $(length(case.layout.L))d problem ($(length(case.layout.zones)) zones)"
    sim = Simulation(Problem(case, solver); folder)
    p   = sim.problem
    @info ic_log(p.mesh, p.mpts, p.time, p.solver)
    mkpath(sim.folder)
    p.solver.plot.status && plot_initial_fields(p, sim.folder)
    try
        save_case(joinpath(sim.folder, "$(basename(sim.folder))_case.toml"), case, solver; folder)
    catch e
        @warn "Case not written as TOML: $(sprint(showerror, e))"
    end
    rm(output_file(sim); force=true)                         # a new simulation starts a new output
    return save_simulation(sim)
end
function setup_simulation(case::Case)
    return setup_simulation(case, Solver(case))
end

"""
    case_workflows(sim) -> Vector{Function}
    case_workflows(workflow, lstps) -> Vector{Function}

The workflow functions of a workflow name (one of `get_workflows()`), the quasi-static ones in `lstps` load
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
