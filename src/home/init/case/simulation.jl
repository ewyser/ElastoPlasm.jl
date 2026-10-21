export setup_simulation, case_workflows

"""
    case_solver_options(case) -> NamedTuple

The `get_solver` keyword arguments of a case: its solver overrides (merged per field over the
defaults by `get_solver`), its boundary conditions as `bcs.dirichlet` (roller on the sides along y of
an extruded problem), and `solution` from the workflow (`"implicit"` for the quasi-static ones) unless
the overrides set it.
"""
function case_solver_options(s::Case)
    dir = isnothing(s.width) ? s.bcs : vcat(s.bcs[1:1, :], [:roller :roller], s.bcs[2:2, :])
    kw  = merge(s.solver, (; bcs=(; dirichlet=dir)))
    haskey(kw, :solution) || (kw = merge(kw, (; solution=quasistatic(s.run.workflow) ? "implicit" : "explicit")))
    return kw
end

"""
    setup_simulation(case) -> String

Build a case's simulation, `Simulation(Problem(case); workflow, lstps, folder)` from its run settings,
and save it (`save_simulation`) with the case's TOML (`save_case`) beside it in its folder, the initial
fields plotted when `solver.plot.status`; returns the `.jld2` path, for `execute!`.
"""
function setup_simulation(s::Case)
    @info "Setting up mesh & material point system for a $(length(domain(s)[1]))d problem ($(length(s.zones)) zones)"
    sim = Simulation(Problem(s); s.run.workflow, s.run.lstps, s.run.folder)
    p   = sim.problem
    @info ic_log(p.mesh, p.mpts, p.time, p.solver)
    mkpath(sim.folder)
    p.solver.plot.status && plot_initial_fields(p, sim.folder)
    try
        save_case(joinpath(sim.folder, "$(basename(sim.folder))_case.toml"), s)
    catch e
        @warn "Case not written as TOML: $(sprint(showerror, e))"
    end
    return save_simulation(sim)
end

"""
    case_workflows(sim) -> Vector{Function}
    case_workflows(workflow, lstps) -> Vector{Function}

The workflow functions of a workflow name (one of `WORKFLOWS`), the quasi-static ones in `lstps` load
steps: what `execute!` runs by default.
"""
case_workflows(sim::Simulation) = case_workflows(sim.workflow, sim.lstps)
function case_workflows(w::AbstractString, lstps::Integer)
    w == "elastodynamic! + elastoplastic!" && return Function[elastodynamic!, elastoplastic!]
    w == "elastodynamic!"                  && return Function[elastodynamic!]
    w == "elastoplastic!"                  && return Function[elastoplastic!]
    plastic = w == "elastoquasistatic! (plastic)"
    return Function[(a...) -> elastoquasistatic!(a...; lstps, plastic)]
end
