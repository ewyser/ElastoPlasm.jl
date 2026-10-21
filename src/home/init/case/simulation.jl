export setup_simulation, run_simulation, case_workflows

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
    setup_simulation(case; fid=case.run.fid) -> String

Build the `Problem` of a case (`Problem(case)`) and save it for `elastoplasm`/`elastoplasm!`, the
initial fields plotted when `solver.plot.status`. Writes `<name>_simulation.jld2` (`save_problem`)
with the case's TOML (`save_case`) beside it into the case's folder `run.dir` (`name` its last part),
or, when that is empty, into `dump/<fid>/` (`name` the last part of `fid`); returns the `.jld2` path.
"""
function setup_simulation(s::Case; fid::AbstractString=s.run.fid)
    @info "Setting up mesh & material point system for a $(length(domain(s)[1]))d problem ($(length(s.zones)) zones)"
    problem = Problem(s)
    dir     = rstrip(normpath(abspath(s.run.dir)), ['/', '\\'])   # "…/case/" → "…/case"
    paths   = isempty(s.run.dir) ? mkpaths(fid, self.sys.dump; interactive=false) :
                                   mkpaths(basename(dir), dirname(dir); interactive=false)
    problem.solver.plot.status && plot_initial_fields(problem, paths)
    name    = isempty(s.run.dir) ? basename(fid) : basename(dir)
    try
        save_case(joinpath(paths[:dat], "$(name)_case.toml"), s)
    catch e
        @warn "Case not written as TOML: $(sprint(showerror, e))"
    end
    return save_problem(joinpath(paths[:dat], "$(name)_simulation.jld2"), problem; paths)
end

"""
    case_workflows(case) -> Vector{Function}

The workflows of a case's run (`case.run.workflow`, one of `WORKFLOWS`), for
`elastoplasm!(jld2; workflows=...)`; the quasi-static ones with `case.run.lstps` load steps.
"""
function case_workflows(s::Case)
    w = s.run.workflow
    w == "elastodynamic! + elastoplastic!" && return Function[elastodynamic!, elastoplastic!]
    w == "elastodynamic!"                  && return Function[elastodynamic!]
    w == "elastoplastic!"                  && return Function[elastoplastic!]
    plastic, lstps = w == "elastoquasistatic! (plastic)", s.run.lstps
    return Function[(a...) -> elastoquasistatic!(a...; lstps=lstps, plastic=plastic)]
end

"""
    run_simulation(case; kwargs...) -> NamedTuple
    run_simulation(file; kwargs...)

Build and run a case as its run settings say: `elastoplasm!(setup_simulation(case; kwargs...);
workflows=case_workflows(case))`. Takes a case file (`load_case`) as well.
"""
run_simulation(s::Case; kwargs...) = elastoplasm!(setup_simulation(s; kwargs...); workflows=case_workflows(s))
run_simulation(file::AbstractString; kwargs...) = run_simulation(load_case(file); kwargs...)
