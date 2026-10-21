export Simulation

"""
    Simulation(problem; workflow, lstps=RUN.lstps, folder=RUN.folder)

A problem ready to run: the `problem`, the `workflow` it is run with (one of `WORKFLOWS`, the
quasi-static ones in `lstps` load steps; by default `RUN.workflow`, or `"elastoquasistatic!"` for an
implicit solver, which the explicit workflows cannot run) and the `folder` its files go to (a name under
`dump/`, or a path). A simulation's `.jld2` holds one (`save_simulation`/`load_simulation`); `execute!`
runs it.
"""
struct Simulation{P<:Problem}
    problem  :: P
    workflow :: String
    lstps    :: Int
    folder   :: String
end
function Simulation(problem::Problem; workflow::AbstractString=problem.solver isa ImplicitSolver ? "elastoquasistatic!" : RUN.workflow,
                    lstps::Integer=RUN.lstps, folder::AbstractString=RUN.folder)
    workflow ∈ WORKFLOWS || throw(ArgumentError("workflow \"$workflow\" is not one of $(WORKFLOWS)"))
    quasistatic(workflow) || problem.solver isa ExplicitSolver ||
        throw(ArgumentError("workflow \"$workflow\" needs an explicit solver (solution = \"explicit\")"))
    lstps ≥ 1 || throw(ArgumentError("lstps ≥ 1 is needed"))
    dir = rstrip(isabspath(folder) ? normpath(folder) : joinpath(self.sys.dump, folder), ['/', '\\'])
    return Simulation(problem, String(workflow), Int(lstps), String(dir))
end
