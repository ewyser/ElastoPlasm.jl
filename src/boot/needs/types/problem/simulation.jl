export Simulation

"""
    Simulation(problem; workflow=RUN.workflow, lstps=RUN.lstps, folder=RUN.folder)

A problem ready to run: the `problem`, the `workflow` it is run with (one of `WORKFLOWS`, the
quasi-static ones in `lstps` load steps) and the `folder` its files go to (a name under `dump/`, or a
path). A simulation's `.jld2` holds one (`save_simulation`/`load_simulation`); `execute!` runs it.
"""
struct Simulation{P<:Problem}
    problem  :: P
    workflow :: String
    lstps    :: Int
    folder   :: String
end
function Simulation(problem::Problem; workflow::AbstractString=RUN.workflow, lstps::Integer=RUN.lstps,
                    folder::AbstractString=RUN.folder)
    workflow ∈ WORKFLOWS || throw(ArgumentError("workflow \"$workflow\" is not one of $(WORKFLOWS)"))
    lstps ≥ 1 || throw(ArgumentError("lstps ≥ 1 is needed"))
    dir = rstrip(isabspath(folder) ? normpath(folder) : joinpath(self.sys.dump, folder), ['/', '\\'])
    return Simulation(problem, String(workflow), Int(lstps), String(dir))
end
