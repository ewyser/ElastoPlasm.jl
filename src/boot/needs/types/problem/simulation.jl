export Simulation

"""
    Simulation(problem; folder="problem")

A problem ready to run, and the `folder` its files go to (a name under `dump/`, or a path; `setup_simulation` gives it
the case's): its workflow and load steps are its solver's run settings (`problem.solver.run`). A simulation's `.jld2` holds one
(`save_simulation`/`load_simulation`); `execute!` runs it.
"""
struct Simulation{P<:Problem}
    problem ::P
    folder  ::String
end
function Simulation(problem::Problem; folder::AbstractString=RUN.folder)
    dir = rstrip(isabspath(folder) ? normpath(folder) : joinpath(self.sys.dump, folder), ['/', '\\'])
    return Simulation(problem, String(dir))
end
