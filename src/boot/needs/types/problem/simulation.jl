export Simulation

"""
    Simulation(problem; folder="problem")

A problem ready to run, and the folder its files go to: its workflow and load steps are its solver's run settings
(`problem.solver.run`). A simulation's `.jld2` holds one (`save_simulation`/`load_simulation`); `execute!` runs it.

# Fields
- `problem::P`: the `Problem`, which `execute!` runs in place.
- `folder::String`: where its files go (the `.jld2`, the case's TOML, plots, the HDF5 output), an absolute path: a name
  given is under `dump/`. `setup_simulation` gives it the case's (`case_folder`); a loaded simulation's is its file's.

`P` is the problem's type.
"""
struct Simulation{P}
    problem ::P
    folder  ::String
end
function Simulation(problem::Problem; folder::AbstractString=get_run().folder)
    dir = rstrip(isabspath(folder) ? normpath(folder) : joinpath(self.sys.dump, folder), ['/', '\\'])
    return Simulation(problem, String(dir))
end
