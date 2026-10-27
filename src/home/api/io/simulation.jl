export save_simulation, load_simulation

"""
    save_simulation(sim, file=<folder>/<name>_simulation.jld2) -> file
    load_simulation(file) -> Simulation

Write a `Simulation` to `file` (by default in its folder, `name` the folder's last part) and read it
back; the file holds one `"simulation"`. A loaded simulation's folder is the file's, so a run folder can
be moved or copied.
"""
function save_simulation(sim::Simulation, file::AbstractString=joinpath(mkpath(sim.folder), "$(basename(sim.folder))_simulation.jld2"))
    jldsave(file; simulation=sim)
    return file
end
function load_simulation(file::AbstractString)
    s = jldopen(file, "r") do f
        haskey(f, "simulation") || throw(ArgumentError("$file was written before Simulation: set it up again, " *
                                                       "setup_simulation(load_case(<name>_case.toml))"))
        f["simulation"]
    end
    return Simulation(s.problem, s.workflow, s.lstps, dirname(abspath(file)))
end
