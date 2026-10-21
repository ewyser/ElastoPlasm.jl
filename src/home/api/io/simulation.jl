export save_simulation, load_simulation

"""
    save_simulation(sim) -> file
    load_simulation(file) -> Simulation

Write a `Simulation` to `<folder>/<name>_simulation.jld2` (`name` the folder's last part) and read it
back; the file holds one `"simulation"`.
"""
function save_simulation(sim::Simulation)
    p = sim.problem
    @info ic_log(p.mesh, p.mpts, p.time, p.solver)
    file = joinpath(mkpath(sim.folder), "$(basename(sim.folder))_simulation.jld2")
    jldsave(file; simulation=sim)
    return file
end
function load_simulation(file::AbstractString)
    return jldopen(file, "r") do f
        haskey(f, "simulation") || throw(ArgumentError("$file was written before Simulation: set it up again, " *
                                                       "setup_simulation(load_case(<name>_case.toml))"))
        f["simulation"]
    end
end
