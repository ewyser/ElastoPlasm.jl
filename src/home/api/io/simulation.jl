export save_simulation, load_simulation

"""
    save_simulation(sim, file=<folder>/<name>_simulation.jld2) -> file
    load_simulation(file) -> Simulation

Write a `Simulation` to `file` (by default in its folder, `name` the folder's last part) and read it
back; the file holds one `"simulation"`. A loaded simulation's folder is the file's, so a run folder can
be moved or copied. The solver is stored as its configuration (`stored_solver`: no kernels, no devices) and rebuilt
on loading (`build_solver`): the file opens on any machine, which picks its own device.
"""
function save_simulation(sim::Simulation, file::AbstractString=joinpath(mkpath(sim.folder), "$(basename(sim.folder))_simulation.jld2"))
    p = sim.problem
    jldsave(file; simulation=Simulation(Problem(p.mpts, p.mesh, p.basis, p.time, stored_solver(p.solver)), sim.workflow,
                                        sim.lstps, sim.folder))
    return file
end
function load_simulation(file::AbstractString)
    s = jldopen(file, "r") do f
        haskey(f, "simulation") || throw(ArgumentError("$file was written before Simulation: set it up again, " *
                                                       "setup_simulation(load_case(<name>_case.toml))"))
        f["simulation"]
    end
    s isa Simulation || throw(ArgumentError("$file was written by an earlier ElastoPlasm whose types differ: set it up " *
                                            "again, setup_simulation(load_case(<name>_case.toml))"))
    p      = s.problem
    solver = build_solver(solver_config(p.solver), length(p.mesh.x₀))   # its kernels and device, on this machine
    return Simulation(Problem(p.mpts, p.mesh, p.basis, p.time, solver), s.workflow, s.lstps, dirname(abspath(file)))
end
