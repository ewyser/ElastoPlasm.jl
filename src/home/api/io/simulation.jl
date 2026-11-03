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
    jldsave(file; simulation=Simulation(with_solver(p, stored_solver(p.solver)), sim.folder))
    return file
end
function load_simulation(file::AbstractString)
    s = jldopen(file, "r") do f
        haskey(f, "simulation") || throw(ArgumentError("$file was written before Simulation: set it up again, " *
                                                       "setup(load_case(<name>_case.toml))"))
        f["simulation"]
    end
    s isa Simulation || throw(ArgumentError("$file was written by an earlier ElastoPlasm whose types differ: set it up " *
                                            "again, setup(load_case(<name>_case.toml))"))
    p      = s.problem
    solver = build_solver(solver_config(p.solver), length(p.mesh.x₀))   # its kernels and device, on this machine
    return Simulation(with_solver(p, solver), dirname(abspath(file)))
end
# a built problem with another solver of the same configuration: the one stored, or the one rebuilt on loading
function with_solver(p::Problem{T1,T2,D,P,M,B}, solver::Solver{T1,T2,D}) where {T1,T2,D,P,M,B}
    return Problem{T1,T2,D,P,M,B,typeof(solver)}(p.mpts, p.mesh, p.basis, p.time, solver)
end
