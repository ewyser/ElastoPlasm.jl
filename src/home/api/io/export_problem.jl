export export_problem

"""
    export_problem(problem, basis, solver, paths; path, file) -> String

Export simulation setup to JLD2 file.

`problem` (bundling `mesh`/`mpts`/`time`) is persisted as a single `ic["problem"]`
entry rather than three separate `ic["mesh"]`/`ic["mpts"]`/`ic["time"]` keys — `Basis`
stays its own top-level `ic["basis"]` key, since it's deliberately independent of
`Problem` (see `MechanicalProblem`'s docstring).

# Arguments
- `problem::MechanicalProblem`: Bundles `mesh`, `mpts`, `time` (see `setup_problem`).
- `basis::Basis`: Topology container linking `problem.mesh` and `problem.mpts`.
- `solver`: Solver instance (e.g. `ExplicitSolver`)
- `paths`: Path configuration
- `path::String`: Directory path for output (default: " ")
- `file::String`: Filename without extension (default: "simulation_setup")

# Returns
- `String`: Path to created JLD2 file

# Example
```julia
sim_file = export_problem(problem, basis, solver, paths;
                        path=pwd(), file="my_simulation")
```
"""
function export_problem(problem::MechanicalProblem,basis::Basis,solver::S,paths; path::String = " ", file::String = "simulation_setup") where {T1<:Integer,T2<:Real,D, S<:AbstractSolver{T1,T2,D}}
    mesh,mpts,time = problem.mesh,problem.mpts,problem.time
    # display summary
    @info ic_log(mesh,mpts,time,solver)
    misc = (;
        prefix = "$(length(mesh.prprt.nel)-1)d_$(solver.transfer.trsfr)"
    )
    # create jld2 file
    sim = joinpath(path,"$(file).jld2")
    jldopen(sim, "w") do fid
        # create initial condition (ic) jld2 group
        ic = JLD2.Group(fid,"ic")
        ic["problem"] = problem
        ic["basis"]   = basis
        # create configuration (cfg) jld2 group
        cfg = JLD2.Group(fid,"cfg")
        cfg["solver"] = solver
        cfg["paths"]  = paths
        cfg["misc" ]  = misc
    end
    return sim
end
