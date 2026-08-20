export export_setup

"""
    export_setup(mesh, mpts, basis, cmpr, time, solver, paths; path, file) -> String

Export simulation setup to JLD2 file.

# Arguments
- `mesh::Mesh`: Mesh data structure
- `mpts::Point`: Material point data structure
- `basis::Basis`: Topology container linking `mesh` and `mpts`
- `cmpr`: Compression data
- `time::Time`: Time data structure
- `solver`: Solver instance (e.g. `ExplicitSolver`)
- `paths`: Path configuration
- `path::String`: Directory path for output (default: " ")
- `file::String`: Filename without extension (default: "simulation_setup")

# Returns
- `String`: Path to created JLD2 file

# Example
```julia
sim_file = export_setup(mesh, mpts, basis, cmpr, time, solver, paths;
                        path=pwd(), file="my_simulation")
```
"""
function export_setup(mesh::Mesh,mpts::Point,basis::Basis,cmpr,time::Time,solver::S,paths; path::String = " ", file::String = "simulation_setup") where {T1<:Integer,T2<:Real,D, S<:AbstractSolver{T1,T2,D}}
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
        ic["mesh"]  = mesh
        ic["mpts"]  = mpts
        ic["basis"] = basis
        ic["cmpr"]  = cmpr
        ic["time"]  = time
        # create configuration (cfg) jld2 group
        cfg = JLD2.Group(fid,"cfg")
        cfg["solver"] = solver
        cfg["paths"]  = paths
        cfg["misc" ]  = misc
    end
    return sim
end
