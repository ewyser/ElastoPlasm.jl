export export_setup

"""
    export_setup(mesh, mpts, cmpr, time, instr, paths, misc; path, file) -> String

Export simulation setup to JLD2 file.

# Arguments
- `mesh::Mesh`: Mesh data structure
- `mpts::Point`: Material point data structure
- `cmpr`: Compression data
- `time::Time`: Time data structure
- `instr::Instruction`: Instruction configuration
- `paths`: Path configuration
- `misc`: Miscellaneous data
- `path::String`: Directory path for output (default: " ")
- `file::String`: Filename without extension (default: "simulation_setup")

# Returns
- `String`: Path to created JLD2 file

# Example
```julia
sim_file = export_setup(mesh, mpts, cmpr, time, instr, paths, misc; 
                        path=pwd(), file="my_simulation")
```
"""
function export_setup(mesh::Mesh,mpts::Point,cmpr,time::Time,solver::S,paths; path::String = " ", file::String = "simulation_setup") where {T1<:Integer,T2<:Real,D<:AbstractDimension, S<:AbstractSolver{T1,T2,D}}
    # display summary
    @info ic_log(mesh,mpts,time,solver)
    misc = (;
        prefix = "$(mesh.prprt.dim)d_$(solver.fwrk.trsfr)"
    )
    # create jld2 file
    sim = joinpath(path,"$(file).jld2")
    jldopen(sim, "w") do fid
        # create initial condition (ic) jld2 group
        ic = JLD2.Group(fid,"ic")
        ic["mesh"] = mesh
        ic["mpts"] = mpts
        ic["cmpr"] = cmpr
        ic["time"] = time
        # create configuration (cfg) jld2 group
        cfg = JLD2.Group(fid,"cfg")
        cfg["solver"] = solver
        cfg["paths"]  = paths
        cfg["misc" ]  = misc
    end
    return sim
end
