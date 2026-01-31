export set_config, export_setup

"""
    set_config(ic::NamedTuple, bc::NamedTuple, instr::Dict, tC::Real) -> NamedTuple

Generates a configuration `NamedTuple` that consolidates the necessary instructions and parameters to accurately execute the `DISPATCH!()` and `cORE!()/cORES!()` methods.

# Arguments
- `ic::NamedTuple`: The initial condition named tuple, containing data required to initialize the solver.
- `bc::NamedTuple`: The boundary condition named tuple, which provides the necessary boundary data for the solver.
- `instr::Dict`: A dictionary of instructions that guide the configuration process, ensuring that all required parameters are correctly set.
- `tC::Real`: A real number representing the current simulation time or a critical time parameter used in the configuration.

# Returns
- `NamedTuple`: A named tuple containing the compiled configuration, including instructions and parameters that are essential for the correct execution of the `DISPATCH!()` and `cORE!()/cORES!()` methods.

# Description
The `setConfig` function aggregates the initial conditions (`ic`), boundary conditions (`bc`), and additional instructions (`instr`) into a comprehensive configuration that is tailored for the execution of the `DISPATCH!()` and `cORE!()/cORES!()` methods. This configuration ensures that all necessary parameters are accounted for and properly initialized, facilitating a smooth and accurate simulation process.

- **Initial and Boundary Conditions**: The function combines the initial and boundary conditions, ensuring that both are fully integrated into the configuration. This step is crucial for initializing the solver with the correct state and boundary data.
- **Instruction Handling**: The provided `instr::Dict` is used to set up various operational parameters, ensuring that the solver behaves according to the specified guidelines. These instructions may include settings for computational options, solver-specific parameters, and other critical execution details.
- **Time Parameter**: The `tC` parameter is included to manage timing-related aspects of the simulation, such as the current simulation time or a key time-based configuration.


"""
function set_config(ic::NamedTuple,time::Vector,instr::Dict,tC::Any)
    T1,T2 = instr[:dtype][:T0]

    nproc = length(instr[:backend])
    mpi   = Vector{NamedTuple}(undef, nproc)
    glob  = get_glob((ic.dim[1],ic.dim[2]),nproc)
    for k ∈ 1:nproc
        mpi[k] = (; 
            status    = if nproc > 1 true else false end,
            root      = 0,
            me        = 0,
            info      = missing,
            comm      = missing,
            comm_size = 1,
            neighb    = glob[:neighb][k],
            index     = glob[:index][k],
            tplgy     = glob[:tplgy],
            globdim   = T1.((ic.dim[1],ic.dim[2])),
        )
    end


    return config::String
end

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
function export_setup(mesh::Mesh,mpts::Point,cmpr,time::Time,instr::Instruction,paths,misc; path::String = " ", file::String = "simulation_setup")
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
        cfg["instr"] = instr
        cfg["paths"] = paths
        cfg["misc" ] = misc
    end
    return sim
end
