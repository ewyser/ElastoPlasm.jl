export set_paths,set_config,set_config!,config_plot,export_setup

"""
    set_paths(new_dir::String, path::String; interactive::Bool=true) -> Dict{Symbol, String}

Create and manage subdirectories within a specified base path for simulation output, with optional interactive selection of which folders to generate.

# Arguments
- `new_dir::String`: Name of the main directory to create within the base path.
- `path::String`: The base directory in which `new_dir` will be created.
- `interactive::Bool=true`: If true, prompts the user to select which subdirectories to create; if false, creates all default subdirectories automatically.

# Returns
- `Dict{Symbol, String}`: Mapping subdirectory names (as symbols) to their absolute paths.

# Example
```julia
paths = set_paths("results", pwd(); interactive=false)
println(paths[:plot])  # Absolute path to the plot subdirectory
```

# Notes
- Ensures the main directory exists within the base path.
- Interactively or automatically creates subdirectories (e.g., `plot`, `dat`).
- Cleans up specific file types in newly created subdirectories if they already existed.
- Logs information about generated and cleaned paths.
"""
function set_paths(new_dir::String, path::String; interactive::Bool=true)
    paths    = Dict{Symbol,Any}()
    msg,msg! = ["Generating additional paths:"], ["Deleting at:"]
    options  = ["plot", "dat"]
    root     = joinpath(path, new_dir)
    if interactive
        select = request("Select folder(s) you'd like to generate:", MultiSelectMenu(options))
        for (k, name) ∈ enumerate(select)
            subdir = joinpath(root, options[name])
            paths[Symbol(options[name])] = subdir
            if !isdir(subdir)
                mkpath(subdir)
                push!(msg, "\n\e[32m+ $(trunc_path(paths[Symbol(options[name])]))\e[0m")
            end
        end 
    else
        # point all options to the root (but don't create subfolders)
        for name ∈ options
            paths[Symbol(name)] = root
            if !isdir(paths[Symbol(name)])
                mkpath(paths[Symbol(name)])
                push!(msg, "\n\e[32m+ $(trunc_path(paths[Symbol(name)]))\e[0m")
            end
        end
    end

    if length(msg) > 1
        @info join(msg)
    end
    if length(msg!) > 1
        @warn join(msg!)
    end

    return paths
end

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


export config_plot
"""
    config_plot(; titf=12, gf=12, tickf=10, lf=10, lw=2, fs=:box, l=nothing, g=false)

Configure the default styling parameters for plots, allowing customization of fonts, line widths, frame style, labels, and grid visibility.

# Keyword Arguments
- `titf`: Font size for the plot title (default: `12`).
- `gf`: Font size for guide elements, such as axis labels and legends (default: `12`).
- `tickf`: Font size for tick labels (default: `10`).
- `lf`: Font size for legend text (default: `10`).
- `lw`: Line width for plot lines (default: `2`).
- `fs`: Frame style for the plot, e.g., `:box` or `:none` (default: `:box`).
- `l`: Label for the plot, or `nothing` for no label (default: `nothing`).
- `g`: Grid display option; `true` to show grid, `false` to hide (default: `false`).

# Returns
- `Nothing`. Configures plot defaults.

# Example
```julia
config_plot(titf=14, gf=12, tickf=8, lf=12, lw=3, fs=:box, l="My Plot", g=true)
```
"""
function config_plot(; titf=12, gf=12, tickf=10, lf=10, lw=2, fs=:box, l=nothing, g=false)
    default(
        fontfamily  = "Computer Modern",
        titlefont   = titf, 
        guidefont   = gf,  
        tickfont    = tickf, 
        legendfont  = lf,
        linewidth   = lw,
        framestyle  = fs,
        label       = l,
        grid        = g,
    )
    return nothing
end

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