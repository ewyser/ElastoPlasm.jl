export set_paths, config_plot

"""
    set_paths(new_dir::String, path::String; interactive::Bool=true) -> Dict{Symbol, String}

Create and manage subdirectories within a specified base path for simulation output, with optional interactive selection of which folders to generate.

# Arguments
- `new_dir::String`: Name of the main directory to create within the base path.
- `path::String`: The base directory in which `new_dir` will be created.
- `interactive::Bool=true`: If true, prompts the user to select which subdirectories to create; if false, creates all default subdirectories automatically.

# Returns
- `Dict{Symbol, String}`: Dictionary mapping subdirectory names (as symbols) to their absolute paths.

# Behavior
- Ensures the main directory exists within the base path.
- Interactively or automatically creates subdirectories (e.g., `plot`, `dat`).
- Cleans up specific file types in newly created subdirectories if they already existed.
- Logs information about generated and cleaned paths.

# Example
```julia
paths = set_paths("results", pwd(); interactive=false)
println(paths[:plot])  # Absolute path to the plot subdirectory
```
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
- `nothing`

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