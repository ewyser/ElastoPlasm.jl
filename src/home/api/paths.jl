export set_paths

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
