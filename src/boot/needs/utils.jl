"""
    tree(sucess, prefix="\n\t", level=0, max_level=1) -> Vector{String}

Formats a list of strings into a tree-like structure for display.

# Arguments
- `sucess`: List of strings to format.
- `prefix="\n\t"`: String prefix for each line.
- `level=0`: Current tree depth.
- `max_level=1`: Maximum tree depth.

# Returns
- `Vector{String}`: Tree-formatted strings.

# Example
```julia
tree(["boot", "home"])
```
"""
function tree(sucess, prefix="\n\t", level=0, max_level=1)
    if level > max_level
        return nothing
    end
    n,printout = length(sucess),[]
    for (i, name) ∈ enumerate(sucess)
        connector = i == n ? "└── " : "├── "
        push!(printout,prefix*connector*name)
    end
    return printout
end

export clean_dump
"""
    clean_dump()

Empty `dump/`, the default folder of simulations, their plots and outputs. Nothing else empties it:
`setup_simulation` only overwrites the files of the simulation it writes (and removes its stale output).
"""
function clean_dump()
    dir   = self.sys.dump
    files = isdir(dir) ? readdir(dir; join=true) : String[]
    isempty(files) && return nothing
    foreach(f -> rm(f; recursive=true), files)
    @info join(["Emptied $(trunc_path(dir; anchor=basename(dir)))/:"; ["\e[31m- $(trunc_path(f))\e[0m" for f ∈ files]], "\n")
    return nothing
end

"""
    trunc_path(full_path::AbstractString; anchor::AbstractString="ElastoPlasm.jl") -> String

Returns the subpath of `full_path` starting from the directory name `anchor`.

# Arguments
- `full_path`: The full absolute or relative path.
- `anchor`: The folder name from which you want to keep the rest of the path.

# Returns
- `String`: Truncated path string.

# Example
```julia
trunc_path("C:/Users/lili8/Documents/GitHub/ElastoPlasm.jl/dump/slump", "ElastoPlasm.jl")
# => "ElastoPlasm.jl/dump/slump"
```
"""
function trunc_path(full_path::AbstractString; anchor::AbstractString="ElastoPlasm.jl")
    parts = splitpath(full_path)
    idx = findfirst(==(anchor), parts)
    return isnothing(idx) ? full_path : joinpath(parts[idx:end]...)
end

"""
    get_version() -> String

Return the current project version as a string, as specified in the Julia project file.

# Returns
- `String`: The project version.

# Example
```julia
v = get_version()
println(v)
```
"""
function get_version()
    return string(Pkg.project().version)
end