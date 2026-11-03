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

Return ElastoPlasm's own version (its Project.toml's, whichever project is active) as a string.

# Returns
- `String`: The project version.

# Example
```julia
v = get_version()
println(v)
```
"""
function get_version()
    return string(pkgversion(@__MODULE__))
end

"""
    adapt_fields(to, x, P) -> P(...)

`x` with each field adapted to `to` (Adapt.jl: a KernelAbstractions backend moves arrays to its device, `Array`
brings them back), rebuilt by `P`: the container's partial type (`Point`, `Mesh{T1,T2,D}`, `Basis{T1,T2,D,NN}`, …),
which pins the type parameters no field carries and takes the arrays' types from the adapted fields. The adapt rules
of the mesh, point and basis containers call it. A plain-bits field (a scalar, `nothing`, an `SVector` such as
`mesh.x₀`, `MeshProperties`) is kept as it is: Adapt would turn an `SVector` into an `Array`, or a device array.
"""
@inline function adapt_fields(to, x, ::Type{P}) where {P}
    return P(ntuple(i -> adapt_field(to, getfield(x, i)), Val(fieldcount(typeof(x))))...)
end
@inline function adapt_field(to, f)
    return isbits(f) ? f : Adapt.adapt(to, f)
end

# what an option's value stands for, in a table keyed by the option's values (`get_constitutive_laws()`'s);
# `option` names it in the error
function table_value(table::NamedTuple, name::AbstractString, option::AbstractString)
    haskey(table, Symbol(name)) || throw(ArgumentError("$option = \"$name\" is not one of $(join(keys(table), ", "))"))
    return getfield(table, Symbol(name))
end
