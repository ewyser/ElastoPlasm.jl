

"""
    get_vals(mesh, mpts, it) -> Vector{Tuple{String, Any}}

Return a summary of mesh and simulation state for progress display or logging.

# Arguments
- `mesh`: Mesh data structure, must contain `nel`.
- `mpts`: Material point data structure, must contain `nmp`.
- `it`: Current iteration or step count.

# Returns
- `Vector{Tuple{String, Any}}`: List of (label, value) pairs summarizing the state.

# Example
```julia
vals = get_vals(mesh, mpts, 10)
println(vals)
```
"""
function get_vals(mesh,mpts,it)
    # save vals
    vals = [
            ("nel,np",(round(Int64,prod(mesh.nel[1:end-1])),mpts.nmp)),
            ("iteration(s)",it),
            ]
    return vals
end