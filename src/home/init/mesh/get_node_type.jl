"""
    get_coords(ndim::T1, L::Vector{T2}, h::Vector{T2}; ghosts::Vector{T2}=[T2(0.0)]) where {T1,T2}

Generate nodal coordinates and mesh topology for 1D, 2D, or 3D domains.

# Arguments
- `ndim::T1`: Number of spatial dimensions.
- `L::Vector{T2}`: Length of the domain in each spatial direction.
- `h::Vector{T2}`: Element size in each direction.
- `ghosts::Vector{T2}`: Optional, size of ghost cells to add at boundaries (default: `[T2(0.0)]`).

# Returns
- `x::Matrix{T2}`: Matrix of nodal coordinates (each column is a node).
- `nel::Vector{Int}`: Number of elements in each direction and total.
- `nno::Vector{Int}`: Number of nodes in each direction and total.

# Example
```julia
x, nel, nno = get_coords(2, [1.0, 1.0], [0.1, 0.1])
```
"""
function get_node_type(ndim::T1,nno::Vector{T1}) where {T1}
    if ndim == 1
        node_type = fill(T1(3), nno[1])
    elseif ndim == 2
        nx, ny = nno[1], nno[2]
        type = fill(T1(3), nx, ny)
        type[1    , :] .= T1(1)
        type[2    , :] .= T1(2)
        type[end-1, :] .= T1(4)
        type[end  , :] .= T1(1)
        node_type = vec(type)'

        type = fill(T1(3), nx, ny)
        type[:, 1    ] .= T1(1)
        type[:, 2    ] .= T1(2)
        type[:, end-1] .= T1(4)
        type[:, end  ] .= T1(1)
        node_type = vcat(node_type, vec(type)')
    elseif ndim == 3
        nx, ny, nz = nno[1], nno[2], nno[3]
        type = fill(T1(3), nx, ny, nz)
        type[1    , :, : ] .= T1(1)
        type[2    , :, : ] .= T1(2)
        type[end-1, :, : ] .= T1(4)
        type[end  , :, : ] .= T1(1)
        node_type = vec(type)'

        type = fill(T1(3), nx, ny, nz)
        type[:, 1    , : ] .= T1(1)
        type[:, 2    , : ] .= T1(2)
        type[:, end-1, : ] .= T1(4)
        type[:, end  , : ] .= T1(1)
        node_type = vcat(node_type, vec(type)')

        type = fill(T1(3), nx, ny, nz)
        type[:, :, 1    ] .= T1(1)
        type[:, :, 2    ] .= T1(2)
        type[:, :, end-1] .= T1(4)
        type[:, :, end  ] .= T1(1)
        node_type = vcat(node_type, vec(type)')
    end

    return node_type
end