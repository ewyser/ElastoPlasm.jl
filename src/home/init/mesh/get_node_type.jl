"""
    get_node_type(ndim::T1, nno::Vector{T1}) where {T1}

Classify mesh nodes by boundary layer: interior nodes vs. the two layers of nodes nearest each
domain edge, along every spatial axis.

# Arguments
- `ndim::T1`: Number of spatial dimensions.
- `nno::Vector{T1}`: Number of nodes in each direction and total.

# Returns
- `node_type::Matrix{T1}`: `ndim × nno[end]` matrix; each row classifies nodes along one axis as
  `1` (outermost layer), `2` (second layer), `4` (second-to-last layer), or `3` (interior).

# Example
```julia
node_type = get_node_type(2, [11, 11, 121])
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