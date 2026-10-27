"""
    get_node_type(ndim, nno) -> Matrix

The position of each node relative to the domain edges, per axis (`ndim × nno[end]`): `1` the
outermost layer, `2` the second, `4` the second-to-last, `3` interior. `BSplineBasis` uses it for
its boundary-modified splines.
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