# Multiple dispatch for get_element_nodes based on neighbs type
function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, stencils::NTuple{1,UnitRange{T}}; neighbors::Vector{T}=T[]) where {T}
    nx = T(nno[1])
    i0 = el
    for i ∈ stencils[1]
        ii = i0 + i
        if 1 <= ii <= nx
            push!(neighbors, T(ii))
        else
            push!(neighbors, zero(T))
        end
    end
    return neighbors
end

function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, stencils::NTuple{2,UnitRange{T}}; neighbors::Vector{T}=T[]) where {T}
    nx, ny = T(nno[1]), T(nno[2])
    j0 = div(el-1, nel[1]) + 1
    i0 = el - nel[1]*(j0-1)
        for dj ∈ stencils[2]
            for di ∈ stencils[1]
                ii = i0 + di
                jj = j0 + dj
                if 1 <= ii <= nx && 1 <= jj <= ny
                    node = T(ii + nx*(jj-1))
                    push!(neighbors, node)
                else
                    push!(neighbors, zero(T))
                end
            end
        end
    return neighbors
end

function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, stencils::NTuple{3,UnitRange{T}}; neighbors::Vector{T}=T[]) where {T}
    nx, ny, nz = T(nno[1]), T(nno[2]), T(nno[3])
    k0 = div(el-1, nel[1]*nel[2]) + 1
    rem2 = el - nel[1]*nel[2]*(k0-1)
    j0 = div(rem2-1, nel[1]) + 1
    i0 = rem2 - nel[1]*(j0-1)
    for dk ∈ stencils[3]
        for dj ∈ stencils[2]
            for di ∈ stencils[1]
                ii = i0 + di
                jj = j0 + dj
                kk = k0 + dk
                if 1 <= ii <= nx && 1 <= jj <= ny && 1 <= kk <= nz
                    node = T(ii + nx*(jj-1) + nx*ny*(kk-1))
                    push!(neighbors, node)
                else
                    push!(neighbors, zero(T))
                end
            end
        end
    end
    return neighbors
end

"""
    get_element_to_nodes(nel, nno, stencils) -> Vector{SVector{NN}}

For each element of the structured grid (`nel`, `nno`: counts per axis, then the total), the
indices of the `NN` nodes its shape functions reach: the node offsets `stencils` (one range per
axis, e.g. `-1:2` for cubic B-splines) around the element, 0 where they fall outside the grid.
"""
function get_element_to_nodes(nel::Vector{T}, nno::Vector{T}, stencils::NTuple{D,UnitRange{T}}) where {T, D}
    NN  = prod(length.(stencils))
    e2n = Vector{SVector{NN,T}}(undef, nel[end])
    for el ∈ T.(collect(1:nel[end]))
        e2n[el] = SVector{NN,T}(get_nodes(el, nno, nel, stencils))
    end
    return e2n
end
