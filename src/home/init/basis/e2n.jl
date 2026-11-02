# Multiple dispatch for get_element_nodes based on neighbs type
function get_nodes(el::T1, nno::Vector{T1}, nel::Vector{T1}, stencils::NTuple{1,UnitRange{T1}}; neighbors::Vector{T1}=T1[]) where {T1}
    nx = T1(nno[1])
    i0 = el
    for i ∈ stencils[1]
        ii = i0 + i
        if 1 <= ii <= nx
            push!(neighbors, T1(ii))
        else
            push!(neighbors, zero(T1))
        end
    end
    return neighbors
end

function get_nodes(el::T1, nno::Vector{T1}, nel::Vector{T1}, stencils::NTuple{2,UnitRange{T1}}; neighbors::Vector{T1}=T1[]) where {T1}
    nx, ny = T1(nno[1]), T1(nno[2])
    j0 = div(el-1, nel[1]) + 1
    i0 = el - nel[1]*(j0-1)
        for dj ∈ stencils[2]
            for di ∈ stencils[1]
                ii = i0 + di
                jj = j0 + dj
                if 1 <= ii <= nx && 1 <= jj <= ny
                    node = T1(ii + nx*(jj-1))
                    push!(neighbors, node)
                else
                    push!(neighbors, zero(T1))
                end
            end
        end
    return neighbors
end

function get_nodes(el::T1, nno::Vector{T1}, nel::Vector{T1}, stencils::NTuple{3,UnitRange{T1}}; neighbors::Vector{T1}=T1[]) where {T1}
    nx, ny, nz = T1(nno[1]), T1(nno[2]), T1(nno[3])
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
                    node = T1(ii + nx*(jj-1) + nx*ny*(kk-1))
                    push!(neighbors, node)
                else
                    push!(neighbors, zero(T1))
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
function get_element_to_nodes(nel::Vector{T1}, nno::Vector{T1}, stencils::NTuple{D,UnitRange{T1}}) where {T1,D}
    NN  = prod(length.(stencils))
    e2n = Vector{SVector{NN,T1}}(undef, nel[end])
    for el ∈ T1.(collect(1:nel[end]))
        e2n[el] = SVector{NN,T1}(get_nodes(el, nno, nel, stencils))
    end
    return e2n
end
