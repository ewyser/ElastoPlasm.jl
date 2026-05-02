# Multiple dispatch for get_element_nodes based on neighbs type
function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, nbs::neighbs{T,D}; neighbors::Vector{T}=T[]) where {T<:Integer,D<:OneDimension}
    nx = T(nno[1])
    i0 = el
    for i ∈ nbs.stencils[1]
        ii = i0 + i
        if 1 <= ii <= nx
            push!(neighbors, T(ii))
        else
            push!(neighbors, zero(T))
        end
    end
    return neighbors
end

function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, nbs::neighbs{T,D}; neighbors::Vector{T}=T[]) where {T<:Integer,D<:TwoDimension}
    nx, ny = T(nno[1]), T(nno[2])
    j0 = div(el-1, nel[1]) + 1
    i0 = el - nel[1]*(j0-1)
    for dj ∈ nbs.stencils[2]
        for di ∈ nbs.stencils[1]
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

function get_nodes(el::T, nno::Vector{T}, nel::Vector{T}, nbs::neighbs{T,D}; neighbors::Vector{T}=T[]) where {T<:Integer,D<:ThreeDimension}
    nx, ny, nz = T(nno[1]), T(nno[2]), T(nno[3])
    k0 = div(el-1, nel[1]*nel[2]) + 1
    rem2 = el - nel[1]*nel[2]*(k0-1)
    j0 = div(rem2-1, nel[1]) + 1
    i0 = rem2 - nel[1]*(j0-1)
    for dk ∈ nbs.stencils[3]
        for dj ∈ nbs.stencils[2]
            for di ∈ nbs.stencils[1]
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
    e2n(nel::Vector{T}, nno::Vector{T}, nbs::neighbs{T,D}) where {T<:Integer, D<:AbstractDimension}

Construct the element-to-node connectivity matrix for a structured mesh using a generic neighbor stencil.

# Arguments
- `nel::Vector{T}`: Number of elements in each direction and total.
- `nno::Vector{T}`: Number of nodes in each direction and total.
- `nbs::neighbs{T,D}`: Neighbor stencil object (see `neighbs` struct).

# Returns
- `e2n::Matrix{T}`: Element-to-node connectivity matrix (each column contains the node indices for one element, padded with zeros for out-of-bounds neighbors).

# Example
```julia
# 2D mesh: 10x10 elements, 11x11 nodes
nno = [11, 11, 121]
nel = [10, 10, 100]
stencil = neighbs((Int(-1):Int(2), Int(-1):Int(2)))  # 4x4 stencil
e2n_mat = get_element_to_nodes(nel, nno, stencil)
```
"""
function get_element_to_nodes(nel::Vector{T}, nno::Vector{T}, nbs::neighbs{T,D}) where {T<:Integer, D<:AbstractDimension}
    e2n = zeros(T, nbs.nn, nel[end]) 
    for el ∈ T.(collect(1:nel[end]))
        println
        e2n[:, el] .= get_nodes(el, nno, nel, nbs)
    end
    return e2n
end
