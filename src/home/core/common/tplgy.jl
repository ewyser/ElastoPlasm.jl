@inline function element_to_nodes_topology(xp::SVector{1,T2},x₀::SVector{1,T2},e2n::Vector{SVector{NN,T1}},h::SVector{1,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1]-x₀[1],h[1]) #clamp(fld(xp[1]-x₀[1],h[1]),0,nel[1]-1)
    el  = round(T1,1+elx)
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{2,T2},x₀::SVector{2,T2},e2n::Vector{SVector{NN,T1}},h::SVector{2,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1] - x₀[1], h[1]) #clamp(fld(xp[1] - x₀[1], h[1]), 0, nel[1] - 1)
    ely = fld(xp[2] - x₀[2], h[2]) #clamp(fld(xp[2] - x₀[2], h[2]), 0, nel[2] - 1)
    el  = round(T1, 1 + elx + nel[1] * ely)  # x varies fastest, (nx, ny) ordering
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{3,T2},x₀::SVector{3,T2},e2n::Vector{SVector{NN,T1}},h::SVector{3,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1] - x₀[1], h[1]) #clamp(fld(xp[1] - x₀[1], h[1]), 0, nel[1] - 1)
    ely = fld(xp[2] - x₀[2], h[2]) #clamp(fld(xp[2] - x₀[2], h[2]), 0, nel[2] - 1)
    elz = fld(xp[3] - x₀[3], h[3]) #clamp(fld(xp[3] - x₀[3], h[3]), 0, nel[3] - 1)
    el  = round(T1, 1 + elx + nel[1] * ely + nel[1] * nel[2] * elz)
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@kernel inbounds = true function p2e2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN}) where {T1,T2,D,NN}
    p = @index(Global)
    if p ≤ mpts.nmp
        el, e2n = element_to_nodes_topology(mpts.x[p], mesh.x₀, basis.e2n, mesh.prprt.h, mesh.prprt.nel)
        # Assign mpts-to-node connectivity
        basis.p2n[p] = e2n
        # Store element index in mpts
        basis.p2e[p] = el
    end
end

"""
    build_el2p(p2e::Vector{T1}, nel::T1) -> (ptr::Vector{T1}, idx::Vector{T1})

Counting-sort `p2e` (particle → its element, already populated by `p2e2n` every step)
into a CSR-style element→particle bucket list: particles in element `el` are
`idx[ptr[el]:ptr[el+1]-1]`. `O(nmp+nel)`, sequential (not a `@kernel`) — replaces the
old `Basis.e2p`/`p2p` dense `nmp×nel`/`nmp×nmp` matrices `nonlocal.jl` used to rebuild
and scan every timestep; see `nonlocal.jl`'s docstring for why this bucket + `e2e`
(genuinely sparse element neighbor structure) is enough to bound the neighbor search to
`O(nmp × k)` instead of `O(nmp²)`.
"""
function build_el2p(p2e::Vector{T1}, nel::T1) where {T1}
    nmp    = length(p2e)
    counts = zeros(T1, nel)
    for p ∈ 1:nmp
        counts[p2e[p]] += 1
    end
    ptr = Vector{T1}(undef, nel+1)
    ptr[1] = 1
    for el ∈ 1:nel
        ptr[el+1] = ptr[el] + counts[el]
    end
    cursor = copy(ptr[1:nel])
    idx    = Vector{T1}(undef, nmp)
    for p ∈ 1:nmp
        el         = p2e[p]
        idx[cursor[el]] = p
        cursor[el] += 1
    end
    return ptr, idx
end
