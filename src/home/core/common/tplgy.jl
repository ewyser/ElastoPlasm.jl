@inline function element_to_nodes_topology(xp::SVector{1,T2},x₀::SVector{1,T2},e2n::Vector{SVector{NN,T1}},h::SVector{1,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1]-x₀[1],h[1])#clamp(fld(xp[1]-x₀[1],h[1]),0,nel[1]-1)
    el  = round(T1,1+elx)
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{2,T2},x₀::SVector{2,T2},e2n::Vector{SVector{NN,T1}},h::SVector{2,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1] - x₀[1], h[1])#clamp(fld(xp[1] - x₀[1], h[1]), 0, nel[1] - 1)
    ely = fld(xp[2] - x₀[2], h[2])#clamp(fld(xp[2] - x₀[2], h[2]), 0, nel[2] - 1)
    el  = round(T1, 1 + elx + nel[1] * ely)  # x varies fastest, (nx, ny) ordering
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{3,T2},x₀::SVector{3,T2},e2n::Vector{SVector{NN,T1}},h::SVector{3,T2},nel::Vector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1] - x₀[1], h[1])#clamp(fld(xp[1] - x₀[1], h[1]), 0, nel[1] - 1)
    ely = fld(xp[2] - x₀[2], h[2])#clamp(fld(xp[2] - x₀[2], h[2]), 0, nel[2] - 1)
    elz = fld(xp[3] - x₀[3], h[3])#clamp(fld(xp[3] - x₀[3], h[3]), 0, nel[3] - 1)
    el  = round(T1, 1 + elx + nel[1] * ely + nel[1] * nel[2] * elz)
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@kernel inbounds = true function p2e2n(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN}) where {T1,T2,D,NN,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp
        el, e2n = element_to_nodes_topology(mpts.x[p], mesh.x₀, basis.e2n, mesh.prprt.h, mesh.prprt.nel)
        # Assign mpts-to-node connectivity
        basis.p2n[p] = e2n
        # Store element index in mpts
        basis.p2e[p] = el
    end
end