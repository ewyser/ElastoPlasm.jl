@inline function element_to_nodes_topology(xp::SVector{1,T2},x₀::SVector{1,T2},e2n::AbstractVector{SVector{NN,T1}},h::SVector{1,T2},nel::AbstractVector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1]-x₀[1],h[1]) #clamp(fld(xp[1]-x₀[1],h[1]),0,nel[1]-1)
    el  = round(T1,1+elx)
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{2,T2},x₀::SVector{2,T2},e2n::AbstractVector{SVector{NN,T1}},h::SVector{2,T2},nel::AbstractVector{T1}) where {T1,T2,NN}
    # Compute element indices
    elx = fld(xp[1] - x₀[1], h[1]) #clamp(fld(xp[1] - x₀[1], h[1]), 0, nel[1] - 1)
    ely = fld(xp[2] - x₀[2], h[2]) #clamp(fld(xp[2] - x₀[2], h[2]), 0, nel[2] - 1)
    el  = round(T1, 1 + elx + nel[1] * ely)  # x varies fastest, (nx, ny) ordering
    # Assign mp-to-no connectivity
    p2n = e2n[el]
    return el,p2n
end
@inline function element_to_nodes_topology(xp::SVector{3,T2},x₀::SVector{3,T2},e2n::AbstractVector{SVector{NN,T1}},h::SVector{3,T2},nel::AbstractVector{T1}) where {T1,T2,NN}
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
    element_particles(mpts, mesh, basis, solver) -> Adjacency

The particles of each element (`p2e`, set by `p2e2n` every step) as an `Adjacency`: those in element `el` are
`rowvals(a)[nzrange(a, el)]`, in ascending order. Built by kernels where the particles live: count the particles per
element (`el2p_count`, atomic), sum the counts on the host (`O(nel)`, one copy each way), put each particle in its
element's next free slot (`el2p_fill`, atomic), then sort each element's few indices (`el2p_sort`), so the order is
the same on every device and thread count, that of the sequential counting sort this replaced, and the nonlocal sums
over it stay deterministic. With `basis.e2e` it bounds `nonlocal.jl`'s neighbour search to `O(nmp × k)` instead of
`O(nmp²)` (see its docstring).
"""
function element_particles(mpts::Point{T1},mesh::Mesh{T1},basis::Basis{T1},solver::AbstractSolver) where {T1}
    bk     = get_backend(mpts.x)
    nel    = Int(mesh.prprt.nel[end])
    counts = fill!(similar(basis.p2e, nel), zero(T1))
    solver.cairn.update.el2p_count!(counts, basis; ndrange=mpts.nmp);sync(bk)
    starts = vcat(one(T1), one(T1) .+ cumsum(Array(counts)))          # ptr: each element's first slot
    ptr    = copyto!(similar(basis.p2e, nel + 1), starts)
    cursor = copyto!(similar(basis.p2e, nel), starts[1:nel])
    idx    = similar(basis.p2e, Int(mpts.nmp))
    solver.cairn.update.el2p_fill!(idx, cursor, basis; ndrange=mpts.nmp);sync(bk)
    solver.cairn.update.el2p_sort!(idx, ptr; ndrange=nel);sync(bk)
    return Adjacency(ptr, idx)
end
@kernel inbounds = true function el2p_count(counts, basis::Basis{T1}) where {T1}
    p = @index(Global)
    if p ≤ length(basis.p2e)
        @atom counts[basis.p2e[p]] += one(T1)
    end
end
@kernel inbounds = true function el2p_fill(idx, cursor, basis::Basis{T1}) where {T1}
    p = @index(Global)
    if p ≤ length(basis.p2e)
        i      = (@atom cursor[basis.p2e[p]] += one(T1)) - one(T1)   # the slot taken (@atom returns the new value)
        idx[i] = T1(p)
    end
end
# insertion sort of an element's indices: a handful per element
@kernel inbounds = true function el2p_sort(idx, ptr)
    el = @index(Global)
    if el < length(ptr)
        lo, hi = ptr[el], ptr[el+1] - one(eltype(ptr))
        for i ∈ lo+one(lo):hi
            q, j = idx[i], i - one(i)
            while j ≥ lo && idx[j] > q
                idx[j+one(j)] = idx[j]
                j -= one(j)
            end
            idx[j+one(j)] = q
        end
    end
end
