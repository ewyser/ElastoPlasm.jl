"""
    nonlocal_pq(W, mpts::Point{T1,T2}, basis::Basis{T1,T2}, ls::T2, el2p::Adjacency) where {T1,T2}
    nonlocal_qp(W, mpts::Point{T1,T2}, basis::Basis{T1,T2}, ls::T2, el2p::Adjacency) where {T1,T2}

Apply nonlocal averaging for regularization of plastic strain at material points.
`nonlocal_pq` accumulates each yielding particle's own normalization weight `W[p]`;
`nonlocal_qp` (run after `nonlocal_pq` has filled `W`) uses it to average neighboring
particles' plastic strain into `mpts.s.ϵpII[p][2]`. Two separate kernels rather than one
kernel branching on a `type::String` — this codebase's convention elsewhere is one
small kernel per concern (see `elast!`, or `ignite`/`mapsto`'s per-phase `Cairn` table
entries), and a `String` compare on every particle on every call was pure overhead once
this could just be two dispatch-table entries instead.

Bounds the neighbor search to `O(nmp × k)` (`k` = local particle density within `ls`,
resolution-independent) instead of scanning all `nmp` particles per yielding particle:
candidate neighbors of `p` are found via `basis.e2e` (genuinely sparse element-neighbor
structure, already bounded by `ls` at `e2e` construction time) intersected with the
particles of each element, `el2p`, built once per step by `element_particles`
(below) from `basis.p2e`. This replaces the old `Basis.e2p`/`p2p` dense
`nmp×nel`/`nmp×nmp` matrices and the freshly-`spzeros`-allocated, scalar-indexed `w`
weight matrix (an `O(nnz)`-per-write `SparseMatrixCSC` `setindex!` pattern) with nothing
persistent beyond `W` — pairwise weights are cheap enough (`sqrt`+`exp`) to recompute in
both passes rather than cache.

Each particle only ever reads its own neighbors and writes its own `W[p]`/`ϵpII[p]` —
never another particle's slot: by symmetry of `ω(d(p,q))`, `W[q]` is correctly computed
by `q`'s own pass whenever `q` itself yields (the only case `W[q]` is read), so there is
no need for `p`'s pass to write into `W[q]` as the old algorithm did. This keeps each
kernel embarrassingly parallel per-particle, no cross-particle writes or atomics needed.

# Arguments
- `W`: Vector of per-particle normalization weights.
- `mpts::Point{T1,T2}`: Material point data structure.
- `basis::Basis{T1,T2}`: Topology (`e2e`, `p2e`).
- `ls::T2`: Nonlocal length scale.
- `el2p`: the particles of each element (`element_particles`), an `Adjacency`.

# Returns
- Updates `W` (`nonlocal_pq`) or `mpts.s.ϵpII` (`nonlocal_qp`) in-place.
"""
@views @kernel inbounds = true function nonlocal_pq(W,mpts::Point{T1,T2},basis::Basis{T1,T2},ls::T2,el2p::Adjacency) where {T1,T2}
    p = @index(Global)

    if p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        acc = T2(0.0)
        el  = basis.p2e[p]
        for k ∈ nzrange(basis.e2e, el)
            el′ = rowvals(basis.e2e)[k]
            for i ∈ nzrange(el2p, el′)
                q  = rowvals(el2p)[i]
                d  = norm(mpts.x[p] - mpts.x[q])
                acc += d/ls*exp(-(d/ls)^2)
            end
        end
        W[p] = acc
    end
end

@views @kernel inbounds = true function nonlocal_qp(W,mpts::Point{T1,T2},basis::Basis{T1,T2},ls::T2,el2p::Adjacency) where {T1,T2}
    p = @index(Global)

    if p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        if !(W[p] > T2(1e-16))   # no neighbour weight: keep the local value
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], mpts.s.ϵpII[p][1])
        else
            acc = mpts.s.ϵpII[p][2]
            el  = basis.p2e[p]
            for k ∈ nzrange(basis.e2e, el)
                el′ = rowvals(basis.e2e)[k]
                for i ∈ nzrange(el2p, el′)
                    q  = rowvals(el2p)[i]
                    d  = norm(mpts.x[p] - mpts.x[q])
                    ω₀ = d/ls*exp(-(d/ls)^2)
                    acc += (ω₀/W[p])*mpts.s.ϵpII[q][1]
                end
            end
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], acc)
        end
    end
end

"""
    element_particles(mpts, mesh, basis, solver) -> Adjacency

The particles of each element (`p2e`, set by `p2e2n` every step, `tplgy.jl`) as an `Adjacency`: those in element `el` are
`rowvals(a)[nzrange(a, el)]`, in ascending order. Built by kernels where the particles live: count the particles per
element (`el2p_count`, atomic), sum the counts on the host (`O(nel)`, one copy each way), put each particle in its
element's next free slot (`el2p_fill`, atomic), then sort each element's few indices (`el2p_sort`), so the order is
the same on every device and thread count, that of the sequential counting sort this replaced, and the nonlocal sums
over it stay deterministic. With `basis.e2e` it bounds the neighbour search of `nonlocal_pq`/`nonlocal_qp` (above)
to `O(nmp × k)` instead of `O(nmp²)` (see their docstring).
"""
function element_particles(mpts::Point{T1},mesh::Mesh{T1},basis::Basis{T1},solver::AbstractSolver) where {T1}
    bk     = exec_backend(solver)
    nel    = Int(mesh.prprt.nel[end])
    counts = fill!(similar(basis.p2e, nel), zero(T1))
    solver.cairn.update.el2p_count!(counts, basis; ndrange=mpts.nmp)
    sync(bk)
    starts = vcat(one(T1), one(T1) .+ cumsum(Array(counts)))          # ptr: each element's first slot
    ptr    = copyto!(similar(basis.p2e, nel + 1), starts)
    cursor = copyto!(similar(basis.p2e, nel), starts[1:nel])
    idx    = similar(basis.p2e, Int(mpts.nmp))
    solver.cairn.update.el2p_fill!(idx, cursor, basis; ndrange=mpts.nmp)
    sync(bk)
    solver.cairn.update.el2p_sort!(idx, ptr; ndrange=nel)
    sync(bk)
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
