#= =#
"""
    nonlocal(W, mpts::Point{T1,T2}, basis::Basis{T1,T2}, ls::T2, ptr::Vector{T1}, idx::Vector{T1}, type::String) where {T1,T2}

Apply nonlocal averaging for regularization of plastic strain at material points.

Bounds the neighbor search to `O(nmp × k)` (`k` = local particle density within `ls`,
resolution-independent) instead of scanning all `nmp` particles per yielding particle:
candidate neighbors of `p` are found via `basis.e2e` (genuinely sparse element-neighbor
structure, already bounded by `ls` at `e2e` construction time) intersected with the
`(ptr,idx)` CSR element→particle bucket list built once per step by `build_el2p`
(`tplgy.jl`) from `basis.p2e`. This replaces the old `Basis.e2p`/`p2p` dense
`nmp×nel`/`nmp×nmp` matrices and the freshly-`spzeros`-allocated, scalar-indexed `w`
weight matrix (an `O(nnz)`-per-write `SparseMatrixCSC` `setindex!` pattern) with nothing
persistent beyond `W` — pairwise weights are cheap enough (`sqrt`+`exp`) to recompute in
both passes rather than cache.

Each particle only ever reads its own neighbors and writes its own `W[p]`/`ϵpII[p]` —
never another particle's slot: by symmetry of `ω(d(p,q))`, `W[q]` is correctly computed
by `q`'s own pass whenever `q` itself yields (the only case `W[q]` is read), so there is
no need for `p`'s pass to write into `W[q]` as the old algorithm did. This keeps the
kernel embarrassingly parallel per-particle, no cross-particle writes or atomics needed.

# Arguments
- `W`: Vector of per-particle normalization weights.
- `mpts::Point{T1,T2}`: Material point data structure.
- `basis::Basis{T1,T2}`: Topology (`e2e`, `p2e`).
- `ls::T2`: Nonlocal length scale.
- `ptr`,`idx`: CSR element→particle bucket list from `build_el2p`.
- `type::String`: Operation type ("p->q" or "p<-q").

# Returns
- Updates `W` and `mpts.s.ϵpII` in-place.
"""
@views @kernel inbounds = true function nonlocal(W,mpts::Point{T1,T2},basis::Basis{T1,T2},ls::T2,ptr::Vector{T1},idx::Vector{T1},type::String) where {T1,T2}
    p = @index(Global)

    if type == "p->q" && p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        acc = T2(0.0)
        el  = basis.p2e[p]
        for k ∈ nzrange(basis.e2e, el)
            el′ = rowvals(basis.e2e)[k]
            for i ∈ ptr[el′]:ptr[el′+1]-T1(1)
                q  = idx[i]
                d  = norm(mpts.x[p] - mpts.x[q])
                acc += d/ls*exp(-(d/ls)^2)
            end
        end
        W[p] = acc
    elseif type == "p<-q" && p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        if isapprox(W[p]>T2(1e-16),T2(0.0),atol=T2(1e-16))
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], mpts.s.ϵpII[p][1])
        else
            acc = mpts.s.ϵpII[p][2]
            el  = basis.p2e[p]
            for k ∈ nzrange(basis.e2e, el)
                el′ = rowvals(basis.e2e)[k]
                for i ∈ ptr[el′]:ptr[el′+1]-T1(1)
                    q  = idx[i]
                    d  = norm(mpts.x[p] - mpts.x[q])
                    ω₀ = d/ls*exp(-(d/ls)^2)
                    acc += (ω₀/W[p])*mpts.s.ϵpII[q][1]
                end
            end
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], acc)
        end
    end
end
