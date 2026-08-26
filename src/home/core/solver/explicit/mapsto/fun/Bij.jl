"""
    Bij(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, basis::Basis{...,TR}) where {TR<:ApicTransfer}

Update the per-particle affine velocity matrix `basis.transfer.Bᵢⱼ` (APIC scheme,
n-dim). `TR<:StdTransfer`/`TpicTransfer` get no-op methods below — both schemes have no
`Bᵢⱼ` to update, and dispatching a no-op lets `mapsto.jl` call this kernel
unconditionally instead of behind an `if transfer.trsfr=="apic"` check.

# Arguments
- `mpts::Point{T1,T2,D}`: Material point data structure.
- `mesh::Mesh{T1,T2,D}`: Mesh data structure.
- `basis::Basis{...}`: Topology container; `basis.transfer.Bᵢⱼ` is written in place.

# Returns
- Updates `basis.transfer.Bᵢⱼ` in-place.
"""
@kernel inbounds = true function Bij(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN,K,TR}) where {T1,T2,D,NN,K,TR<:ApicTransfer}
    p = @index(Global)
    if p ≤ mpts.nmp
        # Bᵢⱼ update
        Bᵢⱼ,xₚ = zero(eltype(basis.transfer.Bᵢⱼ)), mpts.x[p]
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N    = basis.N[nn,p]
            δ    = mesh.x[no] - xₚ
            Bᵢⱼ  = Bᵢⱼ + N * (mesh.s.v[no] * δ')
        end
        basis.transfer.Bᵢⱼ[p] = Bᵢⱼ
    end
end
@kernel inbounds = true function Bij(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN,K,TR}) where {T1,T2,D,NN,K,TR<:Union{StdTransfer,TpicTransfer}}
    p = @index(Global)
    # no-op: std/tpic have no Bᵢⱼ to update
end
