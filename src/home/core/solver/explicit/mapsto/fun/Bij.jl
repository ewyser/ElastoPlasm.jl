"""
    update_Bᵢⱼ(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2) where {T1,T2}

Update material point velocities and positions from mesh nodes (PIC scheme, n-dim).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function Bij(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        # Bᵢⱼ update
        Bᵢⱼ,xₚ = zero(eltype(mpts.Bᵢⱼ)), mpts.x[p]
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N    = basis.N[nn,p]
            δ    = mesh.x[no] - xₚ
            Bᵢⱼ  = Bᵢⱼ + N * (mesh.s.v[no] * δ')
        end
        mpts.Bᵢⱼ[p] = Bᵢⱼ
    end
end