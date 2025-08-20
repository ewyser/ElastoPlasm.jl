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
@kernel inbounds = true function Bij(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp    
        # Bᵢⱼ update
        Bᵢⱼ = zeros(T2,mesh.prprt.dim,mesh.prprt.dim)
        for nn ∈ 1:mesh.prprt.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            Bᵢⱼ.+= mpts.ϕ∂ϕ[nn,p,1] .* (mesh.v[:,no] * mpts.Δnp[nn,:,p]')
        end
        mpts.Bᵢⱼ[:,:,p].= Bᵢⱼ
    end
end