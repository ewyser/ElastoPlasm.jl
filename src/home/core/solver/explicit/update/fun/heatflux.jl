"""
    heatflux(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2})

Update material point velocities and positions from thermal-type mesh nodes using PIC-FLIP scheme.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function heatflux(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # compute heat flux
        for dim ∈ 1:(length(mesh.prprt.nel)-1)
            dQᵢ = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = basis.p2n[p][nn]
                if iszero(no) continue end
                dQᵢ += mpts.ϕ∂ϕ[nn,p,dim+1]*mesh.T[no]
            end    
            mpts.t.q[dim,p] = -mpts.t.k[p]*dQᵢ
        end
    end
end