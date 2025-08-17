"""
    nd_n2p(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},dt::T2,C_pf::T2)

Update material point velocities and positions from solid-type mesh nodes using PIC-FLIP scheme.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure for solid phase.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        for dim ∈ 1:mesh.prprt.dim
            # pic update
            δvPIC = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δvPIC += mpts.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no]
            end
            # flip update
            δaFLIP = δvFLIP = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δaFLIP += mpts.ϕ∂ϕ[nn,p,1]*mesh.a[dim,no]
                δvFLIP += mpts.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no]
            end
        # picflip update for material point's velocity and position
        mpts.s.v[dim,p] = C_pf*(mpts.s.v[dim,p]+dt*δaFLIP) + (T2(1.0)-C_pf)*δvPIC
        mpts.x[dim,p]  += dt*δvPIC
        # find maximum velocity component over mps
        @atom mpts.vmax[dim] = max(mpts.vmax[dim],abs(mpts.s.v[dim,p]))
        end
    end  
end
"""
    nd_n2p(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,C_pf::T2)

Update material point velocities and positions from thermal-type mesh nodes using PIC-FLIP scheme.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        for dim ∈ 1:mesh.prprt.dim
            # pic update
            δTPIC = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δTPIC += mpts.ϕ∂ϕ[nn,p,1]*mesh.T[no]
            end
            # flip update
            δTFLIP = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δTFLIP += mpts.ϕ∂ϕ[nn,p,1]*mesh.dT[no]
            end
            # picflip update for material point's temperature
            if dim == 1
                mpts.t.T[p] = C_pf*(mpts.t.T[p]+dt*δTFLIP) + (T2(1.0)-C_pf)*δTPIC
            end 
        end
    end  
end