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
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvxPIC = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvxPIC += N*mesh.s.v[1,no]
        end
        # flip update
        δaxFLIP = T2(0.0)
        δvxFLIP = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δaxFLIP += N*mesh.s.a[1,no]
            δvxFLIP += N*mesh.s.v[1,no]
        end
        # picflip update for material point's velocity and position
        vx_new = C_pf*(mpts.s.v[p][1]+dt*δaxFLIP) + (T2(1.0)-C_pf)*δvxPIC
        mpts.s.v[p] = SVector{1,T2}(vx_new)
        mpts.s.u[p] = SVector{1,T2}(dt*δvxPIC)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(vx_new))
    end  
end
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvxPIC = δvyPIC = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvxPIC += N*mesh.s.v[1,no]
            δvyPIC += N*mesh.s.v[2,no]
        end
        # flip update
        δaxFLIP = δayFLIP = T2(0.0)
        δvxFLIP = δvyFLIP = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δaxFLIP += N*mesh.s.a[1,no]
            δayFLIP += N*mesh.s.a[2,no]
            δvxFLIP += N*mesh.s.v[1,no]
            δvyFLIP += N*mesh.s.v[2,no]
        end
        # picflip update for material point's velocity and position
        v = mpts.s.v[p]
        vx_new = C_pf*(v[1]+dt*δaxFLIP) + (T2(1.0)-C_pf)*δvxPIC
        vy_new = C_pf*(v[2]+dt*δayFLIP) + (T2(1.0)-C_pf)*δvyPIC
        mpts.s.v[p] = SVector{2,T2}(vx_new, vy_new)
        mpts.s.u[p] = SVector{2,T2}(dt*δvxPIC, dt*δvyPIC)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(vx_new))
        @atom mpts.vmax[2] = max(mpts.vmax[2], abs(vy_new))
    end  
end
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvxPIC = δvyPIC = δvzPIC = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvxPIC += N*mesh.s.v[1,no]
            δvyPIC += N*mesh.s.v[2,no]
            δvzPIC += N*mesh.s.v[3,no]
        end
        # flip update
        δaxFLIP = δayFLIP = δazFLIP = T2(0.0)
        δvxFLIP = δvyFLIP = δvzFLIP = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δaxFLIP += N*mesh.s.a[1,no]
            δayFLIP += N*mesh.s.a[2,no]
            δazFLIP += N*mesh.s.a[3,no]
            δvxFLIP += N*mesh.s.v[1,no]
            δvyFLIP += N*mesh.s.v[2,no]
            δvzFLIP += N*mesh.s.v[3,no]
        end
        # picflip update for material point's velocity and position
        v = mpts.s.v[p]
        vx_new = C_pf*(v[1]+dt*δaxFLIP) + (T2(1.0)-C_pf)*δvxPIC
        vy_new = C_pf*(v[2]+dt*δayFLIP) + (T2(1.0)-C_pf)*δvyPIC
        vz_new = C_pf*(v[3]+dt*δazFLIP) + (T2(1.0)-C_pf)*δvzPIC
        mpts.s.v[p] = SVector{3,T2}(vx_new, vy_new, vz_new)
        mpts.s.u[p] = SVector{3,T2}(dt*δvxPIC, dt*δvyPIC, dt*δvzPIC)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(vx_new))
        @atom mpts.vmax[2] = max(mpts.vmax[2], abs(vy_new))
        @atom mpts.vmax[3] = max(mpts.vmax[3], abs(vz_new))
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
        C_pf = T2(0.0)
        mpts.t.T[p] = C_pf*(mpts.t.T[p]+dt*δTFLIP) + (T2(1.0)-C_pf)*δTPIC
        #mpts.t.T[p] = δTPIC
    end  
end