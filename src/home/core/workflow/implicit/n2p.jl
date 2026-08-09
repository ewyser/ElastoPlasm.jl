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
@kernel inbounds = true function n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvx = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N*mesh.s.v[1,no]
        end
        mpts.s.v[1,p] = δvx
        # find maximum velocity component over mpts
        @atom mpts.vmax[1] = max(mpts.vmax[1],abs(δvx))
    end  
end
@kernel inbounds = true function n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N*mesh.s.v[1,no]
            δvy += N*mesh.s.v[2,no]
        end
        mpts.s.v[p] = SVector{2,T2}(δvx, δvy)
        mpts.s.u[p] = SVector{2,T2}(dt*δvx, dt*δvy)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(δvx))
        @atom mpts.vmax[2] = max(mpts.vmax[2], abs(δvy))
    end  
end
@kernel inbounds = true function n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D<:ThreeDimension}
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
        vx = C_pf*(v[1]+dt*δaxFLIP) + (T2(1.0)-C_pf)*δvxPIC
        vy = C_pf*(v[2]+dt*δayFLIP) + (T2(1.0)-C_pf)*δvyPIC
        vz = C_pf*(v[3]+dt*δazFLIP) + (T2(1.0)-C_pf)*δvzPIC
        mpts.s.v[p] = SVector{3,T2}(vx, vy, vz)
        mpts.s.u[p] = SVector{3,T2}(dt*δvxPIC, dt*δvyPIC, dt*δvzPIC)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(vx))
        @atom mpts.vmax[2] = max(mpts.vmax[2], abs(vy))
        @atom mpts.vmax[3] = max(mpts.vmax[3], abs(vz))
    end  
end
