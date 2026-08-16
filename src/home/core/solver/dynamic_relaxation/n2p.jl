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
@kernel inbounds = true function n2p(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvx = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N*mesh.s.v[no][1]
        end
        mpts.s.v[1,p] = δvx
        # find maximum velocity component over mpts
        @atom mpts.vmax[1] = max(mpts.vmax[1],abs(δvx))
    end  
end
@kernel inbounds = true function n2p(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N*mesh.s.v[no][1]
            δvy += N*mesh.s.v[no][2]
        end
        mpts.s.v[p] = SVector{2,T2}(δvx, δvy)
        mpts.s.u[p] = SVector{2,T2}(dt*δvx, dt*δvy)
        @atom mpts.vmax[1] = max(mpts.vmax[1], abs(δvx))
        @atom mpts.vmax[2] = max(mpts.vmax[2], abs(δvy))
    end  
end
@kernel inbounds = true function n2p(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        δvxPIC = δvyPIC = δvzPIC = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvxPIC += N*mesh.s.v[no][1]
            δvyPIC += N*mesh.s.v[no][2]
            δvzPIC += N*mesh.s.v[no][3]
        end
        # flip update
        δaxFLIP = δayFLIP = δazFLIP = T2(0.0)
        δvxFLIP = δvyFLIP = δvzFLIP = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δaxFLIP += N*mesh.s.a[no][1]
            δayFLIP += N*mesh.s.a[no][2]
            δazFLIP += N*mesh.s.a[no][3]
            δvxFLIP += N*mesh.s.v[no][1]
            δvyFLIP += N*mesh.s.v[no][2]
            δvzFLIP += N*mesh.s.v[no][3]
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
