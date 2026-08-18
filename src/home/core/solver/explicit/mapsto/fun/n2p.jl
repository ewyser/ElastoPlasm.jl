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
# Single generic kernel replaces the former 1D/2D/3D variants.
# SVector arithmetic lets us accumulate whole nodal vectors in one pass,
# eliminating the redundant δvFLIP variables and the duplicate loop.
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},dt::T2,C_pf::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        TV    = eltype(mesh.s.v)
        vPIC  = zero(TV)
        aFLIP = zero(TV)
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N      = basis.N[p][nn]
            vPIC  += N * mesh.s.v[no]
            aFLIP += N * mesh.s.a[no]
        end
        v_new = C_pf * (mpts.s.v[p] + dt * aFLIP) + (T2(1.0) - C_pf) * vPIC
        mpts.s.v[p] = v_new
        mpts.s.u[p] = dt * vPIC
        for d ∈ 1:length(TV)
            @atom mpts.vmax[d] = max(mpts.vmax[d], abs(v_new[d]))
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
# Single loop: T and dT share nodes so PIC+FLIP can be merged.
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp
        δTPIC  = T2(0.0)
        δTFLIP = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            w = mpts.ϕ∂ϕ[nn,p,1]
            δTPIC  += w * mesh.T[no]
            δTFLIP += w * mesh.dT[no]
        end
        C_pf = T2(0.0)
        mpts.t.T[p] = C_pf*(mpts.t.T[p]+dt*δTFLIP) + (T2(1.0)-C_pf)*δTPIC
    end
end