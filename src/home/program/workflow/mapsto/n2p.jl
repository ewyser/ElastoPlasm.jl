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

    if p ≤ mpts.nmp 
        # compute heat flux
        for dim ∈ 1:mesh.prprt.dim
            dQᵢ = T2(0.0)
            for nn ∈ 1:mesh.prprt.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                dQᵢ += mpts.ϕ∂ϕ[nn,p,dim+1]*mesh.T[no]
            end    
            mpts.t.q[dim,p] = -mpts.t.k[p]*dQᵢ
        end
    end
end
"""
    n2p(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},dt::T2,instr::NamedTuple)

Map solid-type mesh node solution back to material points using the selected transfer kernel.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure for solid phase.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function n2p(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(mpts,mesh,dt,T2(instr[:fwrk][:C_pf]); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal velocities
    if instr[:fwrk][:musl]
        # reset nodal quantities
        fill!(mesh.mv,T2(0.0))
        fill!(mesh.v ,T2(0.0))
        # accumulate material point contributions
        instr[:cairn][:mapsto][:augm].p2n!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # solve for nodal incremental displacement
        instr[:cairn][:mapsto][:augm].solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    # (for APIC) compute Bᵢⱼ for material points
    if instr[:fwrk][:trsfr] == "apic"
        instr[:cairn][:mapsto][:map].Bᵢⱼ!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    end
    #=
    instr[:cairn][:mapsto][:map].n2p!(mpts,mesh.t,dt,T2(instr[:fwrk][:C_pf]); ndrange=mpts.nmp);sync(CPU())
    =#
    return nothing
end
"""
    n2p(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,instr::NamedTuple)

Map solid-type mesh node solution back to material points using the selected transfer kernel.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function n2p(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(mpts,mesh,dt,T2(instr[:fwrk][:C_pf]); ndrange=mpts.nmp);sync(CPU())
    return nothing
end