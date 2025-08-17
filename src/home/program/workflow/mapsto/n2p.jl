"""
    nd_n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,C_pf::T2) where {T1,T2}

Update material point velocities and positions from mesh nodes using PIC-FLIP scheme.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function picflip_n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,C_pf::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        for dim ∈ 1:mesh.dim
            # pic update
            δvPIC = T2(0.0)
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δvPIC += mpts.ϕ∂ϕ[nn,p,1]*mesh.s.v[dim,no]
            end
            # flip update
            δaFLIP = δvFLIP = T2(0.0)
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δaFLIP += mpts.ϕ∂ϕ[nn,p,1]*mesh.s.a[dim,no]
                δvFLIP += mpts.ϕ∂ϕ[nn,p,1]*mesh.s.v[dim,no]
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
    n2p(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2, instr::NamedTuple) where {T1,T2}

Map mesh node solution back to material points using the selected transfer kernel.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(ndrange=mpts.nmp,mpts,mesh,dt,T2(instr[:fwrk][:C_pf]));sync(CPU())
    # (for APIC) compute Bᵢⱼ for material points
    if instr[:fwrk][:trsfr] == "apic"
        instr[:cairn][:mapsto][:map].Bᵢⱼ!(ndrange=mpts.nmp,mpts,mesh);sync(CPU())
    end
    # (if musl) reproject nodal velocities
    if instr[:fwrk][:musl]
        # reset nodal quantities
        fill!(mesh.s.mv,T2(0.0))
        fill!(mesh.s.v ,T2(0.0))
        # accumulate material point contributions
        instr[:cairn][:mapsto][:augm].p2n!(ndrange=mpts.nmp,mpts,mesh);sync(CPU())
        # solve for nodal incremental displacement
        instr[:cairn][:mapsto][:augm].solve!(ndrange=mesh.nno[end],mesh);sync(CPU())
    elseif instr[:fwrk][:trsfr] == "apic"
        instr[:cairn][:mapsto][:map].Bᵢⱼ!(ndrange=mpts.nmp,mpts,mesh);sync(CPU())
    end
    return nothing
end