"""
    augm_momentum(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}) where {T1,T2}

Accumulate material point momentum to mesh nodes for DM augmentation.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates mesh fields in-place.
"""
@kernel inbounds = true function augm_momentum(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp
        for dim ∈ 1:mesh.dim 
            # accumulation
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                @atom mesh.p[dim,no]+= mpts.ϕ∂ϕ[nn,p,1]*(mpts.s.m[p]*mpts.s.v[dim,p])
            end
        end
    end
end
"""
    augm_velocity(mesh::Mesh{T1,T2}) where {T1,T2}

Update mesh node velocities for DM augmentation.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates mesh velocities in-place.
"""
@kernel inbounds = true function augm_velocity(mesh::Mesh{T1,T2}) where {T1,T2}
    no = @index(Global)
    if no≤mesh.nno[end] 
        if iszero(mesh.mᵢ[no])
            nothing         
        else
            for dim ∈ 1:mesh.dim       
                # apply boundary contidions || forward euler solution
                if mesh.bcs.status[dim,no]
                    mesh.v[dim,no] = T2(0.0)                                         
                else
                    mesh.v[dim,no] = mesh.p[dim,no]*(T2(1.0)/mesh.mᵢ[no])  
                end
            end
        end
    end
end
"""
    augm_displacement(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2) where {T1,T2}

Update material point displacements from mesh node velocities for DM augmentation.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.

# Returns
- Updates material point displacements in-place.
"""
@views @kernel inbounds = true function augm_displacement(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    # flip update
    if p≤mpts.nmp
        for dim ∈ 1:mesh.dim 
            Δu = T2(0.0)
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                Δu += dt*(mpts.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no])
            end
            mpts.s.u[dim,p]+= Δu
        end
    end
end
"""
    augm(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2, instr::NamedTuple) where {T1,T2}

Perform DM augmentation: accumulate, solve, and update displacements.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function augm(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # initialize for DM
    mesh.p.= T2(0.0)
    mesh.v.= T2(0.0)
    # accumulate material point contributions
    instr[:cairn][:mapsto][:augm].p2n!(ndrange=mpts.nmp,mpts,mesh);sync(CPU())
    # solve for nodal incremental displacement
    instr[:cairn][:mapsto][:augm].solve!(ndrange=mesh.nno[end],mesh);sync(CPU())
    # update material point's displacement
    instr[:cairn][:mapsto][:augm].Δu!(ndrange=mpts.nmp,mpts,mesh,dt);sync(CPU())
    return nothing
end