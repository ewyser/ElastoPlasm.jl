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
                @atom mesh.mv[dim,no]+= mpts.ϕ∂ϕ[nn,p,1]*((mpts.s.ρ[p]*mpts.Ω[p])*mpts.s.v[dim,p])
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
                    mesh.v[dim,no] = mesh.mv[dim,no]*(T2(1.0)/mesh.mᵢ[no])  
                end
            end
        end
    end
end