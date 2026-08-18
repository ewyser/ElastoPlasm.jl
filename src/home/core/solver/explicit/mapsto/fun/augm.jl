"""
    augm_p2n(mpts::Point{T1,T2}, mesh::MeshSolidPhase{T1,T2})

Accumulate material point momentum to mesh nodes for DM augmentation.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure for solid phase.

# Returns
- Updates mesh fields in-place.
"""
@kernel inbounds = true function augm_p2n(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms = mpts.s.ρ[p]*mpts.Ω[p]
        mv = ms * mpts.s.v[p]
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N = basis.N[nn,p]
            @atom mesh.s.mv[1,no] += N * mv[1]
            @atom mesh.s.mv[2,no] += N * mv[2]
        end
    end
end
@kernel inbounds = true function augm_p2n(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},basis::Basis{T1,T2,3}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms = mpts.s.ρ[p]*mpts.Ω[p]
        mv = ms * mpts.s.v[p]
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N = basis.N[nn,p]
            @atom mesh.s.mv[1,no] += N * mv[1]
            @atom mesh.s.mv[2,no] += N * mv[2]
            @atom mesh.s.mv[3,no] += N * mv[3]
        end
    end
end
"""
    augm_solve(mesh::MeshSolidPhase{T1,T2}) where {T1,T2}

Update mesh node velocities for DM augmentation.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates mesh velocities in-place.
"""
@kernel inbounds = true function augm_solve(mesh::Mesh{T1,T2,D}) where {T1,T2,D}
    no = @index(Global)
    if no ≤ mesh.prprt.nno[end] 
        if iszero(mesh.s.m[no])
            nothing         
        else
            m⁻¹ = (T2(1.0)/mesh.s.m[no])
            TV = eltype(mesh.s.v)
            mesh.s.v[no] = TV(ntuple(Val(length(TV))) do d
                mesh.s.bcs.status[d,no] ? T2(0) : mesh.s.mv[d,no]*m⁻¹
            end)
        end
    end
end











"""
    augm_p2n(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2})

Accumulate material point heat capacity to mesh nodes for DM augmentation.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.

# Returns
- Updates mesh fields in-place.
"""
@kernel inbounds = true function augm_p2n(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # accumulation
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            @atom mesh.mcT[no]+= mpts.ϕ∂ϕ[nn,p,1] * mpts.s.ρ[p]*mpts.Ω[p] * mpts.t.c[p] * mpts.t.T[p]
        end
    end
end
@kernel inbounds = true function augm_solve(mesh::MeshThermalPhase{T1,T2}) where {T1,T2}
    no = @index(Global)
    if no ≤ mesh.prprt.nno[end] 
        if iszero(mesh.cᵢ[no])
            nothing         
        else
            # apply boundary contidions || forward euler solution
            if mesh.bcs.status[1,no]
                mesh.T[no] = T2(20.0)
            else
                mesh.T[no]  = mesh.mcT[no] * (T2(1.0)/mesh.cᵢ[no])
            end
        end
    end
end