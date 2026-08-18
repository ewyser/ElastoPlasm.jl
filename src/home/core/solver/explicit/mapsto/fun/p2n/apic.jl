# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# APIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    apic_1d_p2n(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, g::Vector{T2}) where {T1,T2}

Project 1D material point data to mesh nodes (APIC scheme).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.

# Returns
- Updates mesh fields in-place.
"""
@views @kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},basis::Basis{T1,T2,1},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx   = mpts.s.σᵢ[p][1]
        for nn ∈ 1:mesh.prprt.nn
            no        = basis.p2n[p][nn]
            δx        = mesh.x[no][1]-mpts.x[p][1]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(1))
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p])
            else
                D⁻¹ = SMatrix{1,1,T2}(I)
            end
            # accumulation
             mesh.m[no]     += N * ms
             mesh.mv[:,no] .+= N .* ms .* (mpts.s.v[p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * δx)
             mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω       = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx,σyy,σxy = mpts.s.σᵢ[p][1]       ,mpts.s.σᵢ[p][2]     ,mpts.s.σᵢ[p][3]
        for nn ∈ 1:mesh.prprt.nn
            no        = basis.p2n[p][nn]
            δx, δy    = mesh.x[no][1]-mpts.x[p][1], mesh.x[no][2]-mpts.x[p][2]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(2))
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p])
            else
                D⁻¹ = SMatrix{2,2,T2}(I)
            end
            # accumulation
             mesh.s.m[no]     += N * ms
             mesh.s.mv[:,no] .+= N .* ms .* (mpts.s.v[p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * vcat(δx,δy))
             mesh.s.oobf[1,no]-= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
             mesh.s.oobf[2,no]-= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,3},mesh::MeshSolidPhase{T1,T2,3},basis::Basis{T1,T2,3},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms  ,Ω        = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx ,σyy ,σzz = mpts.s.σᵢ[p][1]       ,mpts.s.σᵢ[p][2] ,mpts.s.σᵢ[p][3]
        σyx ,σzy ,σzx = mpts.s.σᵢ[p][6]       ,mpts.s.σᵢ[p][4] ,mpts.s.σᵢ[p][5]
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no            = basis.p2n[p][nn]
            N,∂Nx,∂Ny,∂Nz = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3],mpts.ϕ∂ϕ[nn,p,4]
            # accumulation
            if iszero(no) continue end
            @atom mesh.m[no]    += N * ms
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p])
            else
                D⁻¹ = SMatrix{3,3,T2}(I)
            end
            mesh.mv[:,no] .+= N .* ms .* (mpts.s.v[p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * mpts.Δnp[nn,:,p])
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (ms * g[3])
        end
    end
end






























