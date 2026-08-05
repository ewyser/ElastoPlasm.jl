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
@views @kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx   = mpts.s.σᵢ[1,p]       
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            δx        = mesh.x[1,no]-mpts.x[1,p]
            if iszero(no) continue end
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p]) 
            else
                D⁻¹ =  mpts.δᵢⱼ
            end
            # accumulation
             mesh.m[no]     += N * ms
             mesh.mv[:,no] .+= N .* ms .* (mpts.s.v[:,p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * δx)
             mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω       = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx,σyy,σxy = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p]     ,mpts.s.σᵢ[3,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            δx, δy    = mesh.x[1,no]-mpts.x[1,p], mesh.x[2,no]-mpts.x[2,p]
            if iszero(no) continue end
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p]) 
            else
                D⁻¹ =  mpts.δᵢⱼ
            end
            # accumulation
             mesh.s.m[no]     += N * ms
             mesh.s.mv[:,no] .+= N .* ms .* (mpts.s.v[:,p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * vcat(δx,δy))
             mesh.s.oobf[1,no]-= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
             mesh.s.oobf[2,no]-= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms  ,Ω        = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        σxx ,σyy ,σzz = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p] ,mpts.s.σᵢ[3,p]
        σyx ,σzy ,σzx = mpts.s.σᵢ[6,p]       ,mpts.s.σᵢ[4,p] ,mpts.s.σᵢ[5,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no            = mpts.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3],mpts.ϕ∂ϕ[nn,p,4]
            # accumulation
            if iszero(no) continue end
            @atom mesh.m[no]    += N * ms
            if abs(det(mpts.Dᵢⱼ[:,:,p])) > T2(1e-12)
                D⁻¹ = inv(mpts.Dᵢⱼ[:,:,p]) 
            else
                D⁻¹ =  mpts.δᵢⱼ
            end
            mesh.mv[:,no] .+= N .* ms .* (mpts.s.v[:,p] .+ mpts.Bᵢⱼ[:,:,p] * D⁻¹ * mpts.Δnp[nn,:,p])
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (ms * g[3])
        end
    end
end






























