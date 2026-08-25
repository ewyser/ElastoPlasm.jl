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
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},basis::Basis{T1,T2,1},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        xp    = mpts.x[p]
        vp    = mpts.s.v[p]
        σxx   = get_voigt(mpts.s.σᵢ[p])[1]
        Bp    = mpts.Bᵢⱼ[p]
        Dp    = mpts.Dᵢⱼ[p]
        D⁻¹   = abs(det(Dp)) > T2(1e-12) ? inv(Dp) : SMatrix{1,1,T2}(I)
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(1))
            δ     = mesh.x[no] - xp
            mv    = N * ms * (vp + Bp * D⁻¹ * δ)
            # accumulation
            @atom mesh.s.m[no]   += N * ms
            @atom mesh.s.mv[no]  += mv[1]
            @atom mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω       = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        xp          = mpts.x[p]
        vp          = mpts.s.v[p]
        σ           = get_voigt(mpts.s.σᵢ[p])
        σxx,σyy,σxy = σ[1]                  ,σ[2]                ,σ[3]
        Bp          = mpts.Bᵢⱼ[p]
        Dp          = mpts.Dᵢⱼ[p]
        D⁻¹         = abs(det(Dp)) > T2(1e-12) ? inv(Dp) : SMatrix{2,2,T2}(I)
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(2))
            δ     = mesh.x[no] - xp
            mv    = N * ms * (vp + Bp * D⁻¹ * δ)
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += mv[1]
            @atom mesh.s.mv[2,no]  += mv[2]
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function apic_p2n(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},basis::Basis{T1,T2,3},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms  ,Ω        = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        xp            = mpts.x[p]
        vp            = mpts.s.v[p]
        σ             = get_voigt(mpts.s.σᵢ[p])
        σxx ,σyy ,σzz = σ[1]                  ,σ[2]            ,σ[3]
        σyx ,σzy ,σzx = σ[6]                  ,σ[4]            ,σ[5]
        Bp            = mpts.Bᵢⱼ[p]
        Dp            = mpts.Dᵢⱼ[p]
        D⁻¹           = abs(det(Dp)) > T2(1e-12) ? inv(Dp) : SMatrix{3,3,T2}(I)
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(3))
            ∂Nx,∂Ny,∂Nz = ∂N[1], ∂N[2], ∂N[3]
            δ     = mesh.x[no] - xp
            mv    = N * ms * (vp + Bp * D⁻¹ * δ)
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += mv[1]
            @atom mesh.s.mv[2,no]  += mv[2]
            @atom mesh.s.mv[3,no]  += mv[3]
            @atom mesh.s.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.s.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.s.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (ms * g[3])
        end
    end
end






























