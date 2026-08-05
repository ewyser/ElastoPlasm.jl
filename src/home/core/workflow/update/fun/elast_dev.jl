# Deviatoric stress update for u-P formulation: σ = σ' - p·δᵢⱼ (geomechanics: p > 0 = compression)
# σ' = σ'_n + 2G·ε'  where ε' uses the 3D deviatoric split (valid for plane strain: εv = εxx+εyy)

@kernel inbounds = true function elast_dev(mpts::Point{T1,T2,D},Gc::T2,p_n::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        Δεxx = mpts.s.ΔFᵢⱼ[1,1,p] - T2(1)
        Δεyy = mpts.s.ΔFᵢⱼ[2,2,p] - T2(1)
        Δεxy = T2(0.5) * (mpts.s.ΔFᵢⱼ[1,2,p] + mpts.s.ΔFᵢⱼ[2,1,p])
        # plane strain: volumetric = εxx+εyy (εzz=0), divide by 3 for deviatoric
        Δεv   = (Δεxx + Δεyy) / T2(3)
        Δεxx -= Δεv;  Δεyy -= Δεv
        twoG  = T2(2) * Gc
        p_cur = mpts.s.p[p]
        # strip p_n from σ_n to get σ'_n (σ = σ' - p_geo · I  ⇒  σ' = σ + p_geo)
        mpts.s.σᵢ[1,p] += p_n[p]
        mpts.s.σᵢ[2,p] += p_n[p]
        # add deviatoric increment and subtract new pressure
        mpts.s.σᵢ[1,p] += twoG * Δεxx - p_cur
        mpts.s.σᵢ[2,p] += twoG * Δεyy - p_cur
        mpts.s.σᵢ[3,p] += twoG * Δεxy
    end
end

@kernel inbounds = true function elast_dev(mpts::Point{T1,T2,D},Gc::T2,p_n::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        Δεxx = mpts.s.ΔFᵢⱼ[1,1,p] - T2(1)
        Δεyy = mpts.s.ΔFᵢⱼ[2,2,p] - T2(1)
        Δεzz = mpts.s.ΔFᵢⱼ[3,3,p] - T2(1)
        Δεxy = T2(0.5) * (mpts.s.ΔFᵢⱼ[1,2,p] + mpts.s.ΔFᵢⱼ[2,1,p])
        Δεyz = T2(0.5) * (mpts.s.ΔFᵢⱼ[2,3,p] + mpts.s.ΔFᵢⱼ[3,2,p])
        Δεxz = T2(0.5) * (mpts.s.ΔFᵢⱼ[1,3,p] + mpts.s.ΔFᵢⱼ[3,1,p])
        Δεv   = (Δεxx + Δεyy + Δεzz) / T2(3)
        Δεxx -= Δεv;  Δεyy -= Δεv;  Δεzz -= Δεv
        twoG  = T2(2) * Gc
        p_cur = mpts.s.p[p]
        mpts.s.σᵢ[1,p] += p_n[p];  mpts.s.σᵢ[2,p] += p_n[p];  mpts.s.σᵢ[3,p] += p_n[p]
        mpts.s.σᵢ[1,p] += twoG * Δεxx - p_cur
        mpts.s.σᵢ[2,p] += twoG * Δεyy - p_cur
        mpts.s.σᵢ[3,p] += twoG * Δεzz - p_cur
        mpts.s.σᵢ[4,p] += twoG * Δεyz
        mpts.s.σᵢ[5,p] += twoG * Δεxz
        mpts.s.σᵢ[6,p] += twoG * Δεxy
    end
end
