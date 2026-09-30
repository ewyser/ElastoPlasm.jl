@inline function materialParam(ϕ::T2,ψ::T2,c::T2,D::Int) where {T2}
    if D == 2   # plane strain: plane-strain match to Mohr-Coulomb
        η   = T2(3.0)*tan(ϕ)/(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
        ηB  = T2(3.0)*tan(ψ)/(sqrt(T2(9.0)+T2(12.0)*tan(ψ)*tan(ψ)))
        ξ   = T2(3.0)*c     /(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
    else
        η   = T2(6.0)*sin(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ)))
        ηB  = T2(6.0)*sin(ψ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ψ)))
        ξ   = T2(6.0)*c*cos(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ)))
    end
    return η,ηB,ξ
end

"""
    _druckerprager_return_map(σ::AbstractStress, ϵpII::MVector{2,T}, cmp) -> (Δλ, p, dev)

Closed-form Drucker-Prager return mapping on the stored pressure/deviator split — shared by both
`drucker_prager` methods below, since the yield-surface algebra doesn't care whether the stress is
Cauchy or Kirchhoff. Returns `σ`'s own `(p, dev)` (elastic) unless yielding, in which case the
pressure is returned along the flow direction and the deviator scaled radially onto the cone.
`ϵpII` (the non-local-regularized hardening history) is mutated in place, matching the
pre-existing `MVector` scratch convention.
"""
@inline function _druckerprager_return_map(σ::AbstractStress{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel{T,D}) where {T,D}
    Δλ       = T(0.0)
    ψ        = T(0.0*π/180.0)
    ϕ₀,c₀,cᵣ = cmp.ϕ₀,cmp.c₀,cmp.cᵣ

    # closed-form solution return-mapping for D-P
    c = c₀+cmp.Hp*ϵpII[2]
    if c<cᵣ
        c = cᵣ
    end
    P,τII    = -σ.p,get_τII(σ)   # mean stress (positive in tension) and deviator norm √J₂
    η,ηB,ξ   = materialParam(ϕ₀,ψ,c,D)
    σm       = ξ/η
    fs       = τII+η*P-ξ
    p,dev    = σ.p,σ.dev
    if fs>T(0.0) && P<σm
        Δλ       = fs/(cmp.Gc+cmp.Kc*η*ηB)
        Pn,τn    = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
        p,dev    = -Pn,σ.dev*(τn/τII)
        ϵpII[1] += Δλ*sqrt(T(1/3)+T(2/9)*ηB^2)
    end
    # tension cutoff sits at the apex (σt = σm = ξ/η, Huang et al. 2015 Eq. 19 with σt = σt_max),
    # so the apex σm·I is the only admissible state for P ≥ σm (fs ≤ 0 forces τ = 0, ft ≤ 0 forces
    # mean ≤ σm): return there, Δλ from Eq. 54, ϵpII increment from Eq. 56
    if P≥σm
        Δλ       = (P-σm)/cmp.Kc
        p,dev    = -σm,zero(σ.dev)
        ϵpII[1] += sqrt(T(2.0))*Δλ/T(3.0)
    end
    return Δλ,p,dev
end

"""
    drucker_prager(τᵢ::KirchhoffStress, ϵᵢ::LogarithmicStrain, ϵpII, cmp) -> (ϵᵢ, τᵢ, Δλ, ϵpII)
    drucker_prager(σᵢ::CauchyStress, ϵᵢ::InfinitesimalStrain, ϵpII, cmp) -> (ϵᵢ, σᵢ, Δλ, ϵpII)

Drucker-Prager return mapping, dispatched on stress/strain type — one method per
`(stress,strain)` pair, sharing the closed-form `_druckerprager_return_map` core rather
than duplicating it.

The finite-strain (`KirchhoffStress`/`LogarithmicStrain`) method additionally rebuilds the
elastic strain when yielding, with `_elastic_strain` (`stress.jl`). The infinitesimal-strain (`CauchyStress`/`InfinitesimalStrain`) method never
updates the strain — infinitesimal strain is tracked incrementally by `elast!` itself,
not re-derived from the return-mapped stress.
"""
@inline function drucker_prager(τᵢ::KirchhoffStress{T},ϵᵢ::LogarithmicStrain{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {T}
    Δλ,p,dev = _druckerprager_return_map(τᵢ,ϵpII,cmp)
    if Δλ > T(0.0)
        τᵢ = KirchhoffStress(p,dev)
        ϵᵢ = _elastic_strain(τᵢ,cmp)
    end
    return ϵᵢ,τᵢ,Δλ,ϵpII
end
@inline function drucker_prager(σᵢ::CauchyStress{T},ϵᵢ::InfinitesimalStrain{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {T}
    Δλ,p,dev = _druckerprager_return_map(σᵢ,ϵpII,cmp)
    if Δλ > T(0.0)
        σᵢ = CauchyStress(p,dev)
    end
    return ϵᵢ,σᵢ,Δλ,ϵpII
end

"""
    retmap(mpts::Point{T1,T2,D,CM,ST}) where {CM<:DruckerPrager, ST<:LogarithmicStrain}
    retmap(mpts::Point{T1,T2,D,CM,ST}) where {CM<:DruckerPrager, ST<:InfinitesimalStrain}

Drucker-Prager plastic corrector, dispatched on both `CM` and `ST` — mirrors `elast!`'s
existing `ST`-only dispatch pattern, extended to a second axis so `DP.jl`/`J2.jl` can
both contribute methods to one shared kernel name (`retmap`) instead of each exposing
separately-named `finite_*`/`infinitesimal_*` kernels that `init_update` had to pick
between via a `material.plastic`/`material.elastic` string branch.
"""
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM<:DruckerPrager,ST<:LogarithmicStrain}
    p = @index(Global)
    if p≤mpts.nmp
        # reset the plastic multiplier on every step: it is *the* activity gate read by
        # `nonlocal.jl` (`mpts.s.Δλ[p] != 0`), so leaving last step's value in place on an
        # elastic step keeps a particle spuriously "yielding" forever.
        mpts.s.Δλ[p] = T2(0.0)
        ϵᵢ,τᵢ,Δλ,ϵpII = drucker_prager(mpts.s.τᵢ[p],mpts.s.ϵᵢ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.ϵᵢ[p]  = ϵᵢ
            mpts.s.τᵢ[p]  = τᵢ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        end
    end
end
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM<:DruckerPrager,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p≤mpts.nmp
        mpts.s.Δλ[p] = T2(0.0)
        ϵᵢ,σᵢ,Δλ,ϵpII = drucker_prager(mpts.s.σᵢ[p],mpts.s.ϵᵢ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.σᵢ[p]  = σᵢ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        end
    end
end
