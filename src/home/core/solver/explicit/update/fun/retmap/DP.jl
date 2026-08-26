@inline function σTr(σ0::SVector{3,T}) where {T}
    P   = (σ0[1]+σ0[2])/T(2.0)
    τ0  = σ0 .- SVector{3,T}(P,P,zero(T))
    τII = sqrt(T(0.5) * (τ0[1]^2 + τ0[2]^2) + τ0[3]^2)
    return P,τ0,τII
end
@inline function σTr(σ0::SVector{6,T}) where {T}
    P   = (σ0[1]+σ0[2]+σ0[3])/T(3.0)
    τ0  = σ0 .- SVector{6,T}(P,P,P,zero(T),zero(T),zero(T))
    τII = sqrt(T(0.5)*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    return P,τ0,τII
end
@inline function materialParam(ϕ::T2,ψ::T2,c::T2,nstr::T1) where {T1,T2}
    if nstr == T1(3)
        η   = T2(3.0)*tan(ϕ)/(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
        ηB  = T2(3.0)*tan(ψ)/(sqrt(T2(9.0)+T2(12.0)*tan(ψ)*tan(ψ)))
        ξ   = T2(3.0)*c     /(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
    elseif nstr == T1(6)
        η   = T2(6.0)*sin(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ)))
        ηB  = T2(6.0)*sin(ψ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ψ)))
        ξ   = T2(6.0)*c*cos(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ)))
    end
    return η,ηB,ξ
end

@inline function σn(Pn::T,τ0::SVector{3,T},τn::T,τII::T) where {T}
    return τ0 .* (τn / τII) .+ SVector{3,T}(Pn,Pn,zero(T))
end
@inline function σn(Pn::T,τ0::SVector{6,T},τn::T,τII::T) where {T}
    return τ0 .* (τn / τII) .+ SVector{6,T}(Pn,Pn,Pn,zero(T),zero(T),zero(T))
end

"""
    _druckerprager_return_map(σᵢ::SVector, ϵpII::MVector{2,T}, cmp) -> (Δλ, σᵢ)

Closed-form Drucker-Prager return mapping, in Voigt space — shared by both
`drucker_prager` methods below, since the yield-surface algebra itself doesn't care
whether the input Voigt vector came from a `CauchyStress` or a `KirchhoffStress`.
Returns `σᵢ` unchanged (elastic) unless yielding, in which case it's the return-mapped
deviator+pressure reconstructed via `σn`. `ϵpII` (the non-local-regularized hardening
history) is mutated in place, matching the pre-existing `MVector` scratch convention.
"""
@inline function _druckerprager_return_map(σᵢ::SVector{N,T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {N,T}
    Δλ       = T(0.0)
    ψ,nstr   = T(0.0*π/180.0),length(σᵢ)
    ϕ₀,c₀,cᵣ = cmp.ϕ₀,cmp.c₀,cmp.cᵣ

    # closed-form solution return-mapping for D-P
    c = c₀+cmp.Hp*ϵpII[2]
    if c<cᵣ
        c = cᵣ
    end
    P,τ0,τII = σTr(σᵢ)
    η,ηB,ξ   = materialParam(ϕ₀,ψ,c,nstr)
    σm,τP    = ξ/η,ξ-η*(ξ/η)
    fs,ft    = τII+η*P-ξ,P-σm
    αP,h     = sqrt(T(1.0)+η^2)-η,τII-τP-(sqrt(T(1.0)+η^2))*(P-σm)
    σout     = σᵢ
    if fs>T(0.0) && P<σm || h>T(0.0) && P≥σm
        Δλ       = fs/(cmp.Gc+cmp.Kc*η*ηB)
        Pn,τn    = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
        σout     = σn(Pn,τ0,τn,τII)
        ϵpII[1] += Δλ*sqrt(T(1/3)+T(2/9)*ηB^2)
    end
    if h≤T(0.0) && P≥σm
        Δλ       = (P-σm)/cmp.Kc
        Pn       = σm-P
        σout     = σn(Pn,τ0,T(0.0),τII)
        ϵpII[1] += sqrt(T(2.0))*Δλ/T(3.0)
    end
    return Δλ,σout
end

"""
    drucker_prager(τᵢⱼ::KirchhoffStress, ϵᵢⱼ::LogarithmicStrain, ϵpII, cmp) -> (ϵᵢⱼ, τᵢⱼ, Δλ, ϵpII)
    drucker_prager(σᵢⱼ::CauchyStress, ϵᵢⱼ::InfinitesimalStrain, ϵpII, cmp) -> (ϵᵢⱼ, σᵢⱼ, Δλ, ϵpII)

Drucker-Prager return mapping, dispatched on stress/strain type — one method per
`(stress,strain)` pair, sharing the closed-form `_druckerprager_return_map` core rather
than duplicating it. Both methods do the algebra in Voigt space (`get_voigt`) using the
pre-existing `σTr`/`σn` helpers — the stored `(p,dev)` split is representational rather
than canonically trace-free (see the `AbstractTensor` docstring), so the invariants are
re-derived from the Voigt view.

The finite-strain (`KirchhoffStress`/`LogarithmicStrain`) method additionally updates
the strain when yielding, built directly by `LogarithmicStrain(cmp.Del\\τᵢ)` — the
compliance solve `Del\\τᵢ` is already an engineering-Voigt strain vector, and that
constructor (`tensor.jl`) does the engineering→tensor conversion and vol/dev split in
one step. The infinitesimal-strain (`CauchyStress`/`InfinitesimalStrain`) method never
updates the strain — infinitesimal strain is tracked incrementally by `elast!` itself,
not re-derived from the return-mapped stress.
"""
@inline function drucker_prager(τᵢⱼ::KirchhoffStress{S,T,L},ϵᵢⱼ::LogarithmicStrain{S,T,L},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {S,T,L}
    Δλ, τᵢ  = _druckerprager_return_map(get_voigt(τᵢⱼ),ϵpII,cmp)
    if Δλ > T(0.0)
        τᵢⱼ = KirchhoffStress(τᵢ)
        ϵᵢⱼ = LogarithmicStrain(cmp.Del\τᵢ)
    end
    return ϵᵢⱼ,τᵢⱼ,Δλ,ϵpII
end
@inline function drucker_prager(σᵢⱼ::CauchyStress{S,T,L},ϵᵢⱼ::InfinitesimalStrain{S,T,L},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {S,T,L}
    Δλ,σᵢ = _druckerprager_return_map(get_voigt(σᵢⱼ),ϵpII,cmp)
    if Δλ > T(0.0)
        σᵢⱼ = CauchyStress(σᵢ)
    end
    return ϵᵢⱼ,σᵢⱼ,Δλ,ϵpII
end

"""
    retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {CM<:DruckerPrager, ST<:LogarithmicStrain}
    retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {CM<:DruckerPrager, ST<:InfinitesimalStrain}

Drucker-Prager plastic corrector, dispatched on both `CM` and `ST` — mirrors `elast!`'s
existing `ST`-only dispatch pattern, extended to a second axis so `DP.jl`/`J2.jl` can
both contribute methods to one shared kernel name (`retmap`) instead of each exposing
separately-named `finite_*`/`infinitesimal_*` kernels that `init_update` had to pick
between via a `plast.constitutive`/`strain.deform` string branch.
"""
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM<:DruckerPrager,TM,TV,TS,ST<:LogarithmicStrain}
    p = @index(Global)
    if p≤mpts.nmp
        # reset the plastic multiplier on every step: it is *the* activity gate read by
        # `nonlocal.jl` (`mpts.s.Δλ[p] != 0`), so leaving last step's value in place on an
        # elastic step keeps a particle spuriously "yielding" forever.
        mpts.s.Δλ[p] = T2(0.0)
        ϵᵢⱼ,τᵢⱼ,Δλ,ϵpII = drucker_prager(mpts.s.τᵢⱼ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.ϵᵢⱼ[p]  = ϵᵢⱼ
            mpts.s.τᵢⱼ[p]  = τᵢⱼ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        end
    end
end
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM<:DruckerPrager,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p≤mpts.nmp
        mpts.s.Δλ[p] = T2(0.0)
        ϵᵢⱼ,σᵢⱼ,Δλ,ϵpII = drucker_prager(mpts.s.σᵢⱼ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.σᵢⱼ[p]  = σᵢⱼ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        end
    end
end
