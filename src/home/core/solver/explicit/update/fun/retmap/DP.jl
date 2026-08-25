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
    drucker_prager(τ::AbstractStress, ϵ::AbstractStrain, ϵpII, cmp) -> (ϵ, τ, Δλ, ϵpII)

Closed-form Drucker-Prager return mapping. Consumes/produces the typed
`KirchhoffStress`/`LogarithmicStrain` stored on `PointSolidPhase`, but does the algebra
itself in Voigt space (`get_vector`) using the pre-existing `σTr`/`σn` helpers — the
stored `(p,dev)` split is representational rather than canonically trace-free (see the
`AbstractTensor` docstring), so the invariants are re-derived from the Voigt view.
`σn` naturally returns a pressure `Pn` and a scaled deviator, so the returned stress is
rebuilt from those directly rather than re-split from a Voigt vector.
"""
@inline function drucker_prager(τ::AbstractStress{S,T,L},ϵ::AbstractStrain{S,T,L},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel) where {S,T,L}
        τᵢ     = get_vector(τ)
        ϵᵢⱼ    = get_tensor(ϵ)
        Δλ     = T(0.0)
        ψ,nstr = T(0.0*π/180.0),length(τᵢ)
        ϕ₀,c₀,cᵣ = cmp.ϕ₀,cmp.c₀,cmp.cᵣ

        # closed-form solution return-mapping for D-P
        c = c₀+cmp.Hp*ϵpII[2]
        if c<cᵣ
            c = cᵣ
        end
        P,τ0,τII = σTr(τᵢ)
        η,ηB,ξ   = materialParam(ϕ₀,ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T(1.0)+η^2)-η,τII-τP-(sqrt(T(1.0)+η^2))*(P-σm)
        τout = τ
        if fs>T(0.0) && P<σm || h>T(0.0) && P≥σm
            Δλ       = fs/(cmp.Gc+cmp.Kc*η*ηB)
            Pn,τn    = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            τᵢ       = σn(Pn,τ0,τn,τII)
            τout     = KirchhoffStress(-Pn, τ0 .* (τn / τII))
            ϵpII[1] += Δλ*sqrt(T(1/3)+T(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ      = (P-σm)/cmp.Kc
            Pn      = σm-P
            τᵢ      = (σn(Pn,τ0,T(0.0),τII))
            τout    = KirchhoffStress(-Pn, τ0 .* (T(0.0) / τII))
            ϵpII[1]+= sqrt(T(2.0))*Δλ/T(3.0)
        end
        # update logarithmic strain tensor
        ϵout = ϵ
        if Δλ > T(0.0)
            ϵᵢⱼ  = mutate(cmp.Del\τᵢ,T(0.5),:tensor)
            ϵout = LogarithmicStrain(SMatrix{S,S,T,L}(ϵᵢⱼ))
        end
        return ϵout,τout,Δλ,ϵpII
end

@views @kernel inbounds = true function finite_DP(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp
        # reset the plastic multiplier on every step: it is *the* activity gate read by
        # `nonlocal.jl` (`mpts.s.Δλ[p] != 0`), so leaving last step's value in place on an
        # elastic step kept a particle spuriously "yielding" forever. `infinitesimal_DP`
        # below always did this; `finite_DP` did not — the zeroing was clearly intended,
        # it just sat in the dead commented-out block this kernel used to carry.
        mpts.s.Δλ[p] = T2(0.0)
        ϵᵢⱼ,τᵢ,Δλ,ϵpII = drucker_prager(mpts.s.τᵢ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[1,p],mpts.s.ϵpII[2,p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.ϵᵢⱼ[p]     = ϵᵢⱼ
            mpts.s.τᵢ[p]      = τᵢ
            mpts.s.Δλ[p]      = Δλ
            mpts.s.ϵpII[:,p] .= ϵpII
        end
    end
end
@views @kernel inbounds = true function infinitesimal_DP(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp
        cmp = mpts.s.cmp[p]
        mpts.s.Δλ[p] = T2(0.0)
        σᵢ       = get_vector(mpts.s.σᵢ[p])
        ψ,nstr   = T2(0.0*π/180.0),length(σᵢ)

        # closed-form solution return-mapping for D-P
        c   = cmp.c₀+cmp.Hp*mpts.s.ϵpII[2,p]
        if c<cmp.cᵣ
            c = cmp.cᵣ
        end
        P,τ0,τII = σTr(σᵢ)
        η,ηB,ξ   = materialParam(cmp.ϕ₀,ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T2(1.0)+η^2)-η,τII-τP-(sqrt(T2(1.0)+η^2))*(P-σm)  
        if fs>T2(0.0) && P<σm || h>T2(0.0) && P≥σm 
            Δλ             = fs/(cmp.Gc+cmp.Kc*η*ηB)
            mpts.s.Δλ[p]     = Δλ
            Pn,τn          = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            mpts.s.σᵢ[p]   = CauchyStress(-Pn, τ0 .* (τn / τII))
            mpts.s.ϵpII[1,p]+= Δλ*sqrt(T2(1/3)+T2(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ             = (P-σm)/cmp.Kc
            mpts.s.Δλ[p]     = Δλ
            Pn             = σm-P
            mpts.s.σᵢ[p]   = CauchyStress(-Pn, τ0 .* (T2(0.0) / τII))
            mpts.s.ϵpII[1,p]+= sqrt(T2(2.0))*Δλ/T2(3.0)
        end
    end
end