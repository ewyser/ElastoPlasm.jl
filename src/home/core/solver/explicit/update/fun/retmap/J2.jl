# NOTE ON NAMING: these two `_yield_normal` methods used to be called `get_J2`, which
# collided semantically with `get_J2(::AbstractStress)` in `tensor.jl` — same
# name, different return value (this one returns `(‖ξ‖, n̂)`, the yield-surface normal
# and the norm of the deviator; that one returns the scalar invariant J₂). They were
# deliberately NOT merged during the tensor port: unifying them would have meant either
# silently changing what one call site gets back, or reformulating the return mapping in
# terms of the stored (p,dev) split, which is not canonically trace-free (see the
# `AbstractTensor` docstring). Renamed instead, so the collision is gone and each name
# means one thing.
@inline function _yield_normal(σ0::SVector{3,T}) where {T}
    P  = (σ0[1]+σ0[2])/T(2.0)
    ξ  = σ0 .- SVector{3,T}(P,P,zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+T(2.0)*ξ[3]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function _yield_normal(σ0::SVector{6,T}) where {T}
    P  = (σ0[1]+σ0[2]+σ0[3])/T(3.0)
    ξ  = σ0 .- SVector{6,T}(P,P,P,zero(T),zero(T),zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+ξ[3]^2+T(2.0)*ξ[4]^2+T(2.0)*ξ[5]^2+T(2.0)*ξ[6]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function yield_J2(σ,κ::T) where {T}
    ξn,n = _yield_normal(σ)
    f    = ξn-κ
    return f,n
end

"""
    _vonmises_return_map(σᵢ::SVector, ϵpII::MVector{2,T}, cmp, ftol, ηmax) -> (Δλ, σᵢ)

Closed-form-initialized, CPA (closest-point-projection)-refined J2/von Mises return
mapping, in Voigt space — shared by both `von_mises` methods below, mirroring `DP.jl`'s
`_druckerprager_return_map`. Returns `σᵢ` unchanged (elastic) and `Δλ=0` unless
yielding, in which case it's the CPA-converged deviator and the accumulated plastic
multiplier. On yield, `ϵpII[1]` is set to the accumulated hardening variable `γ0`
(seeded from `ϵpII[2]`, the non-local-regularized value — NOT incremented from
`ϵpII[1]`, unlike `DP.jl`'s convention); `ϵpII[2]` itself is never written here.
"""
@inline function _vonmises_return_map(σᵢ::SVector{N,T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel,ftol::T,ηmax::Int) where {N,T}
    κ     = max(cmp.cᵣ,cmp.c₀+cmp.Hp*ϵpII[2])
    f,n   = yield_J2(σᵢ,κ)
    Δλacc = T(0.0)
    σout  = σᵢ
    if f>T(0.0)
        γ0,σ0,ηit = ϵpII[2],σᵢ,1
        while abs(f)>ftol && ηit<ηmax
            ∂f∂σ   = n
            Δλ     = f/(∂f∂σ'*cmp.Del*∂f∂σ)
            Δσ     = Δλ .* (cmp.Del*∂f∂σ)
            σ0     = σ0 .- Δσ
            γ0    += Δλ
            Δλacc += Δλ
            κ      = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
            f,n    = yield_J2(σ0,κ)
            ηit   += 1
        end
        σout    = σ0
        ϵpII[1] = γ0
    end
    return Δλacc,σout
end

"""
    von_mises(τᵢⱼ::KirchhoffStress, ϵᵢⱼ::LogarithmicStrain, ϵpII, cmp; ftol, ηmax) -> (ϵᵢⱼ, τᵢⱼ, Δλ, ϵpII)
    von_mises(σᵢⱼ::CauchyStress, ϵᵢⱼ::InfinitesimalStrain, ϵpII, cmp; ftol, ηmax) -> (ϵᵢⱼ, σᵢⱼ, Δλ, ϵpII)

J2/von Mises return mapping, dispatched on stress/strain type — mirrors `DP.jl`'s
`drucker_prager`: one method per `(stress,strain)` pair, sharing the CPA core
`_vonmises_return_map` rather than duplicating the iterative loop (Borja (1990); De
Souza Neto (2008)). The finite-strain method additionally rebuilds the strain on yield
via `LogarithmicStrain(cmp.Del\\σ0)` (`Del\\σ0` is already an engineering-Voigt strain
vector — see `tensor.jl`); the infinitesimal-strain method never touches strain,
matching `elast!`'s incremental infinitesimal-strain tracking.
"""
@inline function von_mises(τᵢⱼ::KirchhoffStress{S,T,L},ϵᵢⱼ::LogarithmicStrain{S,T,L},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel;ftol::Real=1e-9,ηmax::Int=20) where {S,T,L}
    Δλ,τᵢ = _vonmises_return_map(get_voigt(τᵢⱼ),ϵpII,cmp,T(ftol),ηmax)
    if Δλ > T(0.0)
        τᵢⱼ = KirchhoffStress(τᵢ)
        ϵᵢⱼ = LogarithmicStrain(cmp.Del\τᵢ)
    end
    return ϵᵢⱼ,τᵢⱼ,Δλ,ϵpII
end
@inline function von_mises(σᵢⱼ::CauchyStress{S,T,L},ϵᵢⱼ::InfinitesimalStrain{S,T,L},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel;ftol::Real=1e-9,ηmax::Int=20) where {S,T,L}
    Δλ,σᵢ = _vonmises_return_map(get_voigt(σᵢⱼ),ϵpII,cmp,T(ftol),ηmax)
    if Δλ > T(0.0)
        σᵢⱼ = CauchyStress(σᵢ)
    end
    return ϵᵢⱼ,σᵢⱼ,Δλ,ϵpII
end

"""
    retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {CM<:VonMises, ST<:LogarithmicStrain}
    retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {CM<:VonMises, ST<:InfinitesimalStrain}

J2/von Mises plastic corrector — contributes the `CM<:VonMises` methods to the shared
`retmap` kernel name (see `DP.jl`'s `retmap` docstring for the `CM`+`ST` dispatch
rationale). `ftol`/`ηmax` control the CPA loop's convergence — they default inside
`von_mises` itself (a plain function, not `@kernel`-decorated) rather than being
kernel-level keyword arguments, since KernelAbstractions.jl doesn't support keyword
arguments on `@kernel` functions for GPU backends (CPU-only tolerates them). The
kernels below just call `von_mises` with no `ftol`/`ηmax` at all, relying on its
defaults — matching the values `finite_J2`/`infinitesimal_J2` always ran with in
practice, since their call site never passed these explicitly either.
"""
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM<:VonMises,TM,TV,TS,ST<:LogarithmicStrain}
    p = @index(Global)
    if p≤mpts.nmp
        ϵᵢⱼ,τᵢⱼ,Δλ,ϵpII = von_mises(mpts.s.τᵢⱼ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.ϵᵢⱼ[p]  = ϵᵢⱼ
            mpts.s.τᵢⱼ[p]  = τᵢⱼ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        else
            mpts.s.Δλ[p] = T2(0.0)
        end
    end
end
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM<:VonMises,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p≤mpts.nmp
        ϵᵢⱼ,σᵢⱼ,Δλ,ϵpII = von_mises(mpts.s.σᵢⱼ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.σᵢⱼ[p]  = σᵢⱼ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        else
            mpts.s.Δλ[p] = T2(0.0)
        end
    end
end
