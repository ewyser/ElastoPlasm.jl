# `_yield_normal` returns `(‖ξ‖, n̂)` for a stress deviator `ξ` (Voigt, shear stored once):
# ‖ξ‖ = √(2J₂) and the unit deviator n̂ = ξ/‖ξ‖. Not to be confused with `get_J2` (`stress.jl`),
# the scalar invariant.
@inline function _yield_normal(ξ::SVector{6,T}) where {T}
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+ξ[3]^2+T(2.0)*ξ[4]^2+T(2.0)*ξ[5]^2+T(2.0)*ξ[6]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    return ξn,ξ./ξn
end
@inline function yield_J2(ξ,κ::T) where {T}
    ξn,n = _yield_normal(ξ)
    return ξn-κ,n
end

"""
    _vonmises_return_map(σ::AbstractStress, ϵpII::MVector{2,T}, cmp, ftol, ηmax) -> (Δλ, dev)

CPA (closest-point-projection) J2/von Mises return mapping on the stored deviator — shared by
both `von_mises` methods below, mirroring `DP.jl`'s `_druckerprager_return_map`. The flow
direction `n̂` is deviatoric, so `Del·∂f/∂σ = 2Gc·n̂` and `∂f/∂σ·Del·∂f/∂σ = 2Gc`: each step
moves only the deviator, radially, and the pressure is unchanged. Returns `σ.dev` (elastic) and
`Δλ=0` unless yielding, in which case it's the converged deviator and the accumulated plastic
multiplier. On yield, `ϵpII[1]` is set to the accumulated hardening variable `γ0` (seeded from
`ϵpII[2]`, the non-local-regularized value — NOT incremented from `ϵpII[1]`, unlike `DP.jl`'s
convention); `ϵpII[2]` itself is never written here.
"""
@inline function _vonmises_return_map(σ::AbstractStress{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel,ftol::T,ηmax::Int) where {T}
    κ     = max(cmp.cᵣ,cmp.c₀+cmp.Hp*ϵpII[2])
    dev   = σ.dev
    f,n   = yield_J2(dev,κ)
    Δλacc = T(0.0)
    if f>T(0.0)
        γ0,ηit = ϵpII[2],1
        while abs(f)>ftol && ηit<ηmax
            Δλ     = f/(T(2.0)*cmp.Gc)
            dev    = dev .- (T(2.0)*cmp.Gc*Δλ) .* n
            γ0    += Δλ
            Δλacc += Δλ
            κ      = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
            f,n    = yield_J2(dev,κ)
            ηit   += 1
        end
        ϵpII[1] = γ0
    end
    return Δλacc,dev
end

"""
    von_mises(τᵢ::KirchhoffStress, ϵᵢ::LogarithmicStrain, ϵpII, cmp; ftol, ηmax) -> (ϵᵢ, τᵢ, Δλ, ϵpII)
    von_mises(σᵢ::CauchyStress, ϵᵢ::InfinitesimalStrain, ϵpII, cmp; ftol, ηmax) -> (ϵᵢ, σᵢ, Δλ, ϵpII)

J2/von Mises return mapping, dispatched on stress/strain type — mirrors `DP.jl`'s
`drucker_prager`: one method per `(stress,strain)` pair, sharing the CPA core
`_vonmises_return_map` rather than duplicating the iterative loop (Borja (1990); De
Souza Neto (2008)). The finite-strain method additionally rebuilds the elastic strain on yield
with `_elastic_strain` (`stress.jl`); the infinitesimal-strain method never touches strain,
matching `elast!`'s incremental infinitesimal-strain tracking.
"""
@inline function von_mises(τᵢ::KirchhoffStress{T},ϵᵢ::LogarithmicStrain{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel;ftol::Real=1e-9,ηmax::Int=20) where {T}
    Δλ,dev = _vonmises_return_map(τᵢ,ϵpII,cmp,T(ftol),ηmax)
    if Δλ > T(0.0)
        τᵢ = KirchhoffStress(τᵢ.p,dev)
        ϵᵢ = _elastic_strain(τᵢ,cmp)
    end
    return ϵᵢ,τᵢ,Δλ,ϵpII
end
@inline function von_mises(σᵢ::CauchyStress{T},ϵᵢ::InfinitesimalStrain{T},ϵpII::MVector{2,T},cmp::AbstractConstitutiveModel;ftol::Real=1e-9,ηmax::Int=20) where {T}
    Δλ,dev = _vonmises_return_map(σᵢ,ϵpII,cmp,T(ftol),ηmax)
    if Δλ > T(0.0)
        σᵢ = CauchyStress(σᵢ.p,dev)
    end
    return ϵᵢ,σᵢ,Δλ,ϵpII
end

"""
    retmap(mpts::Point{T1,T2,D,CM,ST}) where {CM<:VonMises, ST<:LogarithmicStrain}
    retmap(mpts::Point{T1,T2,D,CM,ST}) where {CM<:VonMises, ST<:InfinitesimalStrain}

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
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM<:VonMises,ST<:LogarithmicStrain}
    p = @index(Global)
    if p≤mpts.nmp
        ϵᵢ,τᵢ,Δλ,ϵpII = von_mises(mpts.s.τᵢ[p],mpts.s.ϵᵢ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.ϵᵢ[p]  = ϵᵢ
            mpts.s.τᵢ[p]  = τᵢ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        else
            mpts.s.Δλ[p] = T2(0.0)
        end
    end
end
@kernel inbounds = true function retmap(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM<:VonMises,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p≤mpts.nmp
        ϵᵢ,σᵢ,Δλ,ϵpII = von_mises(mpts.s.σᵢ[p],mpts.s.ϵᵢ[p],MVector{2,T2}(mpts.s.ϵpII[p]),mpts.s.cmp[p])
        if Δλ > T2(0.0)
            mpts.s.σᵢ[p]  = σᵢ
            mpts.s.Δλ[p]   = Δλ
            mpts.s.ϵpII[p] = SVector{2,T2}(ϵpII)
        else
            mpts.s.Δλ[p] = T2(0.0)
        end
    end
end
