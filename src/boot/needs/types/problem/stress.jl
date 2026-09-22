# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stress tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractStress, CauchyStress, KirchhoffStress
export get_tensor, get_voigt, get_J2, get_τII

"""
    AbstractStress{S,T,L}

Second-order stress tensor stored on `PointSolidPhase` (`σᵢⱼ`/`σn`/`τᵢⱼ`) as
`p::T` (pressure, **positive in compression**) + `dev::SMatrix{S,S,T,L}` (deviatoric
part), so the full tensor is `dev - p·I`. `strain.jl`'s `AbstractStrain` follows the
same `(vol, dev)`-split convention for strains (see its docstring for the full note
on `dev` not being canonically trace-free and always reading via `get_voigt`/
`get_tensor`).
"""
abstract type AbstractStress{S,T,L} end

"""
    get_tensor(stress::AbstractStress{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full stress tensor `dev - p·I` (`p` positive in compression).
"""
@inline function get_tensor(stress::AbstractStress{S,T,L}) where {S,T,L}
    return stress.dev - stress.p * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(stress::AbstractStress) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt stress vector `(σxx,σyy,σxy)` / `(σxx,σyy,σzz,σyz,σxz,σxy)` — no shear doubling
(unlike strain).
"""
@inline function get_voigt(stress::AbstractStress{2,T,L}) where {T,L}
    return SVector{3,T}(
        stress.dev[1,1] - stress.p,
        stress.dev[2,2] - stress.p,
        stress.dev[1,2],
    )
end
@inline function get_voigt(stress::AbstractStress{3,T,L}) where {T,L}
    return SVector{6,T}(
        stress.dev[1,1] - stress.p,
        stress.dev[2,2] - stress.p,
        stress.dev[3,3] - stress.p,
        stress.dev[2,3],
        stress.dev[1,3],
        stress.dev[1,2],
    )
end

"""
    get_J2(stress::AbstractStress) -> T

Second deviatoric stress invariant `J₂ = ½ sᵢⱼsᵢⱼ`, re-derived from `get_voigt` (see
`AbstractStrain`/`AbstractStress`'s shared-convention notes) rather than trusting the
stored `dev`.
"""
@inline function get_J2(stress::AbstractStress{2,T,L}) where {T,L}
    σ = get_voigt(stress)
    P = (σ[1] + σ[2]) / T(2.0)
    return T(0.5) * ((σ[1] - P)^2 + (σ[2] - P)^2) + σ[3]^2
end
@inline function get_J2(stress::AbstractStress{3,T,L}) where {T,L}
    σ = get_voigt(stress)
    P = (σ[1] + σ[2] + σ[3]) / T(3.0)
    return T(0.5) * ((σ[1] - P)^2 + (σ[2] - P)^2 + (σ[3] - P)^2) + σ[4]^2 + σ[5]^2 + σ[6]^2
end

"""
    get_τII(stress::AbstractStress) -> T

Second invariant `τII = √J₂` of the deviatoric stress.
"""
@inline function get_τII(stress::AbstractStress)
    return sqrt(get_J2(stress))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CauchyStress
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    CauchyStress{S,T,L} <: AbstractStress{S,T,L}

Cauchy (true) stress. Stored by `mpts.s.σᵢⱼ`/`mpts.s.σn`.
"""
struct CauchyStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    function CauchyStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
    function CauchyStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
end
@adapt_struct CauchyStress

# dev given as a Voigt deviator only (no pressure part)
@inline function CauchyStress(p::T, dev::SVector{3,T}) where {T}
    return CauchyStress(p, SMatrix{2,2,T,4}(
        dev[1], dev[3],
        dev[3], dev[2],
    ))
end
# 3D: SMatrix fills column-major — dev = [dev1 dev4 dev5; dev4 dev2 dev6; dev5 dev6 dev3].
@inline function CauchyStress(p::T, dev::SVector{6,T}) where {T}
    return CauchyStress(p, SMatrix{3,3,T,9}(
        dev[1], dev[4], dev[5],
        dev[4], dev[2], dev[6],
        dev[5], dev[6], dev[3],
    ))
end
# full Voigt stress vector (not a deviator) — 2D, p = -tr(σ)/2
@inline function CauchyStress(σ::SVector{3,T}) where {T}
    p = -(σ[1] + σ[2]) / T(2.0)
    return CauchyStress(p, SMatrix{2,2,T,4}(
        σ[1] + p, σ[3],
        σ[3],     σ[2] + p,
    ))
end
# full Voigt stress vector (not a deviator) — 3D, p = -tr(σ)/3
@inline function CauchyStress(σ::SVector{6,T}) where {T}
    p = -(σ[1] + σ[2] + σ[3]) / T(3.0)
    return CauchyStress(p, SMatrix{3,3,T,9}(
        σ[1] + p, σ[6],     σ[5],
        σ[6],     σ[2] + p, σ[4],
        σ[5],     σ[4],     σ[3] + p,
    ))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# KirchhoffStress
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}

Kirchhoff stress `τ = J·σ`. Stored by `mpts.s.τᵢⱼ`; `transform` divides by `J` to get
`mpts.s.σᵢⱼ`.
"""
struct KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    function KirchhoffStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
    function KirchhoffStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
end
@adapt_struct KirchhoffStress

# dev given as a Voigt deviator only (no pressure part)
@inline function KirchhoffStress(p::T, dev::SVector{3,T}) where {T}
    return KirchhoffStress(p, SMatrix{2,2,T,4}(
        dev[1], dev[3],
        dev[3], dev[2],
    ))
end
# 3D: SMatrix fills column-major — dev = [dev1 dev4 dev5; dev4 dev2 dev6; dev5 dev6 dev3].
@inline function KirchhoffStress(p::T, dev::SVector{6,T}) where {T}
    return KirchhoffStress(p, SMatrix{3,3,T,9}(
        dev[1], dev[4], dev[5],
        dev[4], dev[2], dev[6],
        dev[5], dev[6], dev[3],
    ))
end
# full Voigt stress vector (not a deviator) — 2D, p = -tr(τ)/2
@inline function KirchhoffStress(τ::SVector{3,T}) where {T}
    p = -(τ[1] + τ[2]) / T(2.0)
    return KirchhoffStress(p, SMatrix{2,2,T,4}(
        τ[1] + p, τ[3],
        τ[3],     τ[2] + p,
    ))
end
# full Voigt stress vector (not a deviator) — 3D, p = -tr(τ)/3
@inline function KirchhoffStress(τ::SVector{6,T}) where {T}
    p = -(τ[1] + τ[2] + τ[3]) / T(3.0)
    return KirchhoffStress(p, SMatrix{3,3,T,9}(
        τ[1] + p, τ[6],     τ[5],
        τ[6],     τ[2] + p, τ[4],
        τ[5],     τ[4],     τ[3] + p,
    ))
end

"""
    _trial_elastic_stress(strain::LogarithmicStrain, cmp::AbstractConstitutiveModel) -> KirchhoffStress

Isotropic linear-elastic trial Kirchhoff stress: `p = -Kc·ϵvol` (positive in
compression), `dev = 2·Gc·ϵdev`. Called from `elast.jl`'s `SM<:HenckySolid` kernel
method — which law runs is picked by `Point`'s own `SM` type parameter (see
`AbstractSolid` in `lagrangian.jl`), not by branching here.
"""
@inline function _trial_elastic_stress(strain::LogarithmicStrain{S,T,L}, cmp::AbstractConstitutiveModel{T}) where {S,T,L}
    P   =         - cmp.Kc * strain.vol
    dev =  T(2.0) * cmp.Gc * strain.dev
    return KirchhoffStress(P, dev)
end

"""
    _trial_elastic_stress_improved(strain::LogarithmicStrain, cmp::AbstractConstitutiveModel, n₀) -> KirchhoffStress

Porosity-weighted "improved Hencky" trial Kirchhoff stress (Pretti, Coombs, Augarde,
Marchena Puigvert, Reyna Gutiérrez, *Mechanics of Materials* 192 (2024) 104958, Eqs.
33/34 — elastic part only; no fluid/Terzaghi term, no plastic hardening term). Called
from `elast.jl`'s `SM<:ImprovedHenckySolid` kernel method. `n` is recomputed here from
the elastic volumetric strain via their Eq. (23), `n = 1 - (1-n₀)/exp(ϵᵥᵉ)`,
deliberately **not** read from `Point.n` (which tracks total, not purely-elastic,
deformation — see `planned-improvements.md` for why the two porosity notions are kept
separate). The deviatoric part is unaffected by porosity (paper's own §3, citing
Zytynski et al. 1978: a variable-K/constant-G material is non-hyperelastic otherwise).
Their `p'` is positive in tension; this repo's `p` is positive in compression, hence
the sign flip.
"""
@inline function _trial_elastic_stress_improved(strain::LogarithmicStrain{S,T,L}, cmp::AbstractConstitutiveModel{T}, n₀::T) where {S,T,L}
    ϵᵥᵉ = strain.vol
    n   = T(1.0) - (T(1.0) - n₀) / exp(ϵᵥᵉ)
    P   = -(cmp.Kc * ϵᵥᵉ / n) * (T(1.0) + (ϵᵥᵉ / (T(2.0) * n)) * (n - T(1.0)))
    dev =  T(2.0) * cmp.Gc * strain.dev
    return KirchhoffStress(P, dev)
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Shared
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Base.zero(::Type{CauchyStress{S,T,L}})    where {S,T,L} = CauchyStress(   zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{KirchhoffStress{S,T,L}}) where {S,T,L} = KirchhoffStress(zero(T), zero(SMatrix{S,S,T,L}))

"""
    CauchyStress(τ::KirchhoffStress, J⁻¹) -> CauchyStress
    KirchhoffStress(σ::CauchyStress) -> KirchhoffStress

Cross-type conversions (`τ ↔ σ`), used by `transform` and the u-P `dynamic_relaxation`
path. Defined outside both structs since `CauchyStress`'s version needs `KirchhoffStress`,
which isn't defined until later in the file.
"""
@inline function CauchyStress(τ::KirchhoffStress{S,T,L}, J⁻¹::T) where {S,T,L}
    return CauchyStress(τ.p * J⁻¹, τ.dev * J⁻¹)
end
@inline function KirchhoffStress(σ::CauchyStress{S,T,L}) where {S,T,L}
    return KirchhoffStress(σ.p, σ.dev)
end
