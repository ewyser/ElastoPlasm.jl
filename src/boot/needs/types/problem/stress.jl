# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stress tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractStress, CauchyStress, KirchhoffStress
export get_tensor, get_J2, get_τII

"""
    AbstractStress{T}

Second-order stress tensor stored on `PointSolidPhase` (`σᵢ`/`σn`/`τᵢ`) as its pressure
`p::T = -tr/3` (positive in compression) and its trace-free deviator `dev::SVector{6,T}`, a 3D
Voigt vector `(xx, yy, zz, yz, xz, xy)` in every dimension, like `AbstractStrain` (`strain.jl`) but
without shear doubling: in 2D plane strain `σzz` is carried.

As for strain, tensors meet the split in one place each way: `get_tensor(σ) = dev - p·I` (symmetric
3×3; `get_tensor(σ, Val(D))` is its in-plane `D×D` block, for force assembly) and the 3×3
constructor. Other readers: `get_J2`/`get_τII`.
Code that builds `(p, dev)` directly (the elastic laws, the return maps) passes a trace-free `dev`.
"""
abstract type AbstractStress{T} end

@inline function get_tensor(t::AbstractStress{T}) where {T}
    d, p = t.dev, t.p
    return SMatrix{3,3,T,9}(d[1] - p, d[6], d[5], d[6], d[2] - p, d[4], d[5], d[4], d[3] - p)
end
@inline get_tensor(t::AbstractStress, ::Val{D}) where {D} = get_tensor(t)[SOneTo(D), SOneTo(D)]

"""
    get_J2(stress::AbstractStress) -> T

Second deviatoric stress invariant `J₂ = ½ sᵢⱼsᵢⱼ` of the full 3D stress.
"""
@inline function get_J2(t::AbstractStress{T}) where {T}
    s = t.dev
    return T(0.5) * (s[1]^2 + s[2]^2 + s[3]^2) + s[4]^2 + s[5]^2 + s[6]^2
end

"""
    get_τII(stress::AbstractStress) -> T

Second invariant `τII = √J₂` of the deviatoric stress.
"""
@inline get_τII(t::AbstractStress) = sqrt(get_J2(t))

"""
    CauchyStress{T} <: AbstractStress{T}

Cauchy (true) stress. Stored by `mpts.s.σᵢ`/`mpts.s.σn`.
"""
struct CauchyStress{T} <: AbstractStress{T}
    p  ::T
    dev::SVector{6,T}
end
@adapt_struct CauchyStress

"""
    KirchhoffStress{T} <: AbstractStress{T}

Kirchhoff stress `τ = J·σ`. Stored by `mpts.s.τᵢ`; `transform` divides by `J` to get
`mpts.s.σᵢ`.
"""
struct KirchhoffStress{T} <: AbstractStress{T}
    p  ::T
    dev::SVector{6,T}
end
@adapt_struct KirchhoffStress

# from a symmetric 3×3 stress tensor: `p = -tr/3`, `dev = M + p·I`
@inline function CauchyStress(M::SMatrix{3,3,T,9}) where {T}
    p = -(M[1,1] + M[2,2] + M[3,3]) / T(3.0)
    return CauchyStress(p, SVector{6,T}(M[1,1] + p, M[2,2] + p, M[3,3] + p, M[2,3], M[1,3], M[1,2]))
end
# from a D×D in-plane block with zero out-of-plane components (see the strain counterpart)
@inline CauchyStress(M::SMatrix{D,D,T}) where {D,T} = CauchyStress(_block3(M))

"""
    CauchyStress(τ::KirchhoffStress, J⁻¹) -> CauchyStress

`σ = τ/J`, used by `transform` and the dynamic-relaxation `update`.
"""
@inline CauchyStress(τ::KirchhoffStress{T}, J⁻¹::T) where {T} = CauchyStress(τ.p * J⁻¹, τ.dev * J⁻¹)

# maps an engineering-shear strain deviator to `σ_dev / G`: 2 on the normals, 1 on the shears
@inline _shear_weights(::Type{T}) where {T} = SVector{6,T}(2, 2, 2, 1, 1, 1)

"""
    _trial_elastic_stress(strain::LogarithmicStrain, cmp::AbstractConstitutiveModel) -> KirchhoffStress

Isotropic linear-elastic trial Kirchhoff stress `τ = Kc·tr(ϵ)·I + 2·Gc·dev(ϵ)`: pressure
`p = -Kc·ϵᵥ`, deviator `2Gc` on the normals and `Gc·γ` on the (engineering) shears. Called from `elast.jl`'s
`EL<:HenckySolid` kernel method — which law runs is picked by `Point`'s own `EL` type parameter
(see `AbstractElasticLaw` in `lagrangian.jl`), not by branching here.
"""
@inline function _trial_elastic_stress(strain::LogarithmicStrain{T}, cmp::AbstractConstitutiveModel{T}) where {T}
    return KirchhoffStress(-cmp.Kc * strain.vol, cmp.Gc * strain.dev .* _shear_weights(T))
end

"""
    _elastic_strain(τ::KirchhoffStress, cmp::AbstractConstitutiveModel) -> LogarithmicStrain

Inverse of `_trial_elastic_stress`: `vol = -p/Kc`, and the deviator divided by `2Gc` on the
normals and by `Gc` on the (engineering) shears. Used by the return maps to rebuild the elastic
strain from the returned stress, under improved Hencky too (where it is the linear-Hencky
approximation).
"""
@inline _elastic_strain(τ::KirchhoffStress{T}, cmp::AbstractConstitutiveModel{T}) where {T} =
    LogarithmicStrain(-τ.p / cmp.Kc, τ.dev ./ (cmp.Gc * _shear_weights(T)))

"""
    _trial_elastic_stress_improved(strain::LogarithmicStrain, cmp::AbstractConstitutiveModel, n₀) -> KirchhoffStress

Porosity-weighted "improved Hencky" trial Kirchhoff stress (Pretti, Coombs, Augarde,
Marchena Puigvert, Reyna Gutiérrez, *Mechanics of Materials* 192 (2024) 104958, Eqs.
33/34 — elastic part only; no fluid/Terzaghi term, no plastic hardening term). Called
from `elast.jl`'s `EL<:ImprovedHenckySolid` kernel method. `n` is recomputed here from
the elastic volumetric strain via their Eq. (23), `n = 1 - (1-n₀)/exp(ϵᵥᵉ)`. The
paper's `n` is the total-`J` porosity (Eqs. 22-23, 45) — this equals it only while
plastic flow is isochoric, which holds today (`DP.jl` runs with dilatancy `ψ=0`, J2 is
isochoric); a dilatant flow rule would need the total porosity here instead. The deviatoric part is unaffected by porosity (paper's own §3, citing
Zytynski et al. 1978: a variable-K/constant-G material is non-hyperelastic otherwise).
Their `p'` is positive in tension; this repo's `p` is positive in compression, hence
the sign flip.
"""
@inline function _trial_elastic_stress_improved(strain::LogarithmicStrain{T}, cmp::AbstractConstitutiveModel{T}, n₀::T) where {T}
    ϵᵥᵉ  = strain.vol
    n    = T(1.0) - (T(1.0) - n₀) / exp(ϵᵥᵉ)
    P    = -(cmp.Kc * ϵᵥᵉ / n) * (T(1.0) + (ϵᵥᵉ / (T(2.0) * n)) * (n - T(1.0)))
    return KirchhoffStress(P, cmp.Gc * strain.dev .* _shear_weights(T))
end

Base.zero(::Type{CauchyStress{T}})    where {T} = CauchyStress(zero(T), zero(SVector{6,T}))
Base.zero(::Type{KirchhoffStress{T}}) where {T} = KirchhoffStress(zero(T), zero(SVector{6,T}))
