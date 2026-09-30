# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The elastic step on the strain and stress objects (tensor.jl): the Hencky push-forward
# (trial elastic strain), the linear and improved Hencky laws, the inverse used by the return maps,
# and the isotropic stiffness in Voigt form
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export StressStrainStiffness, get_Del

"""
    _trial_elastic_strain(ΔFᵢⱼ::SMatrix{D,D}, ϵᵢ::LogarithmicStrain) -> LogarithmicStrain

Push the stored logarithmic strain forward through the incremental deformation gradient:
`bᵉ = ΔF·exp(2ϵ)·ΔFᵀ`, `ϵ = ½ ln bᵉ`. Only the in-plane `D×D` block is pushed forward: in plane
strain (2D) and uniaxial strain (1D) `ΔF` has no out-of-plane stretch or shear (`ΔF_zz = 1`), so
`bᵉ` stays block-diagonal and the out-of-plane components are unchanged (`b_zz = exp(2ϵzz)`, hence
`ϵzz` is kept as is); `LogarithmicStrain(ϵᵢⱼ, ϵᵢ)` builds the new strain from the pushed-forward
block, carrying the out-of-plane components (the carry constructor, `tensor.jl`).
"""
@inline function _trial_elastic_strain(ΔFᵢⱼ::SMatrix{D,D,T}, ϵᵢ::LogarithmicStrain{T}) where {D,T}
    # get the in-plane block of the previously-converged logarithmic strain
    ϵᵢⱼ  = get_tensor(ϵᵢ, Val(D))
    λ, n = eigen(Symmetric(ϵᵢⱼ))
    # calculate left Cauchy-Green strain
    bᵢⱼ  = ΔFᵢⱼ * (n * diagm(exp.(T(2.0) * λ)) * n') * ΔFᵢⱼ'
    # update logarithmic strain tensor from the eigen-decomposition of bᵢⱼ
    λ, n = eigen(Symmetric(bᵢⱼ))
    ϵᵢⱼ  = T(0.5) * (n * diagm(log.(λ)) * n')
    return LogarithmicStrain(ϵᵢⱼ, ϵᵢ)
end

# maps an engineering-shear strain deviator to its tensor components (1 on the normals, ½ on the
# shears), so that `σ_dev = 2G·ϵ_dev` reads as the standard isotropic law
@inline _shear_weights(::Type{T}) where {T} = SVector{6,T}(1, 1, 1, 0.5, 0.5, 0.5)

"""
    _trial_elastic_stress(strain::LogarithmicStrain, cmp::AbstractConstitutiveModel) -> KirchhoffStress

Isotropic linear-elastic trial Kirchhoff stress `τ = Kc·tr(ϵ)·I + 2·Gc·dev(ϵ)`: pressure
`p = -Kc·ϵᵥ`, deviator `2Gc` on the normals and `Gc·γ` on the (engineering) shears. Called from `elast.jl`'s
`EL<:HenckySolid` kernel method — which law runs is picked by `Point`'s own `EL` type parameter
(see `AbstractElasticLaw` in `lagrangian.jl`), not by branching here.
"""
@inline function _trial_elastic_stress(strain::LogarithmicStrain{T}, cmp::AbstractConstitutiveModel{T}) where {T}
    return KirchhoffStress(-cmp.Kc * strain.vol, T(2.0) * cmp.Gc * strain.dev .* _shear_weights(T))
end

"""
    _elastic_strain(τ::KirchhoffStress, cmp::AbstractConstitutiveModel) -> LogarithmicStrain

Inverse of `_trial_elastic_stress`: `vol = -p/Kc`, and the deviator divided by `2Gc` on the
normals and by `Gc` on the (engineering) shears. Used by the return maps to rebuild the elastic
strain from the returned stress, under improved Hencky too (where it is the linear-Hencky
approximation).
"""
@inline _elastic_strain(τ::KirchhoffStress{T}, cmp::AbstractConstitutiveModel{T}) where {T} =
    LogarithmicStrain(-τ.p / cmp.Kc, τ.dev ./ (T(2.0) * cmp.Gc * _shear_weights(T)))

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
    return KirchhoffStress(P, T(2.0) * cmp.Gc * strain.dev .* _shear_weights(T))
end

# ------------------------------------------------ stiffness

"""
    StressStrainStiffness{T}()

Volumetric and deviatoric projections on a 6-component Voigt strain `(xx, yy, zz, γyz, γxz, γxy)`
with engineering shear (`γ = 2ϵ`): `vol·ϵ` is `tr(ϵ)` on each normal component, and `dev·ϵ` is the
deviatoric strain, with `½` on the shear diagonal so that `2G·dev` maps `γ` to the shear stress
`G·γ`. Any isotropic stiffness is `K·vol + 2G·dev`; see `get_Del`.
"""
struct StressStrainStiffness{T}
    vol::SMatrix{6,6,T,36}
    dev::SMatrix{6,6,T,36}
end
@inline function StressStrainStiffness{T}() where {T}
    a, b, h = T(2.0/3.0), T(-1.0/3.0), T(0.5)
    vol = SMatrix{6,6,T,36}(
        1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
    )
    dev = SMatrix{6,6,T,36}(
        a, b, b, 0, 0, 0,
        b, a, b, 0, 0, 0,
        b, b, a, 0, 0, 0,
        0, 0, 0, h, 0, 0,
        0, 0, 0, 0, h, 0,
        0, 0, 0, 0, 0, h,
    )
    return StressStrainStiffness{T}(vol, dev)
end

"""
    get_Del(cmp::AbstractConstitutiveModel) -> SMatrix{6,6}

Isotropic elastic stiffness `Del = Kc·vol + 2Gc·dev` of a constitutive model, built on the fly from
its moduli (acts on engineering-shear Voigt strain, gives Voigt stress). Not stored per particle.
The laws above do not use it: they apply `Kc`/`Gc` to the pressure/deviator split directly. It is
kept for code that needs the matrix itself, e.g. a consistent tangent.
"""
@inline function get_Del(cmp::AbstractConstitutiveModel{T}) where {T}
    P = StressStrainStiffness{T}()
    return cmp.Kc * P.vol + T(2.0) * cmp.Gc * P.dev
end
