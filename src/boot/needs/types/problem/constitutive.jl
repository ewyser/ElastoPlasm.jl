abstract type AbstractConstitutiveModel{T2, D} end

export PerfectlyElastic,DruckerPrager,VonMises

"""
    PerfectlyElastic{T2,D} <: AbstractConstitutiveModel{T2,D}

Purely elastic per-particle constitutive constants: shear modulus `Gc` and bulk
modulus `Kc`. Not currently constructed anywhere in `setup_cmp` (no code path in this
repo selects a non-plastic constitutive model, mirroring the already-unwired
`material.plastic ∈ {"MC","camC"}` branches in `update.jl`'s `init_update`) — defined
here as forward-looking scaffolding for whenever a genuine elastic-only workflow exists.
"""
struct PerfectlyElastic{T2,D} <: AbstractConstitutiveModel{T2,D}
    Gc::T2
    Kc::T2
end
@adapt_struct PerfectlyElastic

"""
    DruckerPrager{T2,D} <: AbstractConstitutiveModel{T2,D}

Per-particle static (never mutated after init) constitutive constants for the
Drucker-Prager plastic model: shear modulus `Gc`, bulk modulus `Kc`, softening modulus `Hp`, and initial/residual cohesion + friction angle
`c₀`,`cᵣ`,`ϕ₀`. One instance per material point (`mpts.s.cmp::Vector{DruckerPrager}`),
built by `setup_cmp`. Evolving plastic state (`Δλ`,`ϵpII`,`ϵpV`) stays on
`PointSolidPhase` as flat mutable vectors, not in this bundle.

Built by `setup_cmp` when `material.plastic == "DP"`. `VonMises` (below) is the
sibling built for `material.plastic == "VM"` — the same field set minus `ϕ₀`, since
J2/von Mises yield has no pressure dependence. Both `retmap/DP.jl` and `retmap/J2.jl`
contribute methods to the shared `retmap` kernel, dispatched on `Point`'s `CM`
(`DruckerPrager`/`VonMises`) type parameter together with `ST`
(`LogarithmicStrain`/`InfinitesimalStrain`).

Which elastic trial-stress law (`Hencky`/`ImprovedHencky`) applies is **not** encoded
here — it lives on `Point`'s own `EL<:AbstractElasticLaw` type parameter instead (see
`lagrangian.jl`), since it's a solid-phase/kinematics concern, not a plastic-model
constant.

The elastic stiffness is not stored: `get_Del(cmp)` (`stiffness.jl`) builds the 6×6 Voigt `Del`
from `Kc`/`Gc` where needed. `D` is the problem dimension, which Drucker–Prager uses to pick its 2D
(plane-strain) parameter fit.
"""
struct DruckerPrager{T2,D} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
    ϕ₀ ::T2
end
@adapt_struct DruckerPrager

"""
    VonMises{T2,D} <: AbstractConstitutiveModel{T2,D}

Per-particle static constitutive constants for the J2/von Mises plastic model: shear
modulus `Gc`, bulk modulus `Kc`, softening modulus `Hp`,
and initial/residual cohesion `c₀`,`cᵣ`. Mirrors `DruckerPrager` exactly, minus the
friction angle `ϕ₀` — J2/von Mises yield has no pressure dependence, so there is no
friction angle to store. One instance per material point
(`mpts.s.cmp::Vector{VonMises}`), built by `setup_cmp` when `material.plastic=="VM"`.
"""
struct VonMises{T2,D} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
end
@adapt_struct VonMises
