export PerfectlyElastic,DruckerPrager,VonMises

"""
    PerfectlyElastic{T2,D} <: AbstractConstitutiveModel{T2,D}

Purely elastic per-particle constitutive constants: shear modulus `Gc` and bulk
modulus `Kc`. Not currently constructed anywhere in `setup_cmp` (no code path in this
repo selects a non-plastic constitutive model, mirroring the already-unwired
`plast.constitutive ∈ {"MC","camC"}` branches in `update.jl`'s `init_update`) — defined
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
Drucker-Prager plastic model: shear modulus `Gc`, bulk modulus `Kc`, elastic stiffness
matrix `Del`, softening modulus `Hp`, and initial/residual cohesion + friction angle
`c₀`,`cᵣ`,`ϕ₀`. One instance per material point (`mpts.s.cmp::Vector{DruckerPrager}`),
built by `setup_cmp`. Evolving plastic state (`Δλ`,`ϵpII`,`ϵpV`) stays on
`PointSolidPhase` as flat mutable vectors, not in this bundle.

Built by `setup_cmp` when `plast.constitutive == "DP"`. `VonMises` (below) is the
sibling built for `plast.constitutive == "VM"` — the same field set minus `ϕ₀`, since
J2/von Mises yield has no pressure dependence. Both `retmap/DP.jl` and `retmap/J2.jl`
contribute methods to the shared `retmap` kernel, dispatched on `Point`'s `CM`
(`DruckerPrager`/`VonMises`) type parameter together with `ST`
(`LogarithmicStrain`/`InfinitesimalStrain`).

Deviates from the plan's literal `Del::SMatrix{D,D,T2}` field spec: `Del` is the Voigt-
notation elastic stiffness matrix, sized `nstr×nstr` (3×3 in 2D, 6×6 in 3D — see
`get_elastic_stiffness`), not `D×D` (2×2/3×3) — using `D` directly would silently build
the wrong-sized matrix. `NSTR`/`L` carry the actual Voigt size instead, with `D` kept
as the second parameter purely to satisfy `AbstractConstitutiveModel{T2,D}`.
"""
struct DruckerPrager{T2,D,NSTR,L} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Del::SMatrix{NSTR,NSTR,T2,L}
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
    ϕ₀ ::T2
end
@adapt_struct DruckerPrager

"""
    VonMises{T2,D,NSTR,L} <: AbstractConstitutiveModel{T2,D}

Per-particle static constitutive constants for the J2/von Mises plastic model: shear
modulus `Gc`, bulk modulus `Kc`, elastic stiffness matrix `Del`, softening modulus `Hp`,
and initial/residual cohesion `c₀`,`cᵣ`. Mirrors `DruckerPrager` exactly, minus the
friction angle `ϕ₀` — J2/von Mises yield has no pressure dependence, so there is no
friction angle to store. One instance per material point
(`mpts.s.cmp::Vector{VonMises}`), built by `setup_cmp` when `plast.constitutive=="VM"`.
"""
struct VonMises{T2,D,NSTR,L} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Del::SMatrix{NSTR,NSTR,T2,L}
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
end
@adapt_struct VonMises
