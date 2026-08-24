export PerfectlyElastic,DruckerPrager

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

Both `finite_DP`/`infinitesimal_DP` (`retmap/DP.jl`) and `finite_J2`/`infinitesimal_J2`
(`retmap/J2.jl`) read the same field set off this struct (`J2` simply never reads
`ϕ₀`), so `setup_cmp` always builds `DruckerPrager` regardless of
`plast.constitutive ∈ {"DP","J2"}`.

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
