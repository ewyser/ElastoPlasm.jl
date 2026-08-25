# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain/stress tensor types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractTensor, AbstractStrain, AbstractStress
export InfinitesimalStrain, LogarithmicStrain, CauchyStress, KirchhoffStress
export get_tensor, get_voigt, get_J2, get_τII

"""
    AbstractTensor{T}

Supertype of every per-particle second-order tensor stored on `PointSolidPhase`
(`ϵᵢⱼ`/`ϵn`/`σᵢ`/`σn`/`τᵢ`). Every concrete subtype stores the tensor as a
**volumetric + deviatoric additive split** — `(vol, dev)` for strains, `(p, dev)` for
stresses — rather than as a bare `SVector`/`SMatrix`, so the physical meaning of a
stored field is carried by its type instead of by the name of the array it lives in.

Ported from the abandoned `fancy-implementation` branch's
`src/boot/needs/types/concrete/tensor.jl`, re-implemented against the current
`Point`/`Basis` architecture (that branch's embedded `basis::B`, type-level dimension
and retained connectivity are deliberately *not* carried over — see CLAUDE.md).

# A note on the split: representational, not canonical

The pair `(vol,dev)` / `(p,dev)` is only required to satisfy the reconstruction identity

    get_voigt(x)[i] == x.dev[i,i] ± x.p    (diagonal entries)

It is **not** required that `dev` be exactly trace-free in `S` dimensions. This is
deliberate and mirrors a convention split that already existed in this codebase before
the port:

- `_trial_elastic_stress` (finite strain) uses the *plane-strain 3D* convention
  `p = -3·Kc·tr(ϵ)/3`, i.e. the volumetric part of a genuinely three-dimensional
  material response, even when `S == 2` (the out-of-plane `σzz` is implied, never
  stored). Its `dev` therefore has a non-zero 2×2 trace.
- the Voigt constructors below use the *stored-components* convention `p = -tr(σ)/S`,
  matching `σTr` (`retmap/DP.jl`) and `yield_J2` (`retmap/J2.jl`), which divide by 2 in
  2D and by 3 in 3D.

Consequently `get_J2`/`get_τII` — which genuinely need a trace-free deviator — are
defined to re-derive it from `get_voigt` rather than to trust the stored `dev`, so
they give the same number whichever writer produced the object. Use `get_voigt` /
`get_tensor` as the source of truth when consuming a stored tensor.

`get_voigt` (named for *the* Voigt-notation representation, not a generic "vector"
view) is **one function, dispatched by multiple dispatch onto two conceptually
different conventions** — see the `AbstractStrain`/`AbstractStress` sections below for
why strain and stress do not share a shear convention.

The three abstract types themselves (`AbstractTensor`/`AbstractStrain`/`AbstractStress`)
are declared in `types/abstract.jl`, alongside `AbstractConstitutiveModel`, because
`lagrangian.jl` needs them to constrain `PointSolidPhase`'s type parameters and loads
before this file.
"""
AbstractTensor

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    AbstractStrain{S,T,L} <: AbstractTensor{T}

Second-order strain tensor stored as `vol::T` (volumetric part) + `dev::SMatrix{S,S,T,L}`
(deviatoric part). `S` is the spatial dimension, `L == S*S`. Declared in
`types/abstract.jl`; concrete subtypes and methods live here.
"""
AbstractStrain

Base.getindex(strain::AbstractStrain{S,T,L}, i::Int, j::Int) where {S,T,L} =
    strain.dev[i,j] + (i == j ? strain.vol : zero(T))

"""
    get_tensor(strain::AbstractStrain{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full strain tensor `dev + vol·I`.
"""
@inline function get_tensor(strain::AbstractStrain{S,T,L}) where {S,T,L}
    return strain.dev + strain.vol * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(strain::AbstractStrain) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt-notation strain vector `(ϵxx,ϵyy,γxy)` / `(ϵxx,ϵyy,ϵzz,γyz,γxz,γxy)`, in the
**classical engineering-shear** convention (`γxy = 2ϵxy`) — *not* the tensor-shear
convention `get_tensor` uses for its internal `dev` storage. This is the standard
Voigt convention for strain (as opposed to stress, which has no such doubling): it's
why the elastic stiffness matrix `Del` is built the way it is (`σ = Del·ε_voigt` uses
the ordinary isotropic stiffness matrix only when `ε_voigt` carries this factor of 2),
and it's what makes the Voigt dot product `σ_voigt·ε_voigt` equal the tensor double
contraction `σ:ε`. Round-trips exactly with `LogarithmicStrain(ϵ::SVector)`/
`InfinitesimalStrain(ϵ::SVector)` (below): `get_voigt(LogarithmicStrain(ϵ)) == ϵ`.
"""
@inline function get_voigt(strain::AbstractStrain{2,T,L}) where {T,L}
    return SVector{3,T}(
        strain.dev[1,1] + strain.vol,
        strain.dev[2,2] + strain.vol,
        T(2.0) * strain.dev[1,2],
    )
end
@inline function get_voigt(strain::AbstractStrain{3,T,L}) where {T,L}
    return SVector{6,T}(
        strain.dev[1,1] + strain.vol,
        strain.dev[2,2] + strain.vol,
        strain.dev[3,3] + strain.vol,
        T(2.0) * strain.dev[2,3],
        T(2.0) * strain.dev[1,3],
        T(2.0) * strain.dev[1,2],
    )
end

"""
    InfinitesimalStrain{S,T,L} <: AbstractStrain{S,T,L}

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, the strain measure used under
`strain.deform == "infinitesimal"`. Built by `_infinitesimal_strain`.

Unlike `LogarithmicStrain` this type has no counterpart in `fancy-implementation`
beyond the bare struct definition (that branch left it as unused scaffolding with no
methods at all); its constructors and `_infinitesimal_strain` are new code written
here, modelled on the same `vol`/`dev` pattern and on `elast.jl`'s pre-existing
infinitesimal-strain kernel math.
"""
struct InfinitesimalStrain{S,T,L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
    InfinitesimalStrain{S,T,L}(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
    InfinitesimalStrain(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
end
@adapt_struct InfinitesimalStrain

"""
    LogarithmicStrain{S,T,L} <: AbstractStrain{S,T,L}

Logarithmic (Hencky) elastic strain tensor, the strain measure used under
`strain.deform == "finite"`. Built by `_trial_elastic_strain`.
The `SVector` constructors take the *deviatoric* part in Voigt notation.
"""
struct LogarithmicStrain{S,T,L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
    LogarithmicStrain{S,T,L}(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
    LogarithmicStrain(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
    LogarithmicStrain(vol::T, dev::SVector{3,T}) where {T} =
        new{2,T,4}(vol, SMatrix{2,2,T,4}(dev[1], dev[3], dev[3], dev[2]))
    LogarithmicStrain(vol::T, dev::SVector{6,T}) where {T} =
        new{3,T,9}(vol, SMatrix{3,3,T,9}(dev[1], dev[4], dev[5], dev[4], dev[2], dev[6], dev[5], dev[6], dev[3]))
end
@adapt_struct LogarithmicStrain

"""
    LogarithmicStrain(ϵᵢⱼ::SMatrix{S,S,T,L})

Split a full logarithmic strain tensor into `vol = tr(ϵ)/3` + `dev = ϵ - vol·I` (same
plane-strain 3D convention as `_trial_elastic_strain`).
"""
@inline function LogarithmicStrain(ϵᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    vol = tr(ϵᵢⱼ) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)))
end

"""
    LogarithmicStrain(ϵ::SVector{3,T}) / LogarithmicStrain(ϵ::SVector{6,T})

Build a `LogarithmicStrain` directly from a full **engineering-Voigt** strain vector
(`γxy = 2ϵxy`, the `get_voigt` convention) — the strain-side analogue of
`CauchyStress(σ::SVector)`/`KirchhoffStress(τ::SVector)`, letting Voigt-space return-
mapping algebra (e.g. a compliance-solved `Del\\τ`) wrap into a typed strain in one
step, with no intermediate matrix reshape. `vol = (ϵ[1]+ϵ[2])/3` (2D) /
`(ϵ[1]+ϵ[2]+ϵ[3])/3` (3D) — the same plane-strain-3D convention as the `SMatrix`
constructor above (**not** stress's `tr/S` convention — these are deliberately
different, see the `AbstractTensor` docstring). Exact inverse of `get_voigt`:
`get_voigt(LogarithmicStrain(ϵ)) == ϵ`.
"""
@inline function LogarithmicStrain(ϵ::SVector{3,T}) where {T}
    vol = (ϵ[1] + ϵ[2]) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{2,2,T,4}(ϵ[1] - vol, ϵ[3]/T(2.0), ϵ[3]/T(2.0), ϵ[2] - vol))
end
@inline function LogarithmicStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3]) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol, ϵ[6]/T(2.0), ϵ[5]/T(2.0),
        ϵ[6]/T(2.0), ϵ[2] - vol, ϵ[4]/T(2.0),
        ϵ[5]/T(2.0), ϵ[4]/T(2.0), ϵ[3] - vol,
    ))
end

# zero-initialization used by setup_mpts
Base.zero(::Type{InfinitesimalStrain{S,T,L}}) where {S,T,L} = InfinitesimalStrain(zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{LogarithmicStrain{S,T,L}})   where {S,T,L} = LogarithmicStrain(  zero(T), zero(SMatrix{S,S,T,L}))

"""
    LinearAlgebra.eigen(strain::AbstractStrain)

Eigen-decomposition of the reassembled (symmetric) strain tensor. Lets kernels write
`eigen(mpts.s.ϵn[p])` directly on a stored strain object.
"""
@inline LinearAlgebra.eigen(strain::AbstractStrain{S,T,L}) where {S,T,L} = eigen(Symmetric(get_tensor(strain)))

"""
    _trial_elastic_strain(ΔFᵢⱼ, strain::LogarithmicStrain) -> LogarithmicStrain

Push the stored logarithmic strain forward through the incremental deformation
gradient `ΔFᵢⱼ` and return the trial elastic logarithmic strain.

Numerically identical, operation for operation, to the `_logarithmic_strain` free
function it replaces in `update/fun/elast.jl` (the only change is that the volumetric/
deviatoric split is done here rather than downstream in `_kirchoff_stress`).
"""
@inline function _trial_elastic_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}, strain::LogarithmicStrain{S,T,L}) where {S,T,L}
    # trial left Cauchy-Green tensor
    λ, n = eigen(strain)
    bᵢⱼ  = ΔFᵢⱼ * (n * diagm(exp.(T(2.0) * λ)) * n') * ΔFᵢⱼ'
    # logarithmic strain tensor from the eigen-decomposition of bᵢⱼ
    λ, n = eigen(Symmetric(bᵢⱼ))
    ϵᵢⱼ  = T(0.5) * (n * diagm(log.(λ)) * n')
    # volumetric/deviatoric split (plane-strain 3D convention, see AbstractTensor docs)
    vol  = tr(ϵᵢⱼ) / T(3.0)
    dev  = ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(dev))
end

"""
    _infinitesimal_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}) -> InfinitesimalStrain

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, split into `vol = tr(ϵ)/3` and
`dev = ϵ - vol·I` (same plane-strain 3D convention as `_trial_elastic_strain`, so the
two strain types decompose consistently). `get_tensor`/`get_voigt` reproduce `ϵ`
exactly.
"""
@inline function InfinitesimalStrain(ϵᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    vol = tr(ϵᵢⱼ) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)))
end
@inline _infinitesimal_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L} =
    InfinitesimalStrain(T(0.5) .* (ΔFᵢⱼ + ΔFᵢⱼ') .- SMatrix{S,S,T,L}(I))

"""
    InfinitesimalStrain(ϵ::SVector{3,T}) / InfinitesimalStrain(ϵ::SVector{6,T})

Build an `InfinitesimalStrain` directly from a full **engineering-Voigt** strain vector
(`γxy = 2ϵxy`, the `get_voigt` convention) — same role and convention as
`LogarithmicStrain(ϵ::SVector)` above, for the infinitesimal-strain path. Exact inverse
of `get_voigt`: `get_voigt(InfinitesimalStrain(ϵ)) == ϵ`.
"""
@inline function InfinitesimalStrain(ϵ::SVector{3,T}) where {T}
    vol = (ϵ[1] + ϵ[2]) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{2,2,T,4}(ϵ[1] - vol, ϵ[3]/T(2.0), ϵ[3]/T(2.0), ϵ[2] - vol))
end
@inline function InfinitesimalStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3]) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol, ϵ[6]/T(2.0), ϵ[5]/T(2.0),
        ϵ[6]/T(2.0), ϵ[2] - vol, ϵ[4]/T(2.0),
        ϵ[5]/T(2.0), ϵ[4]/T(2.0), ϵ[3] - vol,
    ))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stress tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    AbstractStress{S,T,L} <: AbstractTensor{T}

Second-order stress tensor stored as `p::T` (pressure, **positive in compression**) +
`dev::SMatrix{S,S,T,L}` (deviatoric part), so the full tensor is `dev - p·I`.
Declared in `types/abstract.jl`; concrete subtypes and methods live here.
"""
AbstractStress

"""
    get_tensor(stress::AbstractStress{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full stress tensor `dev - p·I` (`p` positive in compression).
"""
@inline function get_tensor(stress::AbstractStress{S,T,L}) where {S,T,L}
    return stress.dev - stress.p * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(stress::AbstractStress) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt-notation stress vector `(σxx,σyy,σxy)` / `(σxx,σyy,σzz,σyz,σxz,σxy)` — the exact
component ordering every P2G/return-mapping kernel in this repo already indexes.
Unlike the `AbstractStrain` overload above, stress has only one Voigt convention (no
engineering-vs-tensor shear distinction), so this is a plain reshape with no factor of
2 anywhere.
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
    CauchyStress{S,T,L} <: AbstractStress{S,T,L}

Cauchy (true) stress. This is what `mpts.s.σᵢ`/`mpts.s.σn` store, in both the finite
(pushed forward from `τᵢ` by `transform`) and infinitesimal (Jaumann-rate updated)
paths.

The `SVector` constructors take a full Voigt **stress** vector (not a deviator, unlike
`KirchhoffStress(p,dev::SVector)`) and split it with `p = -tr(σ)/S`, matching `σTr`
(`retmap/DP.jl`). They exist so incremental Voigt-space arithmetic — the Jaumann-rate
update and `elast_fast`'s hand-inlined component math — can stay in `SVector` form
internally and be wrapped once at the point of writing into `mpts.s.σᵢ[p]`.
"""
struct CauchyStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    CauchyStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(p, dev)
    CauchyStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(p, dev)
    CauchyStress(p::T, dev::SVector{3,T}) where {T} =
        new{2,T,4}(p, SMatrix{2,2,T,4}(dev[1], dev[3], dev[3], dev[2]))
    CauchyStress(p::T, dev::SVector{6,T}) where {T} =
        new{3,T,9}(p, SMatrix{3,3,T,9}(dev[1], dev[4], dev[5], dev[4], dev[2], dev[6], dev[5], dev[6], dev[3]))
    function CauchyStress(σ::SVector{3,T}) where {T}
        p = -(σ[1] + σ[2]) / T(2.0)
        new{2,T,4}(p, SMatrix{2,2,T,4}(σ[1] + p, σ[3], σ[3], σ[2] + p))
    end
    function CauchyStress(σ::SVector{6,T}) where {T}
        p = -(σ[1] + σ[2] + σ[3]) / T(3.0)
        new{3,T,9}(p, SMatrix{3,3,T,9}(σ[1] + p, σ[6], σ[5], σ[6], σ[2] + p, σ[4], σ[5], σ[4], σ[3] + p))
    end
end
@adapt_struct CauchyStress

"""
    KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}

Kirchhoff stress `τ = J·σ`. This is what `mpts.s.τᵢ` stores under the finite-strain
path; `transform` divides it by `J` to obtain `mpts.s.σᵢ`.

Two constructor families: `(p, dev)` — used by `_trial_elastic_stress` and by
`retmap/DP.jl`, which naturally produces a returned pressure and deviator separately —
and a full Voigt vector, split like `CauchyStress`'s.
"""
struct KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    KirchhoffStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(p, dev)
    KirchhoffStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(p, dev)
    KirchhoffStress(p::T, dev::SVector{3,T}) where {T} =
        new{2,T,4}(p, SMatrix{2,2,T,4}(dev[1], dev[3], dev[3], dev[2]))
    KirchhoffStress(p::T, dev::SVector{6,T}) where {T} =
        new{3,T,9}(p, SMatrix{3,3,T,9}(dev[1], dev[4], dev[5], dev[4], dev[2], dev[6], dev[5], dev[6], dev[3]))
    function KirchhoffStress(τ::SVector{3,T}) where {T}
        p = -(τ[1] + τ[2]) / T(2.0)
        new{2,T,4}(p, SMatrix{2,2,T,4}(τ[1] + p, τ[3], τ[3], τ[2] + p))
    end
    function KirchhoffStress(τ::SVector{6,T}) where {T}
        p = -(τ[1] + τ[2] + τ[3]) / T(3.0)
        new{3,T,9}(p, SMatrix{3,3,T,9}(τ[1] + p, τ[6], τ[5], τ[6], τ[2] + p, τ[4], τ[5], τ[4], τ[3] + p))
    end
end
@adapt_struct KirchhoffStress

# zero-initialization used by setup_mpts
"""
    _trial_elastic_stress(strain::LogarithmicStrain, Kc, Gc) -> KirchhoffStress

Isotropic linear-elastic trial Kirchhoff stress from a logarithmic strain:
`p = -3·Kc·ϵvol` (positive in compression), `dev = 2·Gc·ϵdev`.

`get_voigt` of the result reproduces, expression for expression, the Voigt vector the
former `_kirchoff_stress` free function returned (`τii = 2·Gc·ϵdev_ii - P`), so the
finite-strain elastic predictor is bit-identical across the port.

Note that in 2D the returned `dev` is *not* 2×2-trace-free: the volumetric split comes
from `tr(ϵ)/3`, i.e. genuine three-dimensional plane-strain physics with an implied
out-of-plane component. See the `AbstractTensor` docstring.
"""
@inline function _trial_elastic_stress(strain::LogarithmicStrain{S,T,L}, Kc::T, Gc::T) where {S,T,L}
    P   = -T(3.0) * Kc * strain.vol
    dev =  T(2.0) * Gc * strain.dev
    return KirchhoffStress(P, dev)
end

Base.zero(::Type{CauchyStress{S,T,L}})    where {S,T,L} = CauchyStress(   zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{KirchhoffStress{S,T,L}}) where {S,T,L} = KirchhoffStress(zero(T), zero(SMatrix{S,S,T,L}))

# cross-type conversions (τ ↔ σ), used by `transform` and by the u-P dynamic_relaxation path
CauchyStress(τ::KirchhoffStress{S,T,L}, J⁻¹::T) where {S,T,L} = CauchyStress(τ.p * J⁻¹, τ.dev * J⁻¹)
KirchhoffStress(σ::CauchyStress{S,T,L}) where {S,T,L} = KirchhoffStress(σ.p, σ.dev)

"""
    get_J2(stress::AbstractStress) -> T

Second deviatoric stress invariant `J₂ = ½ sᵢⱼsᵢⱼ`, with `s` the *trace-free* deviator
re-derived from `get_voigt` (see the `AbstractTensor` docstring on why the stored
`dev` is not trusted here). Reproduces exactly the `τII²` computed by `σTr`
(`retmap/DP.jl`) and the `J2` computed by `yield_J2` (`retmap/J2.jl`), including their
`tr/2` (2D) vs `tr/3` (3D) pressure convention.
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
@inline get_τII(stress::AbstractStress) = sqrt(get_J2(stress))
