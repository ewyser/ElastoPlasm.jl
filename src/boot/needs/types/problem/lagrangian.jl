# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractLagrangian end
abstract type AbstractMaterialPoint{T1, T2} <: AbstractLagrangian end
abstract type AbstractMaterialPointPhase{T1, T2} <: AbstractMaterialPoint{T1,T2} end

export Point,PointSolidPhase,PointFluidPhase,PointThermalPhase
export AbstractElasticLaw,HypoelasticSolid,HenckySolid,ImprovedHenckySolid

"""
    AbstractElasticLaw

Elastic law of the solid phase, carried as `Point`'s `EL` type parameter and picked
from `material.elastic` by `build_solid_phase`. `elast.jl` and `get_dt`'s `_Ktan`
dispatch on it. `ST` can't be derived from `EL` in a field type, so both are separate
parameters, always built together: `HenckySolid`/`ImprovedHenckySolid` with
`LogarithmicStrain`, `HypoelasticSolid` with `InfinitesimalStrain`.
"""
abstract type AbstractElasticLaw end

"""
    HypoelasticSolid <: AbstractElasticLaw

Small-strain hypoelastic solid: Jaumann-rate Cauchy stress update (`elast.jl`'s
`ST<:InfinitesimalStrain` kernel). Always paired with `ST=InfinitesimalStrain`.
"""
struct HypoelasticSolid <: AbstractElasticLaw end

"""
    HenckySolid <: AbstractElasticLaw

Finite-strain solid with the original constant-`Kc` linear-Hencky volumetric law
(`p = -Kc·tr(ϵ)`). Default. Always paired with `ST=LogarithmicStrain`.
"""
struct HenckySolid <: AbstractElasticLaw end

"""
    ImprovedHenckySolid <: AbstractElasticLaw

Finite-strain solid with the porosity-weighted "improved Hencky" volumetric law
(Pretti, Coombs, Augarde, Marchena Puigvert, Reyna Gutiérrez, *Mechanics of
Materials* 192 (2024) 104958, Eqs. 33/34, elastic part only — see `elast.jl`). Always
paired with `ST=LogarithmicStrain`.
"""
struct ImprovedHenckySolid <: AbstractElasticLaw end

"""
    PointSolidPhase{T1,T2,D,CM,ST,EL,L}

Per-particle solid-phase state. Beyond the index/float types `T1`/`T2` and dimension `D`:

- `CM<:AbstractConstitutiveModel` — plastic model (`DruckerPrager`/`VonMises`), from
  `material.plastic`.
- `ST<:AbstractStrain` — `ϵᵢ`/`ϵn` storage, from `material.elastic`:
  `InfinitesimalStrain` for `"hypoelastic"`, `LogarithmicStrain` for the Hencky laws.
- `EL<:AbstractElasticLaw` — elastic law, from `material.elastic`; see
  `AbstractElasticLaw`.
- `L == D*D` — static length of the `D×D` kinematic fields (`SMatrix{D,D,T2,L}`). Julia can't
  compute `D*D` in a field type, so it is carried as a parameter. Kernels never dispatch on it;
  recover a field's type with `eltype` instead.

Stress and strain (`σᵢ`, `σn`, `τᵢ`, `ϵᵢ`, `ϵn`) are full 3D Voigt vectors in every dimension (plane
strain in 2D), see `AbstractStrain`/`AbstractStress`.

`m` is the solid mass `(1−n₀)ρ₀Ω₀`, set once at setup: `deform!` keeps `ρ·Ω` equal to it
(solid mass conservation), so kernels read `m` rather than recomputing `ρ·Ω`.
"""
struct PointSolidPhase{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L} <: AbstractMaterialPointPhase{T1,T2}
    u    ::Vector{SVector{D,T2}}   # displacement per MP
    v    ::Vector{SVector{D,T2}}   # velocity per MP
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    m    ::Vector{T2}   # solid mass per MP, constant (solid mass conservation)
    Δλ   ::Vector{T2}
    ϵpII ::Vector{SVector{2,T2}}
    ϵpV  ::Vector{T2}
    # typed stress tensors (see tensor.jl): pressure `.p` and deviator `.dev`; tensor via `get_tensor`
    σᵢ   ::Vector{CauchyStress{T2}}
    σn   ::Vector{CauchyStress{T2}}
    τᵢ   ::Vector{KirchhoffStress{T2}}
    # tensor in matrix notation (SMatrix{ndim,ndim,T2} per MP)
    ∇vᵢⱼ ::Vector{SMatrix{D,D,T2,L}}
    ∇uᵢⱼ ::Vector{SMatrix{D,D,T2,L}}
    ΔFᵢⱼ ::Vector{SMatrix{D,D,T2,L}}
    Fᵢⱼ  ::Vector{SMatrix{D,D,T2,L}}
    Fn   ::Vector{SMatrix{D,D,T2,L}}
    ϵᵢ   ::Vector{ST}
    ϵn   ::Vector{ST}
    ωᵢⱼ  ::Vector{SMatrix{D,D,T2,L}}
    # per-particle static constitutive-model constants (Gc,Kc,Del,Hp,c₀,cᵣ,ϕ₀, ...)
    cmp  ::Vector{CM}
end
@adapt_struct PointSolidPhase

struct PointFluidPhase{T1,T2,D} <: AbstractMaterialPointPhase{T1,T2}
    # Add concrete fields as needed, e.g.:
    # v    ::Matrix{T2}
end
@adapt_struct PointFluidPhase

struct PointThermalPhase{T1,T2,D} <: AbstractMaterialPointPhase{T1,T2}
    c   ::Vector{T2} # specific heat capacity vector
    k   ::Vector{T2} # thermal conductivity vector
    q   ::Matrix{T2} # heat flux array
    T   ::Vector{T2} # temperature vector
end
@adapt_struct PointThermalPhase

struct Point{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L} <: AbstractMaterialPoint{T1,T2}
    # general information
    ndim ::T1
    nmp  ::T1
    # CFL-related quantity
    vmax ::Vector{T2}
    # connectivity
    nn   ::T1
    # material point properties
    x    ::Vector{SVector{D,T2}}  # coordinates per MP
    ℓ₀   ::Vector{SVector{D,T2}}  # reference domain half-lengths per MP
    ℓ    ::Vector{SVector{D,T2}}  # current domain half-lengths per MP
    n₀   ::Vector{T2}
    n    ::Vector{T2}    
    Ω₀   ::Vector{T2}
    Ω    ::Vector{T2}
    ΔJ   ::Vector{T2}
    J    ::Vector{T2}
    # solid phase
    s    ::PointSolidPhase{T1,T2,D,CM,ST,EL,L}
    # fluid phase
    f    ::Union{Nothing, PointFluidPhase{T1,T2,D}}
    # thermal phase
    t    ::Union{Nothing, PointThermalPhase{T1,T2,D}}
end
@adapt_struct Point
