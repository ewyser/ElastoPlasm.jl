# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractLagrangian end
abstract type AbstractMaterialPoint{T1, T2} <: AbstractLagrangian end
abstract type AbstractMaterialPointPhase{T1, T2} <: AbstractMaterialPoint{T1,T2} end

export Point,PointSolidPhase,PointFluidPhase,PointThermalPhase
export AbstractSolid,LinearSolid,HenckySolid,ImprovedHenckySolid

"""
    AbstractSolid

Dispatch tag for the solid-phase elasto-kinematic formulation, carried as `Point`'s
own `SM` type parameter (see `PointSolidPhase` below). A run-wide choice, resolved
entirely from `solver.material.elastic` (`"linear"`/`"hencky"`/`"improved hencky"` —
there is no separate strain-formulation config key), never per-particle, resolved
once by `build_solid_phase`. `elast.jl`'s three kernel methods dispatch on `SM`
directly instead of branching on a config string, mirroring how `retmap` already
dispatches on `Point`'s `CM` type with no branch anywhere.

Julia cannot make a struct field's type a function of another type parameter (a
field's declared type is elaborated once, generically over the `TypeVar`, not
per-instantiation — confirmed empirically before adding this), so `ST<:AbstractStrain`
(the actual `ϵᵢⱼ`/`ϵn` storage type) stays a real, independent type parameter
alongside `SM` rather than being derived from it. The two are always built together,
consistently, by the single construction site (`build_solid_phase`) — `HenckySolid`/
`ImprovedHenckySolid` always pair with `ST=LogarithmicStrain`, `LinearSolid` always
with `ST=InfinitesimalStrain`.
"""
abstract type AbstractSolid end

"""
    LinearSolid <: AbstractSolid

Infinitesimal-strain solid: Jaumann-rate Cauchy stress update (`elast.jl`'s
`ST<:InfinitesimalStrain` kernel). Always paired with `ST=InfinitesimalStrain`.
"""
struct LinearSolid <: AbstractSolid end

"""
    HenckySolid <: AbstractSolid

Finite-strain solid with the original constant-`Kc` linear-Hencky volumetric law
(`p = -Kc·tr(ϵ)`). Default. Always paired with `ST=LogarithmicStrain`.
"""
struct HenckySolid <: AbstractSolid end

"""
    ImprovedHenckySolid <: AbstractSolid

Finite-strain solid with the porosity-weighted "improved Hencky" volumetric law
(Pretti, Coombs, Augarde, Marchena Puigvert, Reyna Gutiérrez, *Mechanics of
Materials* 192 (2024) 104958, Eqs. 33/34, elastic part only — see `elast.jl`). Always
paired with `ST=LogarithmicStrain`.
"""
struct ImprovedHenckySolid <: AbstractSolid end

"""
    PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SM,SC,SK}

Per-particle solid-phase state. Beyond the historical `TM`/`TV`/`TS` static-array
shape parameters, four trailing parameters carry the *typed tensor* storage
introduced by the tensor port (see `AbstractStrain` in `strain.jl` /
`AbstractStress` in `stress.jl`) plus the solid-formulation dispatch tag:

- `ST<:AbstractStrain` — `ϵᵢⱼ`/`ϵn`. One parameter for both, since the field pair is
  dual-purpose: `InfinitesimalStrain` under `material.elastic="linear"`,
  `LogarithmicStrain` under `material.elastic∈{"hencky","improved hencky"}`. Resolved
  in `build_solid_phase` by the `solver.material.elastic` branch.
- `SM<:AbstractSolid` — `LinearSolid`/`HenckySolid`/`ImprovedHenckySolid`, see
  `AbstractSolid`'s docstring. Always consistent with `ST` by construction.
- `SC<:AbstractStress` — `σᵢⱼ`/`σn`, always `CauchyStress`.
- `SK<:AbstractStress` — `τᵢⱼ`, always `KirchhoffStress`.

They are deliberately *trailing*: nearly every kernel signature in this repo
pattern-matches only `Point{T1,T2,D}` or `Point{T1,T2,D,CM,TM,TV,TS,ST}`, so adding
`SM` right after `ST` (ahead of `SC`/`SK`) left those signatures untouched — only
`elast.jl` needs `SM` explicitly.
"""
struct PointSolidPhase{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SM<:AbstractSolid,SC<:AbstractStress,SK<:AbstractStress} <: AbstractMaterialPointPhase{T1,T2}
    u    ::Vector{TV}   # displacement per MP : SVector{ndim,T2}
    v    ::Vector{TV}   # velocity per MP     : SVector{ndim,T2}
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    Δλ   ::Vector{T2}
    ϵpII ::Vector{SVector{2,T2}}
    ϵpV  ::Vector{T2}
    # typed stress tensors (see stress.jl); Voigt view via `get_voigt`
    σᵢⱼ  ::Vector{SC}
    σn   ::Vector{SC}
    τᵢⱼ  ::Vector{SK}
    # tensor in matrix notation (SMatrix{ndim,ndim,T2} per MP)
    ∇vᵢⱼ ::Vector{TM}
    ∇uᵢⱼ ::Vector{TM}
    ΔFᵢⱼ ::Vector{TM}
    Fᵢⱼ  ::Vector{TM}
    Fn   ::Vector{TM}
    ϵᵢⱼ  ::Vector{ST}
    ϵn   ::Vector{ST}
    ωᵢⱼ  ::Vector{TM}
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

struct Point{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SM<:AbstractSolid,SC<:AbstractStress,SK<:AbstractStress} <: AbstractMaterialPoint{T1,T2}
    # general information
    ndim ::T1
    nmp  ::T1
    # CFL-related quantity
    vmax ::Vector{T2}
    # connectivity
    nn   ::T1
    # material point properties
    x    ::Vector{TV}  # coordinates per MP : SVector{ndim,T2}
    ℓ₀   ::Vector{TV}  # reference domain half-lengths per MP
    ℓ    ::Vector{TV}  # current domain half-lengths per MP
    n₀   ::Vector{T2}
    n    ::Vector{T2}    
    Ω₀   ::Vector{T2}
    Ω    ::Vector{T2}
    ΔJ   ::Vector{T2}
    J    ::Vector{T2}
    # solid phase
    s    ::PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SM,SC,SK}
    # fluid phase
    f    ::Union{Nothing, PointFluidPhase{T1,T2,D}}
    # thermal phase
    t    ::Union{Nothing, PointThermalPhase{T1,T2,D}}
end
@adapt_struct Point
