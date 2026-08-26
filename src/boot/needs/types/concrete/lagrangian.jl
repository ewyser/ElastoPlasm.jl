# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Point,PointSolidPhase,PointFluidPhase,PointThermalPhase

"""
    PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}

Per-particle solid-phase state. Beyond the historical `TM`/`TV`/`TS` static-array
shape parameters, three trailing parameters carry the *typed tensor* storage
introduced by the `tensor.jl` port (see `AbstractTensor`):

- `ST<:AbstractStrain` — `ϵᵢⱼ`/`ϵn`. One parameter for both, since the field pair is
  dual-purpose: `InfinitesimalStrain` under `strain.deform="infinitesimal"`,
  `LogarithmicStrain` under `strain.deform="finite"`. Resolved in `build_solid_phase`
  by the `solver.strain.deform` branch.
- `SC<:AbstractStress` — `σᵢⱼ`/`σn`, always `CauchyStress`.
- `SK<:AbstractStress` — `τᵢⱼ`, always `KirchhoffStress`.

They are deliberately *trailing*: nearly every kernel signature in this repo
pattern-matches only `Point{T1,T2,D}`, so adding them at the end left those
signatures untouched.
"""
struct PointSolidPhase{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SC<:AbstractStress,SK<:AbstractStress} <: AbstractMaterialPointPhase{T1,T2}
    u    ::Vector{TV}   # displacement per MP : SVector{ndim,T2}
    v    ::Vector{TV}   # velocity per MP     : SVector{ndim,T2}
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    Δλ   ::Vector{T2}
    ϵpII ::Vector{SVector{2,T2}}
    ϵpV  ::Vector{T2}
    # typed stress tensors (see concrete/tensor.jl); Voigt view via `get_voigt`
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

struct Point{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SC<:AbstractStress,SK<:AbstractStress} <: AbstractMaterialPoint{T1,T2}
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
    s    ::PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}
    # fluid phase
    f    ::Union{Nothing, PointFluidPhase{T1,T2,D}}
    # thermal phase
    t    ::Union{Nothing, PointThermalPhase{T1,T2,D}}
end
@adapt_struct Point
