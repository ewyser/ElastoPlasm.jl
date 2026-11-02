# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
# The arrays' types are parameters appended after the dispatched ones (`VV`: of SVector{D,T2}, `V`: of T2, `V2`: of
# SVector{2,T2}, `VC`/`VK`: of Cauchy/Kirchhoff stress, `VM`: of SMatrix{D,D,T2,L}, `VS`: of ST, `VCM`: of CM), so the
# struct holds host or device arrays alike (see eulerian.jl).
struct PointSolidPhase{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L,
                       VV,V,V2,VC,VK,VM,VS,VCM}
    u    ::VV   # displacement per MP
    v    ::VV   # velocity per MP
    # mechanical properties
    ρ₀   ::V
    ρ    ::V
    m    ::V    # solid mass per MP, constant (solid mass conservation)
    Δλ   ::V
    ϵpII ::V2
    ϵpV  ::V
    # typed stress tensors (see tensor.jl): pressure `.p` and deviator `.dev`; tensor via `get_tensor`
    σᵢ   ::VC
    σn   ::VC
    τᵢ   ::VK
    # tensor in matrix notation (SMatrix{ndim,ndim,T2} per MP)
    ∇vᵢⱼ ::VM
    ∇uᵢⱼ ::VM
    ΔFᵢⱼ ::VM
    Fᵢⱼ  ::VM
    Fn   ::VM
    ϵᵢ   ::VS
    ϵn   ::VS
    ωᵢⱼ  ::VM
    # per-particle static constitutive-model constants (Gc,Kc,Del,Hp,c₀,cᵣ,ϕ₀, ...)
    cmp  ::VCM
end
function PointSolidPhase{T1,T2,D,CM,ST,EL,L}(u::VV, v::VV, ρ₀::V, ρ::V, m::V, Δλ::V, ϵpII::V2, ϵpV::V, σᵢ::VC, σn::VC,
                                             τᵢ::VK, ∇vᵢⱼ::VM, ∇uᵢⱼ::VM, ΔFᵢⱼ::VM, Fᵢⱼ::VM, Fn::VM, ϵᵢ::VS, ϵn::VS,
                                             ωᵢⱼ::VM, cmp::VCM) where {T1,T2,D,CM,ST,EL,L,VV,V,V2,VC,VK,VM,VS,VCM}
    return PointSolidPhase{T1,T2,D,CM,ST,EL,L,VV,V,V2,VC,VK,VM,VS,VCM}(u, v, ρ₀, ρ, m, Δλ, ϵpII, ϵpV, σᵢ, σn, τᵢ,
                                                                      ∇vᵢⱼ, ∇uᵢⱼ, ΔFᵢⱼ, Fᵢⱼ, Fn, ϵᵢ, ϵn, ωᵢⱼ, cmp)
end
function Adapt.adapt_structure(to, s::PointSolidPhase{T1,T2,D,CM,ST,EL,L}) where {T1,T2,D,CM,ST,EL,L}
    return adapt_fields(to, s, PointSolidPhase{T1,T2,D,CM,ST,EL,L})   # pins T1 and EL, which no field carries
end

struct PointFluidPhase{T1,T2,D}
    # Add concrete fields as needed, e.g.:
    # v    ::Matrix{T2}
end

struct PointThermalPhase{T1,T2,D,V,M}
    c   ::V # specific heat capacity vector
    k   ::V # thermal conductivity vector
    q   ::M # heat flux array
    T   ::V # temperature vector
end
function PointThermalPhase{T1,T2,D}(c::V, k::V, q::M, T::V) where {T1,T2,D,V,M}
    return PointThermalPhase{T1,T2,D,V,M}(c, k, q, T)
end
function Adapt.adapt_structure(to, t::PointThermalPhase{T1,T2,D}) where {T1,T2,D}
    return adapt_fields(to, t, PointThermalPhase{T1,T2,D})
end

struct Point{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L,V,VV,S,F,TH}
    # general information
    ndim ::T1
    nmp  ::T1
    # CFL-related quantity
    vmax ::V
    # material point properties
    x    ::VV  # coordinates per MP
    ℓ₀   ::VV  # reference domain half-lengths per MP
    ℓ    ::VV  # current domain half-lengths per MP
    n₀   ::V
    n    ::V
    Ω₀   ::V
    Ω    ::V
    ΔJ   ::V
    J    ::V
    # solid phase
    s    ::S
    # fluid phase (nothing without one)
    f    ::F
    # thermal phase (nothing without one)
    t    ::TH
end
# the dispatched parameters are the solid phase's
function Point(ndim::T1, nmp::T1, vmax::V, x::VV, ℓ₀::VV, ℓ::VV, n₀::V, n::V, Ω₀::V, Ω::V, ΔJ::V, J::V,
               s::PointSolidPhase{T1,T2,D,CM,ST,EL,L}, f::F, t::TH) where {T1,T2,D,CM,ST,EL,L,V,VV,F,TH}
    return Point{T1,T2,D,CM,ST,EL,L,V,VV,typeof(s),F,TH}(ndim, nmp, vmax, x, ℓ₀, ℓ, n₀, n, Ω₀, Ω, ΔJ, J, s, f, t)
end
function Adapt.adapt_structure(to, mpts::Point)
    return adapt_fields(to, mpts, Point)
end

# The constitutive laws a case names (its options `material.elastic` and `material.plastic`) and the types the points
# carry for them: the strain storage `ST` with the elastic law `EL`, and the constitutive model `CM`. Read one way to
# build the points (`build_solid_phase`, `setup_cmp`) and to list the options (`get_option`), the other to name a
# problem's laws (`laws`).
const ELASTIC_LAWS = ("hencky"          => (LogarithmicStrain, HenckySolid),
                      "improved_hencky" => (LogarithmicStrain, ImprovedHenckySolid),
                      "hypoelastic"     => (InfinitesimalStrain, HypoelasticSolid))
const PLASTIC_LAWS = ("DP" => DruckerPrager,
                      "VM" => VonMises)

"""
    laws(mpts::Point) -> (; elastic, plastic)

The names of the points' constitutive laws, as a case's options give them (`material.elastic`, `material.plastic`),
read from their types (`EL`, `CM`).
"""
function laws(::Point{T1,T2,D,CM,ST,EL}) where {T1,T2,D,CM,ST,EL}
    elastic = first(ELASTIC_LAWS[findfirst(l -> last(l)[2] == EL, ELASTIC_LAWS)])
    plastic = first(PLASTIC_LAWS[findfirst(l -> CM <: last(l), PLASTIC_LAWS)])
    return (; elastic, plastic)
end
