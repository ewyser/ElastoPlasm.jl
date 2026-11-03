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
    PointSolidPhase{T1,T2,D,CM,ST,EL,L,…}

The solid phase of the material points: their mechanical state, one entry per point.

# Fields
- `u::VV`: the displacement, an `SVector{D,T2}` per point: the step's, `dt·v`, on the explicit path; the total over
  the workflow under dynamic relaxation.
- `v::VV`: the velocity [m/s].
- `ρ₀::V`: the grain density [kg/m³]; `ρ::V`: the density, `(1 − n)ρ₀`.
- `m::V`: the solid mass, `(1 − n₀)ρ₀Ω₀` [kg], constant: `deform!` keeps `ρ·Ω` equal to it, so kernels read `m`.
- `Δλ::V`: the plastic multiplier of the step's return map.
- `ϵpII::V2`: the equivalent plastic strain, an `SVector{2,T2}` per point: `[1]` the point's own, `[2]` what its yield
  stress reads, the non-local average with `nonloc.status`, else a copy of `[1]`.
- `ϵpV::V`: the volumetric plastic strain; no return map updates it yet, so it stays 0.
- `σᵢ::VC`: the Cauchy stress, a `CauchyStress` (read with `get_tensor`, or `.p` and `.dev`); `σn::VC`: not read.
- `τᵢ::VK`: the Kirchhoff stress, a `KirchhoffStress`: what the finite-strain laws compute.
- `∇vᵢⱼ::VM`: the velocity gradient, an `SMatrix{D,D,T2,L}`; `∇uᵢⱼ::VM`: not read.
- `ΔFᵢⱼ::VM`, `Fᵢⱼ::VM`: the step's and the total deformation gradient; `Fn::VM`: `Fᵢⱼ` at the start of a load step
  (dynamic relaxation).
- `ϵᵢ::VS`: the elastic strain, an `ST`; `ϵn::VS`: at the start of a load step (dynamic relaxation).
- `ωᵢⱼ::VM`: the spin, the hypoelastic law's Jaumann rotation.
- `cmp::VCM`: the constitutive constants, a `CM` per point.

`T1` and `T2` are the index and float types, `D` the dimension. The parameters code dispatches on:
- `CM<:AbstractConstitutiveModel`: the constitutive model, `DruckerPrager` or `VonMises` (`material.plastic`);
- `ST<:AbstractStrain`: the strain storage, `LogarithmicStrain` for the Hencky laws, `InfinitesimalStrain` for the
  hypoelastic one (`material.elastic`);
- `EL<:AbstractElasticLaw`: the elastic law (`material.elastic`); `get_constitutive_laws()` maps the option values
  to these types;
- `L = D*D`: the length of the `D×D` fields: Julia cannot compute `D*D` in a field type, so it is a parameter, which
  kernels never dispatch on (they read a field's type with `eltype`).
The rest are the arrays' types, host or device arrays alike: `VV` of `SVector{D,T2}`, `V` of `T2`, `V2` of
`SVector{2,T2}`, `VC`/`VK` of Cauchy/Kirchhoff stresses, `VM` of `SMatrix{D,D,T2,L}`, `VS` of `ST`, `VCM` of `CM`.

Stress and strain are full 3D tensors in every dimension (plane strain in 2D), stored as a pressure or volumetric
strain and a trace-free Voigt deviator: see `AbstractStrain`, `AbstractStress`.

# Constructors
    PointSolidPhase{T1,T2,D,CM,ST,EL,L}(u, v, ρ₀, ρ, m, Δλ, ϵpII, ϵpV, σᵢ, σn, τᵢ, ∇vᵢⱼ, ∇uᵢⱼ, ΔFᵢⱼ, Fᵢⱼ, Fn, ϵᵢ, ϵn, ωᵢⱼ, cmp)

The arrays' types come from the arguments. `build_solid_phase` builds one from a case's points.

See also: `Point`, `laws`.
"""
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

"""
    PointFluidPhase{T1,T2,D}

The placeholder of a fluid phase, with no fields yet: a point's `f` is `nothing` (see "Staged path to a
poro-hydro-mechanical solution" in `planned-improvements.md`).
"""
struct PointFluidPhase{T1,T2,D}
    # Add concrete fields as needed, e.g.:
    # v    ::Matrix{T2}
end

"""
    PointThermalPhase{T1,T2,D,V,M}

The thermal phase of the material points, one entry per point (a point's `t`, when the case is thermal).

# Fields
- `c::V`: the specific heat capacity [J/(kg·K)].
- `k::V`: the thermal conductivity [W/(m·K)].
- `q::M`: the heat flux, `D × nmp`.
- `T::V`: the temperature [K].

`T1`, `T2` and `D` are the index type, float type and dimension; `V` and `M` the arrays' types (a vector and a
matrix of `T2`, host or device).

# Constructors
    PointThermalPhase{T1,T2,D}(c, k, q, T)

`build_thermal_phase` builds one from a case's points.
"""
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

"""
    Point{T1,T2,D,CM,ST,EL,L,…}

The material points (Lagrangian): where they are, their domains, volumes and porosity, and their phases, one entry
per point, `p ∈ 1:nmp`. Positions are a vector of `SVector`s: index `x[p]`, never `x[i, :]`.

# Fields
- `ndim::T1`: the dimension `D` (not read).
- `nmp::T1`: the number of points.
- `vmax::V`: the largest point speed per axis, `D` entries, for the next time step (`get_dt`).
- `x::VV`: the positions, an `SVector{D,T2}` per point [m].
- `ℓ₀::VV`, `ℓ::VV`: the initial and current half-lengths of each point's domain per axis (GIMP's, `particle_domain`).
- `n₀::V`, `n::V`: the initial and current porosity.
- `Ω₀::V`, `Ω::V`: the initial and current volume (in 2D, per unit width).
- `ΔJ::V`, `J::V`: the volume ratio of the step, `det ΔF`, and since the start, `det F`.
- `s::S`: the solid phase, a `PointSolidPhase`.
- `f::F`: the fluid phase: `nothing`, none yet (`PointFluidPhase`).
- `t::TH`: the thermal phase, a `PointThermalPhase` for a thermal case, else `nothing`.

`T1`, `T2` and `D` are the index type, float type and dimension. `CM`, `ST`, `EL` and `L` are the solid phase's
(constitutive model, strain storage, elastic law, `D*D`), which kernels dispatch on; `laws(mpts)` names the laws.
`V` and `VV` are the arrays' types (of `T2`, of `SVector{D,T2}`), host or device arrays alike (`adapt`), and `S`, `F`,
`TH` the phases'.

# Constructors
    Point(ndim, nmp, vmax, x, ℓ₀, ℓ, n₀, n, Ω₀, Ω, ΔJ, J, s, f, t)

The type parameters come from the arguments, `CM`, `ST`, `EL` and `L` from the solid phase. `setup_mpts` builds the
points of a case (`case_points`).

See also: `PointSolidPhase`, `Problem`, `laws`.
"""
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

"""
    get_constitutive_laws() -> NamedTuple

The constitutive laws a case can name, shaped like the option section `material`, each law by its option value with
the types the points carry for it:
- `elastic`: `hencky`, `improved_hencky`, `hypoelastic`, each `(; strain, law)`: its strain storage `ST` and elastic
  law `EL`;
- `plastic`: `DP`, `VM`, each its constitutive model `CM`.

The one description of the laws: `build_solid_phase` and `setup_cmp` read it by name (`table_value`), `laws(mpts)` by
type, and `get_option` lists its names, in this order.
"""
function get_constitutive_laws()
    return (; elastic = (; hencky          = (; strain=LogarithmicStrain, law=HenckySolid),
                           improved_hencky = (; strain=LogarithmicStrain, law=ImprovedHenckySolid),
                           hypoelastic     = (; strain=InfinitesimalStrain, law=HypoelasticSolid)),
              plastic = (; DP = DruckerPrager,
                           VM = VonMises))
end

"""
    laws(mpts::Point) -> (; elastic, plastic)

The names of the points' constitutive laws, as a case's options give them (`material.elastic`, `material.plastic`),
read from their types (`EL`, `CM`) in `get_constitutive_laws()`.
"""
function laws(::Point{T1,T2,D,CM,ST,EL}) where {T1,T2,D,CM,ST,EL}
    t = get_constitutive_laws()
    elastic = string(findfirst(l -> l.law == EL, t.elastic))
    plastic = string(findfirst(M -> CM <: M, t.plastic))
    return (; elastic, plastic)
end
