"""
    get_elastic_moduli(E, ν) -> Kc, Gc

Bulk and shear moduli from Young's modulus and Poisson's ratio. The Voigt stiffness `Del` is built
from them on the fly by `get_Del` (`stiffness.jl`).
"""
function get_elastic_moduli(E::T2,ν::T2) where {T2}
    return T2(E/(3.0*(1.0-2.0*ν))), T2(E/(2.0*(1.0+ν)))
end

"""
    setup_material_constants(solver::AbstractSolver{T1,T2,D}; E::T2=1.0e6, ν::T2=0.3, ρ0::T2=2700.0) -> NamedTuple

Compute the raw scalar material constants (elastic + plastic + thermal) for a
simulation, based on dimensionality. This is a **transient, setup-time-only**
NamedTuple (formerly `setup_cmpr`/`cmpr`) — it is consumed by `get_slump`/`get_thermal`/
`get_collision`/`get_column`/`get_collapse`/`mpts_populate` (to build possibly-heterogeneous
per-particle geometry fields, e.g. under GRF) and by `setup_cmp` (to build the
per-particle `cmp::Vector{<:AbstractConstitutiveModel}` that actually lands on
`Point`), but it is never itself persisted to JLD2 and never threaded into any
solver/update/get kernel — those all read `mpts.s.cmp[p]` instead.

Takes `solver` rather than `mesh` deliberately: `AbstractSolver{T1,T2,D}` already
carries `T1`/`T2`/`D`, the only things this function needs, so callers don't have to
build a (possibly throwaway) `Mesh` just to hand this function a dimension.

# Arguments
- `solver::AbstractSolver{T1,T2,D}`: Solver instance, used only for its `T1,T2,D` type parameters.
- `E::T2=1.0e6`: (Optional) Young's modulus (Pa).
- `ν::T2=0.3`: (Optional) Poisson's ratio.
- `ρ0::T2=2700.0`: (Optional) Initial density (kg/m³).

# Returns
- `NamedTuple`: Raw material constants, including elastic, plastic, and thermal properties.

# Example
```julia
mat = setup_material_constants(solver; E=2.0e6, ν=0.25, ρ0=2500.0)
println(mat.Kc)  # Bulk modulus
```
"""
function setup_material_constants(solver::AbstractSolver{T1,T2,D};
    E::T2=T2(1.0e6),
    ν::T2=T2(0.3),
    ρ0::T2= T2(2700.0),
    ) where {T1,T2,D}
    # independant physical constant
    K,G     = get_elastic_moduli(E,ν)                                           # bulk & shear modulus [Pa] from Young's mod. + Poisson's ratio
    c       = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]
    Hp      = -30.0e3                                                           # softening modulus
    # material constants
    mat = (;
        # solid phase
        E   = T2(E),
        ν   = T2(ν),
        Kc  = T2(K),
        Gc  = T2(G),
        Hp  = T2(Hp),
        c0  = T2(c0),
        cr  = T2(cr),
        ϕ0  = T2(ϕ0),
        ϕr  = T2(ϕr),
        ρ0  = T2(ρ0),
        c   = T2(c),
        # thermal phase
        specific_heat_capacity = T2(10.0),   # J/(kg·K)
        thermal_conductivity   = T2(3000.0), # W/(m·K)
        initial_temperature    = T2(293.15), # K, 20°C in Kelvin
    )
    return mat::NamedTuple
end

"""
    setup_cmp(nmp, coh0, cohr, phi; E, ν, Hp, D, constitutive) -> Vector{<:AbstractConstitutiveModel}

Build one per-particle `AbstractConstitutiveModel` instance for every material point,
from the (possibly spatially-heterogeneous, e.g. GRF-perturbed) per-particle cohesion/
friction fields already computed by `get_slump`/`get_thermal`/`get_collision`/
`get_column`/`get_collapse`/`mpts_populate`, plus the uniform elastic constants derived once from
`E`,`ν` via `get_elastic_moduli`. Branches on `constitutive` (typically
`solver.material.plastic`) to build `Vector{DruckerPrager}` (`"DP"`) or
`Vector{VonMises}` (`"VM"`) — never `PerfectlyElastic`, see its docstring for why. An
unrecognized `constitutive` string throws immediately here, at setup time, rather than
deferring to a runtime error the first time the retmap kernel is invoked.

Which elastic trial-stress law applies (`material.elastic`) is resolved separately, in
`build_solid_phase`, onto `Point`'s own `EL<:AbstractElasticLaw` type parameter — not here;
`cmp` stays scoped purely to the plastic-model constants.

# Arguments
- `nmp`: Number of material points.
- `coh0`,`cohr`,`phi`: Per-particle initial cohesion, residual cohesion, friction angle.
  `phi` is ignored when `constitutive=="VM"` (von Mises has no friction angle) — still
  required positionally so callers don't need to conditionally omit it.
- `E`,`ν`: Young's modulus, Poisson's ratio (used to derive `Kc`,`Gc`).
- `Hp`: Softening modulus (uniform across particles).
- `D`: Spatial dimension.
- `constitutive`: `"DP"` or `"VM"`, selects the concrete `AbstractConstitutiveModel`.
"""
function setup_cmp(nmp::Integer,coh0::AbstractVector{T2},cohr::AbstractVector{T2},phi::AbstractVector{T2};
    E::T2,ν::T2,Hp::T2,D::Integer,constitutive::String,
    ) where {T2}
    Kc,Gc = get_elastic_moduli(E,ν)
    if constitutive == "DP"
        return [DruckerPrager{T2,D}(Gc,Kc,Hp,T2(coh0[p]),T2(cohr[p]),T2(phi[p])) for p ∈ 1:nmp]
    elseif constitutive == "VM"
        return [VonMises{T2,D}(Gc,Kc,Hp,T2(coh0[p]),T2(cohr[p])) for p ∈ 1:nmp]
    else
        throw(error("InvalidConstitutiveModel: $(constitutive) (expected \"DP\" or \"VM\")"))
    end
end