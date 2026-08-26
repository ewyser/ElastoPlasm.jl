
"""
    build_solid_phase(T1,T2,D,solver,mat,geom,nmp,xp,vp,ρ0,TM,TV,TS) -> (s, cmp, CM, ST, SC, SK)

Build the per-particle solid-phase state (`PointSolidPhase`) — per-particle
constitutive-model bundle (`cmp`, via `setup_cmp`), and all mechanical state fields,
zero-initialized except `v`/`ρ`. `CM`/`ST`/`SC`/`SK` are also returned since the caller
needs them to build `Point`'s type parameters.

`ST` (the typed strain storage of `ϵᵢⱼ`/`ϵn`, see `concrete/tensor.jl`) is picked by
`solver.strain.deform`: `LogarithmicStrain` for `"finite"`, `InfinitesimalStrain` for
`"infinitesimal"`. `SC`/`SK` (`σᵢⱼ`/`σn` and `τᵢⱼ`) are always `CauchyStress`/
`KirchhoffStress` regardless of `deform`, exactly as both field families existed
unconditionally before the port.
"""
function build_solid_phase(T1,T2,D,solver,mat,geom,nmp,xp,vp,ρ0,TM,TV,TS)
    L = D*D
    if solver.strain.deform == "finite"
        ST    = LogarithmicStrain{D,T2,L}
    elseif solver.strain.deform == "infinitesimal"
        ST    = InfinitesimalStrain{D,T2,L}
    end
    SC = CauchyStress{D,T2,L}
    SK = KirchhoffStress{D,T2,L}

    cmp = setup_cmp(nmp,T2.(vec(copy(geom.coh0))),T2.(vec(copy(geom.cohr))),T2.(vec(copy(geom.phi))); E=T2(mat[:E]),ν=T2(mat[:ν]),Hp=T2(mat[:Hp]),D=Int(D))
    CM  = eltype(cmp)

    s = PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}(
        [zero(TV)  for _ in 1:nmp]                         , # u
        [TV(T2.(vp[:,p])) for p in 1:nmp]                  , # v
        # mechanical properties
        T2.(vec(copy(ρ0)))                                  , # ρ₀
        T2.(vec(copy(ρ0)))                                  , # ρ
        T2.(zeros(nmp))                                     , # Δλ
        [zero(SVector{2,T2}) for _ in 1:nmp]                , # ϵpII
        T2.(zeros(nmp))                                     , # ϵpV
        # typed stress tensors (concrete/tensor.jl)
        [zero(SC) for _ in 1:nmp]                          , # σᵢⱼ
        [zero(SC) for _ in 1:nmp]                          , # σn
        [zero(SK) for _ in 1:nmp]                          , # τᵢⱼ
        # tensor in matrix notation
        [zero(TM) for _ in 1:nmp]                          , # ∇vᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ∇uᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ΔFᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fn
        [zero(ST) for _ in 1:nmp]                          , # ϵᵢⱼ
        [zero(ST) for _ in 1:nmp]                          , # ϵn
        [zero(TM) for _ in 1:nmp]                          , # ωᵢⱼ
        # per-particle constitutive-model constants
        cmp                                                  , # cmp::Vector{CM}
    )
    return s, cmp, CM, ST, SC, SK
end

"""
    build_thermal_phase(T1,T2,D,geom,nmp; thermal::Bool=false) -> Union{Nothing,PointThermalPhase}

Build the per-particle thermal-phase state (`PointThermalPhase`) from
`geom.c`/`geom.k`/`geom.T` when `thermal=true` (see `get_thermal`), otherwise `nothing`.
"""
function build_thermal_phase(T1,T2,D,geom,nmp; thermal::Bool=false)
    return thermal ? PointThermalPhase{T1,T2,D}(
        T2.(vec(copy(geom.c)))                            , # c::Vector{T2} specific heat capacity vector
        T2.(vec(copy(geom.k)))                             , # k::Vector{T2} thermal conductivity vector
        T2.(zeros(D,nmp))                                  , # q::Matrix{T2} heat flux array
        T2.(vec(copy(geom.T)))                            , # T::Vector{T2} temperature vector
    ) : nothing
end

"""
    setup_mpts(mesh::Mesh{T1,T2,D}, solver::S, mat::NamedTuple; geom::NamedTuple=(), thermal::Bool=false) -> Point

Construct the material point system (mpts) for a simulation, initializing all relevant fields.
Connectivity (`p2n`, `p2e`, `e2p`, `p2p`) is built separately by `setup_basis`, after both `Mesh`
and `Point` exist.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object containing geometry.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`) with matching `T1,T2,D`.
- `mat::NamedTuple`: Raw material constants (`setup_material_constants`) — transient, used only to
  build the per-particle `cmp::Vector{<:AbstractConstitutiveModel}` here (via `setup_cmp`),
  not persisted or stored on `Point` itself.
- `geom::NamedTuple=()`: (Optional) Geometry definition (e.g., number of intervals, number of material points, geometry struct).
  When `thermal=true`, must additionally carry `c`/`k`/`T` (per-particle specific heat capacity,
  thermal conductivity, initial temperature — see `get_thermal`).
- `thermal::Bool=false`: (Optional) Build a real `PointThermalPhase` from `geom.c`/`geom.k`/`geom.T`
  instead of leaving `mpts.t` as `nothing`. Only `thermal_problem` sets this.

# Returns
- `Point`: Material point data structure with all fields initialized for simulation.

# Example
```julia
mpts = setup_mpts(mesh, solver, mat; geom=get_slump(mesh, mat, solver))
println(mpts.nmp)  # Number of material points
```

# Notes
- Initializes material point positions, volumes, densities, and all state variables. Phase
  construction itself is delegated to `build_solid_phase`/`build_thermal_phase` above.
- Handles both 2D and 3D cases.
"""
function setup_mpts(mesh::Mesh{T1,T2,D},solver::S,mat::NamedTuple; geom::NamedTuple=(;), thermal::Bool=false) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    props = mesh.prprt
    # non-dimensional constant
    if D == 2
        nstr = 3
    elseif D == 3
        nstr = 6
    else
        error("Unsupported dimension: $D")
    end
    # unpack material geometry
    ni,nmp,xp = geom.ni,geom.nmp,geom.xp
    # scalars & vectors
    n0 = 0.1.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(props.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = fill(mat[:ρ0],nmp)
    # initial velocity (if provided)
    vp = haskey(geom, :vp) ? geom.vp : zeros(size(xp))

    # static array types for the new AoS memory layout
    if D == 2
        TM = SMatrix{2,2,T2,4}
        TV = SVector{2,T2}
        TS = SVector{3,T2}
    else
        TM = SMatrix{3,3,T2,9}
        TV = SVector{3,T2}
        TS = SVector{6,T2}
    end

    # constructor - create components
    s, cmp, CM, ST, SC, SK = build_solid_phase(T1,T2,D,solver,mat,geom,nmp,xp,vp,ρ0,TM,TV,TS)
    t = build_thermal_phase(T1,T2,D,geom,nmp; thermal=thermal)

    mpts = Point{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}(
        # general information
        T1(D)                              , # ndim
        T1(nmp)                              , # nmp
        T2.(zeros(D))                      , # vmax
        # basis-related quantities
        T2.(zeros(props.nn,D,nmp        ))  , # Δnp
        # APIC-related
        [zero(TM) for _ in 1:nmp]          , # Bᵢⱼ
        [zero(TM) for _ in 1:nmp]          , # Dᵢⱼ
        # connectivity
        T1(props.nn)                          , # nn
        # material point properties
        [TV(T2.(xp[:,i])) for i in 1:nmp]           , # x
        [TV(T2.(l0[:,i])) for i in 1:nmp]           , # ℓ₀
        [TV(T2.(l0[:,i])) for i in 1:nmp]           , # ℓ
        T2.(copy(n0))                        , # n₀
        T2.(copy(n0))                        , # n
        T2.(vec(copy(v0)))                   , # Ω₀
        T2.(vec(copy(v0)))                   , # Ω
        T2.(ones(nmp))                       , # ΔJ
        T2.(ones(nmp))                       , # J
        # solid phase
        s                                    , #
        # fluid phase
        nothing                              , #
        # thermal phase
        t                                    , #
    )
    return mpts
end
