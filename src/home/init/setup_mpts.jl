
"""
    setup_mpts(mesh::Mesh{T1,T2,D}, solver::S, mat::NamedTuple; geom::NamedTuple=()) -> Point

Construct the material point system (mpts) for a simulation, initializing all relevant fields.
Connectivity (`p2n`, `p2e`, `e2p`, `p2p`) is built separately by `setup_basis`, after both `Mesh`
and `Point` exist.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object containing geometry.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`) with matching `T1,T2,D`.
- `mat::NamedTuple`: Raw material constants (`setup_matconst`) — transient, used only to
  build the per-particle `cmp::Vector{<:AbstractConstitutiveModel}` here (via `setup_cmp`),
  not persisted or stored on `Point` itself.
- `geom::NamedTuple=()`: (Optional) Geometry definition (e.g., number of intervals, number of material points, geometry struct).

# Returns
- `Point`: Material point data structure with all fields initialized for simulation.

# Example
```julia
mpts = setup_mpts(mesh, solver, mat; geom=get_slump(mesh, mat, solver))
println(mpts.nmp)  # Number of material points
```

# Notes
- Initializes material point positions, volumes, densities, and all state variables.
- Sets up phase properties (solid and liquid).
- Handles both 2D and 3D cases.
"""
function setup_mpts(mesh::Mesh{T1,T2,D},solver::S,mat::NamedTuple; geom::NamedTuple=(;)) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
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
    # constructor - create components
    if solver.strain.deform == "finite"
        elast = FiniteElasticity(T1, T2, nmp, D)
    elseif solver.strain.deform == "infinitesimal"
        elast = LinearElasticity(T1, T2, nmp, D)
    end

    rheo = DruckerPragerRheology{T1,T2,D}(
        T2.(vec(copy(geom.coh0))),  # c₀
        T2.(vec(copy(geom.cohr))),  # cᵣ
        T2.(vec(copy(geom.phi))),   # ϕ
        T2.(zeros(nmp)),            # Δλ
        T2.(zeros(2,nmp)),          # ϵpII
        T2.(zeros(nmp))             # ϵpV
    )

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

    cmp = setup_cmp(nmp,T2.(vec(copy(geom.coh0))),T2.(vec(copy(geom.cohr))),T2.(vec(copy(geom.phi))); E=T2(mat[:E]),ν=T2(mat[:ν]),Hp=T2(mat[:Hp]),D=Int(D))
    CM  = eltype(cmp)

    s = PointSolidPhase{T1,T2,D,typeof(elast),typeof(rheo),CM,TM,TV,TS}(
        [zero(TV)  for _ in 1:nmp]                         , # u
        [TV(T2.(vp[:,p])) for p in 1:nmp]                  , # v
        # mechanical properties
        T2.(vec(copy(ρ0)))                                  , # ρ₀
        T2.(vec(copy(ρ0)))                                  , # ρ
        T2.(vec(copy(geom.coh0)))                           , # c₀
        T2.(vec(copy(geom.cohr)))                           , # cᵣ
        T2.(vec(copy(geom.phi)))                            , # ϕ₀
        T2.(zeros(nmp))                                     , # Δλ
        T2.(zeros(2,nmp))                                   , # ϵpII
        T2.(zeros(nmp))                                     , # ϵpV
        # tensor in voigt notation
        [zero(TS) for _ in 1:nmp]                          , # σᵢ
        [zero(TS) for _ in 1:nmp]                          , # σn
        [zero(TS) for _ in 1:nmp]                          , # τᵢ
        T2.(zeros(nmp))                                     , # P
        # tensor in matrix notation
        [zero(TM) for _ in 1:nmp]                          , # ∇vᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ∇uᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ΔFᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fn
        [zero(TM) for _ in 1:nmp]                          , # ϵᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ϵn
        [zero(TM) for _ in 1:nmp]                          , # ωᵢⱼ
        # component-based fields
        elast                                               , # elast::E
        rheo                                               , # rheo::R
        cmp                                                  , # cmp::Vector{CM}
    )
    #=
    t = PointThermalPhase{T1,T2,D}(
        T2.(vec(copy(geom.c)))                            , # c::Vector{T2} specific heat capacity vector
        T2.(vec(copy(geom.k)))                             , # k::Vector{T2} thermal conductivity vector
        T2.(zeros(D,nmp))                                  , # q::Matrix{T2} heat flux array
        T2.(vec(copy(geom.T)))                            , # T::Vector{T2} temperature vector
    )
    =#

    mpts = Point{T1,T2,D,typeof(elast),typeof(rheo),CM,TM,TV,TS}(
        # general information
        T1(D)                              , # ndim
        T1(nmp)                              , # nmp
        T2.(zeros(D))                      , # vmax
        # basis-related quantities
        T2.(zeros(props.nn,D,nmp        ))  , # Δnp
        # APIC-related
        T2.(zeros(D,D,nmp            ))  , # Bᵢⱼ
        T2.(zeros(D,D,nmp            ))  , # Dᵢⱼ
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
        nothing                              , #
    )
    return mpts
end
