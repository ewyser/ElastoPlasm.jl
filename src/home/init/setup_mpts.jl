
"""
    setup_mpts(mesh::Mesh{T1,T2,D}, instr::Instruction{T1,T2,D}, cmpr::NamedTuple; geom::NamedTuple=()) -> Point{T1,T2}

Construct the material point system (mpts) for a simulation, initializing all relevant fields and connectivity.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object containing geometry and topology.
- `instr::Instruction{T1,T2,D}`: Simulation instruction structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `geom::NamedTuple=()`: (Optional) Geometry definition (e.g., number of intervals, number of material points, geometry struct).

# Returns
- `Point{T1,T2}`: Material point data structure with all fields initialized for simulation.

# Example
```julia
mpts = setup_mpts(mesh, instr, cmpr; geom=get_slump(mesh, cmpr, instr))
println(mpts.nmp)  # Number of material points
```

# Notes
- Initializes material point positions, volumes, densities, and all state variables.
- Sets up connectivity arrays and phase properties (solid and liquid).
- Handles both 2D and 3D cases.
"""
function setup_mpts(mesh::Mesh{T1,T2,D},instr::Instruction{T1,T2,D},cmpr::NamedTuple; geom::NamedTuple=(;)) where {T1,T2,D}
    props = mesh.prprt
    # non-dimensional constant
    if D == TwoDimension
        nstr = 3
    elseif D == ThreeDimension
        nstr = 6
    else
        error("Unsupported dimension type: $D")
    end
    # unpack material geometry
    ni,nmp,xp = geom.ni,geom.nmp,geom.xp 
    # scalars & vectors
    n0 = 0.1.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(props.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = fill(cmpr[:ρ0],nmp)
    # initial velocity (if provided)
    vp = haskey(geom, :vp) ? geom.vp : zeros(size(xp))
    #vp = zeros(size(xp))
    # constructor - create components
    if instr.fwrk.deform == "finite"
        elast = FiniteElasticity(T1, T2, nmp, props.dim)
    elseif instr.fwrk.deform == "infinitesimal"
        elast = LinearElasticity(T1, T2, nmp, props.dim)
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
    if D == TwoDimension
        TM = SMatrix{2,2,T2,4}
        TV = SVector{2,T2}
        TS = SVector{3,T2}
    else
        TM = SMatrix{3,3,T2,9}
        TV = SVector{3,T2}
        TS = SVector{6,T2}
    end

    s = PointSolidPhase{T1,T2,D,typeof(elast),typeof(rheo),TM,TV,TS}(
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
        [TM(I)    for _ in 1:nmp]                          , # bᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # bn
        [zero(TM) for _ in 1:nmp]                          , # ϵᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ϵn
        [zero(TM) for _ in 1:nmp]                          , # ωᵢⱼ
        # component-based fields
        elast                                               , # elast::E
        rheo                                               , # rheo::R
    )
    #=
    t = PointThermalPhase{T1,T2,D}(
        T2.(vec(copy(geom.c)))                            , # c::Vector{T2} specific heat capacity vector
        T2.(vec(copy(geom.k)))                             , # k::Vector{T2} thermal conductivity vector
        T2.(zeros(props.dim,nmp))                            , # q::Matrix{T2} heat flux array
        T2.(vec(copy(geom.T)))                            , # T::Vector{T2} temperature vector
    )
    =#

    basis = if instr.basis.which == "bsmpm"
        BSplineBasis()
    elseif instr.basis.which == "gimpm"
        GimpBasis()
    elseif instr.basis.which == "smpm"
        LinearBasis()
    end
    mpts = Point{T1,T2,D,typeof(basis),typeof(elast),typeof(rheo),TM,TV,TS}(
        # basis
        basis                                , # basis
        # general information
        T1(props.dim)                         , # ndim
        T1(nmp)                              , # nmp
        T2.(zeros(props.dim))                 , # vmax
        # basis-related quantities
        T2.(zeros(props.nn,props.dim,nmp   ))  , # Δnp
        # APIC-related
        T2.(zeros(props.dim,props.dim,nmp  ))  , # Bᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp  ))  , # Dᵢⱼ  
        # connectivity
        T1(props.nn)                          , # nn
        T1.(spzeros(Int,nmp,props.nel[end]))  , # e2p
        T1.(spzeros(Int,nmp,nmp          ))  , # p2p
        T1.(zeros(Int,nmp                ))  , # p2e
        T1.(zeros(Int,props.nn,nmp        ))  , # p2n
        # material point properties
        [TV(T2.(xp[:,i])) for i in 1:nmp]           , # x
        T2.(copy(l0))                        , # ℓ₀
        T2.(copy(l0))                        , # ℓ
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
