
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
function setup_mpts(mesh::Mesh{T1,T2,Flag,Meta,D},instr::Instruction{T1,T2,D},cmpr::NamedTuple; geom::NamedTuple=(;)) where {T1,T2,Flag,Meta,D}
    props = mesh.prprt
    # non-dimensional constant                                                   
    if D<:TwoDimension
        nstr = 3 
    elseif D<:ThreeDimension
        nstr = 6 
    end
    # unpack material geometry
    ni,nmp,xp = geom.ni,geom.nmp,geom.xp 
    # scalars & vectors
    n0 = 0.1.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(props.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = fill(cmpr[:ρ0],nmp)
    # constructor - create components
    if instr.fwrk.deform == "finite"
        elast = FiniteElasticity(T1, T2, nmp, props.dim)
    elseif instr.fwrk.deform == "infinitesimal"
        elast = LinearElasticity(T1, T2, nmp, props.dim)
    end
    
    rheo = DruckerPragerRheology{T1,T2}(
        T2.(vec(copy(geom.coh0))),  # c₀
        T2.(vec(copy(geom.cohr))),  # cᵣ
        T2.(vec(copy(geom.phi))),   # ϕ
        T2.(zeros(nmp)),            # Δλ
        T2.(zeros(2,nmp)),          # ϵpII
        T2.(zeros(nmp))             # ϵpV
    )

    s = PointSolidPhase{T1,T2,typeof(elast),typeof(rheo)}(
        T2.(zeros(size(xp)))                               , # u
        T2.(zeros(size(xp)))                               , # v
        # mechanical properties
        T2.(vec(copy(ρ0)))                                 , # ρ₀
        T2.(vec(copy(ρ0)))                                 , # ρ
        T2.(vec(copy(geom.coh0)))                          , # c₀
        T2.(vec(copy(geom.cohr)))                          , # cᵣ
        T2.(vec(copy(geom.phi)))                           , # ϕ
        T2.(zeros(nmp))                                    , # Δλ
        T2.(zeros(2,nmp))                                  , # ϵpII
        T2.(zeros(nmp))                                    , # ϵpV
        # tensor in voigt notation
        T2.(zeros(nstr,nmp))                               , # σᵢ
        T2.(zeros(nstr,nmp))                               , # τᵢ
        # tensor in matrix notation
        T2.(zeros(props.dim,props.dim,nmp))                  , # ∇vᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp))                  , # ∇uᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp))                  , # ΔFᵢⱼ
        T2.(repeat(Matrix(1.0I,props.dim,props.dim),1,1,nmp)), # Fᵢⱼ 
        T2.(repeat(Matrix(1.0I,props.dim,props.dim),1,1,nmp)), # bᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp))                  , # ϵᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp))                  , # ωᵢⱼ
        T2.(zeros(props.dim,props.dim,nmp))                  , # σJᵢⱼ
        # new component-based fields
        elast                                              , # elast::E
        rheo                                               , # rheo::R
    )
    t = PointThermalPhase{T1,T2}(
        T2.(vec(copy(geom.c)))                            , # c::Vector{T2} specific heat capacity vector
        T2.(vec(copy(geom.k)))                             , # k::Vector{T2} thermal conductivity vector
        T2.(zeros(props.dim,nmp))                            , # q::Matrix{T2} heat flux array
        T2.(vec(copy(geom.T)))                            , # T::Vector{T2} temperature vector
    )
    f = PointFluidPhase{T1,T2}(

    )


    mpts = Point{T1,T2,typeof(elast),typeof(rheo)}(
        # general information
        T1(props.dim)                         , # ndim
        T1(nmp)                              , # nmp
        T2.(zeros(props.dim))                 , # vmax
        # basis-related quantities
        T2.(zeros(props.nn,nmp ,props.dim+1))  , # ϕ∂ϕ
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
        # utils
        T2.(Matrix(1.0I,props.dim,props.dim))  , # δᵢⱼ
        # material point properties
        T2.(copy(xp))                        , # x
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
        f                                    , #
        # thermal phase
        t                                    , #
    )
    return mpts 
end
