
"""
    setup_mpts(mesh::Mesh{T1,T2}, cmp::NamedTuple; define::Tuple=(nothing, nothing)) -> Point{T1,T2}

Construct the material point system (mpts) for a simulation, initializing all relevant fields and connectivity.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh object containing geometry and topology.
- `cmp::NamedTuple`: Constitutive model parameters.
- `define::Tuple=(nothing, nothing)`: (Optional) Tuple containing geometry definition (e.g., number of intervals, number of material points, geometry struct).

# Returns
- `Point{T1,T2}`: Material point data structure with all fields initialized for simulation.

# Example
```julia
mpts = setup_mpts(mesh, cmp; define=(ni, nmp, geom))
println(mpts.nmp)  # Number of material points
```

# Notes
- Initializes material point positions, volumes, densities, and all state variables.
- Sets up connectivity arrays and phase properties (solid and liquid).
- Handles both 2D and 3D cases.
"""
function setup_mpts(mesh::Mesh{T1,T2},cmpr::NamedTuple; geom::NamedTuple=(;)) where {T1,T2}
    props = mesh.prprt
    # non-dimensional constant                                                   
    if props.dim == 2 
        nstr = 3 
    elseif props.dim == 3 
        nstr = 6 
    end
    # unpack material geometry
    ni,nmp,xp = geom.ni,geom.nmp,geom.xp 
    # scalars & vectors
    n0 = 0.1.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(props.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = fill(cmpr[:ρ0],nmp)
    m  = cmpr[:ρ0].*v0
    # constructor
    s = MaterialPointSolidPhase{T1,T2}(
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
    )
    f = MaterialPointFluidPhase{T1,T2}(

    )
    t = MaterialPointThermalPhase{T1,T2}(

    )

    mpts = Point{T1,T2}(
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
