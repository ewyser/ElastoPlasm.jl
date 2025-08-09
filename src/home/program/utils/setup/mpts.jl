
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
function setup_mpts(mesh::Mesh{T1,T2},cmp::NamedTuple;define::Tuple=(nothing,nothing)) where {T1,T2}
    # non-dimensional constant                                                   
    if mesh.dim == 2 
        nstr = 3 
    elseif mesh.dim == 3 
        nstr = 6 
    end
    # material geometry
    ni,nmp,geom = define
    xp = geom.xp 
    # scalars & vectors
    n0 = 0.0.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(mesh.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = ones(nmp) .* cmp[:ρ0]
    m  = (1.0 .- n0).*cmp[:ρ0].*v0
    # constructor
    mpts = (
        ndim = mesh.dim,
        nmp  = nmp,
        vmax = zeros(mesh.dim),
        x    = copy(xp),
        z₀   = copy(xp[end,:]),
        n₀   = copy(n0),
        n    = copy(n0),
        ℓ₀   = copy(l0), 
        ℓ    = copy(l0),
        Ω₀   = vec(copy(v0)),
        Ω    = vec(copy(v0)),
        s = (;
            u    = zeros(size(xp)), 
            v    = zeros(size(xp)),
            p    = zeros(size(xp)),

            ρ₀   = vec(copy(ρ0)),
            ρ    = vec(copy(ρ0)),
            m    = vec(copy(m)),
            c₀   = vec(copy(geom.coh0)),
            cᵣ   = vec(copy(geom.cohr)),
            ϕ    = vec(copy(geom.phi)),            
            Δλ   = zeros(nmp),
            ϵpII = zeros(2,nmp),
            ϵpV  = zeros(nmp), 
            ΔJ   = ones(nmp),
            J    = ones(nmp),
            # tensor in matrix notation
            δᵢⱼ  = Matrix(1.0I,mesh.dim,mesh.dim), 
            ∇vᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
            ∇uᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
            ΔFᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
            Fᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
            Bᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
            ϵᵢⱼ  = zeros(mesh.dim,mesh.dim,nmp),
            ωᵢⱼ  = zeros(mesh.dim,mesh.dim,nmp),
            σJᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
            # tensor in voigt notation
            σᵢ   = zeros(nstr,nmp),
            τᵢ   = zeros(nstr,nmp),
        ),
        l = (;

        ),
        # additional quantities
        ϕ∂ϕ  = zeros(mesh.nn,nmp ,mesh.dim+1   ),
        δnp  = zeros(mesh.nn,mesh.dim,nmp      ),
        # connectivity
        e2p  = spzeros(Int,nmp,mesh.nel[end]),
        p2p  = spzeros(Int,nmp,nmp),
        p2e  = zeros(Int,nmp),
        p2n  = zeros(Int,mesh.nn,nmp),
    )

    s = Solid{T1,T2}(
        T2.(mpts.s.u)    ,
        T2.(mpts.s.v)    ,
        T2.(mpts.s.p)    ,
        # mechanical properties
        T2.(mpts.s.ρ₀)   ,
        T2.(mpts.s.ρ)    ,
        T2.(mpts.s.m)    ,
        T2.(mpts.s.c₀)   ,
        T2.(mpts.s.cᵣ)   ,
        T2.(mpts.s.ϕ)    ,
        T2.(mpts.s.Δλ)   ,
        T2.(mpts.s.ϵpII) ,
        T2.(mpts.s.ϵpV)  ,
        T2.(mpts.s.ΔJ)   ,
        T2.(mpts.s.J)    ,
        # tensor in voigt notation
        T2.(mpts.s.σᵢ)   ,
        T2.(mpts.s.τᵢ)   ,
        # tensor in matrix notation
        T2.(mpts.s.δᵢⱼ)  ,
        T2.(mpts.s.∇vᵢⱼ) ,
        T2.(mpts.s.∇uᵢⱼ) ,
        T2.(mpts.s.ΔFᵢⱼ) ,
        T2.(mpts.s.Fᵢⱼ)  ,
        T2.(mpts.s.Bᵢⱼ)  ,
        T2.(mpts.s.ϵᵢⱼ)  ,
        T2.(mpts.s.ωᵢⱼ)  ,
        T2.(mpts.s.σJᵢⱼ) ,
    )
    l = Liquid{T1,T2}(

    )
    out = Point{T1,T2}(
        # general information
        T1(mpts.ndim) ,
        T1(mpts.nmp)  ,
        T2.(mpts.vmax) ,
        # basis-related quantities
        T2.(mpts.ϕ∂ϕ)  ,
        T2.(mpts.δnp)  ,
        # connectivity
        T1.(mpts.e2p)  ,
        T1.(mpts.p2p)  ,
        T1.(mpts.p2e)  ,
        T1.(mpts.p2n)  ,
        # material point properties
        T2.(mpts.x)    ,
        T2.(mpts.ℓ₀)   ,
        T2.(mpts.ℓ)    ,
        T2.(mpts.n₀)   ,
        T2.(mpts.n)    ,
        T2.(mpts.Ω₀)   ,
        T2.(mpts.Ω)    ,
        # solid phase
        s       ,
        # liquid phase
        l       ,
    )
    return out 
end
