"""
    setup_mesh(nel::Vector{T1}, L::Vector{T2}, instr) -> Mesh

Construct the mesh and associated nodal and topological data structures for a simulation domain.

# Arguments
- `nel::Vector{T1}`: Number of elements in each spatial direction.
- `L::Vector{T2}`: Domain size in each spatial direction.
- `instr`: Simulation instruction dictionary (must include basis and boundary condition settings).

# Returns
- `Mesh`: Mesh object containing geometry, nodal quantities, topology, and boundary conditions.

# Example
```julia
mesh = setup_mesh([40, 10], [64.0, 16.0], instr)
println(mesh.nel)
```

# Notes
- Sets up mesh geometry, nodal coordinates, and boundary conditions.
- Initializes nodal quantities (mass, force, acceleration, etc.) and mesh-to-node topology.
- Handles ghost nodes if required by the basis.
"""
function setup_mesh(nel::Vector{T1},L::Vector{T2},instr) where {T1,T2}
    if instr[:basis][:ghost]
        buffer = T2(2.0)
    else
        buffer = T2(0.0)
    end
    # mesh & boundary conditions                                      
    ndim,nn,h  = get_geom(nel,L,instr)
    xn,nel,nno = get_coords(ndim,L,h; ghosts=buffer.*h)
    status,xB  = get_bc(xn,nno,ndim ; ghosts=buffer.*h, set=instr[:bcs].dirichlet)
    # constructor
    mesh = (;
        dim  = ndim,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h, # mᵢ Mᵢⱼ
        # nodal quantities
        x₀   = vec(minimum(xn,dims=2)     ),
        x    = xn                          ,
        mᵢ   = zeros(nno[end]             ), # lumped mass vector
        Mᵢⱼ  = zeros(nno[end],nno[end]    ), # consistent mass matrix
        oobf = zeros(ndim,nno[end]        ),
        D    = zeros(ndim,nno[end]        ),
        f    = zeros(ndim,nno[end]        ),
        a    = zeros(ndim,nno[end]        ),
        p    = zeros(ndim,nno[end]        ),
        v    = zeros(ndim,nno[end]        ),
        Δu   = zeros(ndim,nno[end]        ),
        ΔJ   = zeros(ndim,nno[end]        ),
        bij  = zeros(ndim,ndim,nno[end]   ),
        # mesh-to-node topology
        e2n  = e2n(ndim,nno,nel,nn        ),
        e2e  = e2e(ndim,nel,h,instr),
        xB   = xB                          ,
    )
    bcs = Boundary{Bool}(
        status
    )
    out = Mesh{T1,T2,Bool,NamedTuple}(
        T1(mesh.dim), 
        T1.(mesh.nel), 
        T1.(mesh.nno), 
        T1(mesh.nn), 
        T2.(mesh.L), 
        T2.(mesh.h), 
        # nodal quantities
        T2.(mesh.x₀), 
        T2.(mesh.x), 
        T2.(mesh.mᵢ), 
        T2.(mesh.Mᵢⱼ),
        T2.(mesh.oobf), 
        T2.(mesh.D), 
        T2.(mesh.f), 
        T2.(mesh.a), 
        T2.(mesh.p), 
        T2.(mesh.v), 
        T2.(mesh.Δu), 
        T2.(mesh.ΔJ), 
        T2.(mesh.bij),
        # mesh-to-node topology
        T1.(mesh.e2n), 
        T1.(mesh.e2e), 
        T2.(mesh.xB), 
        # mesh boundary conditions
        bcs
    )
    return out
end