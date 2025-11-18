
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
function setup_mesh(instr::NamedTuple; geom::NamedTuple=(;))
    # unpack arithmetic precision
    T1,T2 = first(instr[:dtype].T0),last(instr[:dtype].T0) 
    # add ghost points if needed
    if instr[:basis][:ghost]
        buffer = T2(2.0)
    else
        buffer = T2(0.0)
    end
    # mesh & boundary conditions   
    ndim       = geom.ndim                                  
    L,nel,nn,h = geom.L,geom.nel,geom.nn,geom.h
    xn,nel,nno = get_coords(ndim,L,h; ghosts=buffer.*h)
    status,xB  = get_bc(xn,instr; ghosts=buffer.*h)
    # constructor
    prop = MeshProperties{T1,T2}(
        T1(ndim                         ), # dim
        T1.(nel                         ), # nel
        T1.(nno                         ), # nno
        T1(nn                           ), # nn
        T2.(L                           ), # L
        T2.(h                           ), # h
        T2.(xB                          ), # xB
    )
    bcs = MeshBoundary{Bool}(
        status
    )
    s = MeshSolidPhase{T1,T2,Bool}(
        prop,
        bcs,
        T2.(zeros(nno[end]             )), # mᵢ
        T2.(zeros(nno[end],nno[end]    )), # Mᵢⱼ
        T2.(zeros(ndim,nno[end]        )), # oobf
        T2.(zeros(ndim,nno[end]        )), # a
        T2.(zeros(ndim,nno[end]        )), # mv
        T2.(zeros(ndim,nno[end]        )), # v
    )
    t = MeshThermalPhase{T1,T2,Bool}(
        prop,
        bcs,
        T2.(zeros(nno[end]             )), # cᵢ
        T2.(zeros(nno[end]             )), # oobq
        T2.(zeros(nno[end]             )), # Q
        T2.(zeros(nno[end]             )), # mcT
        T2.(zeros(nno[end]             )), # T
    )
    mesh = Mesh{T1,T2,Bool,NamedTuple}(
        prop,
        # nodal quantities
        T2.(vec(minimum(xn,dims=2)     )), # x₀
        T2.(xn                          ), # x
        T2.(zeros(ndim,nno[end]        )), # ΔJ
        # solid phase
        s                                , # solid phase
        # thermal phase
        t                                , # thermal phase
        # connectivity
        T1.(e2n(ndim,nno,nel,nn        )), # e2n
        T1.(e2e(ndim,nel,h,instr       )), # e2e
    )
    return mesh
end