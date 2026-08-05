
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
function setup_mesh(geom::Geometry{T1,T2,D,N},instr::Instruction{T1,T2,D}) where {T1,T2,D,N}
    # Mesh & boundary condition setup   
    ndim       = geom.dim                                  
    L,nel,nn,h = geom.L,geom.nel,geom.nn,geom.h 
    xn,nel,nno = get_coords(geom)
    e2n        = get_element_to_nodes(nel, nno, nn)  # precompute e2n for topology kernel
    node_type  = get_node_type(ndim,nno)
    status,xB  = get_bc(xn,instr; ghosts=geom.ghost*h)
    # Constructors
    prop = MeshProperties{T1,T2,D}(
        T1(ndim                         ), # dim
        T1.(nel                         ), # nel
        T1.(nno                         ), # nno
        T1(nn.nn                        ), # nn
        T2.(L                           ), # L
        T2.(h                           ), # h
        T2.(xB                          ), # xB
    )
    bcs = MeshBoundary(
        status
    )
    s = MeshSolidPhase{T1,T2,D}(
        prop,
        bcs,
        T2.(zeros(nno[end]             )), # m
        T2.(zeros(nno[end],nno[end]    )), # mⱼ
        T2.(zeros(ndim,nno[end]        )), # oobf
        T2.(zeros(nno[end]             )), # oobp
        T2.(zeros(ndim,nno[end]        )), # a
        T2.(zeros(ndim,nno[end]        )), # mv
        T2.(zeros(ndim,nno[end]        )), # v
        T2.(zeros(nno[end]             )), # p
    )
    #=
    t = MeshThermalPhase{T1,T2,D}(
        prop,
        bcs,
        T2.(zeros(nno[end]             )), # cᵢ
        T2.(zeros(nno[end]             )), # oobq
        T2.(zeros(nno[end]             )), # Q
        T2.(zeros(nno[end]             )), # mcT
        T2.(zeros(nno[end]             )), # T
    )
    =#
    mesh = Mesh{T1,T2,D}(
        prop,
        # nodal quantities
        T2.(vec(minimum(xn,dims=2)     )), # x₀
        T2.(xn                          ), # x
        T1.(node_type                   ), # node
        T2.(zeros(nno[end]             )), # ΔJ
        # solid phase
        s                                , # solid phase
        # thermal phase
        nothing                          , # thermal phase
        # connectivity
        T1.(e2n                        ) , # e2n
        T1.(e2e(ndim,nel,h,instr       )), # e2e
    )
    return mesh
end