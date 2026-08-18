
"""
    setup_mesh(geom::Geometry{T1,T2,D}, solver::S) -> Mesh

Construct the mesh and associated nodal data structures for a simulation domain. Connectivity
(`e2n`, `e2e`) is built separately by `setup_basis`, after both `Mesh` and `Point` exist.

# Arguments
- `geom::Geometry{T1,T2,D}`: Geometry struct (see `setup_geometry`/`Geometry`).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), must include basis and boundary condition settings.

# Returns
- `Mesh`: Mesh object containing geometry, nodal quantities, and boundary conditions.

# Example
```julia
mesh = setup_mesh(geom, solver)
println(mesh.prprt.nel)
```

# Notes
- Sets up mesh geometry, nodal coordinates, and boundary conditions.
- Initializes nodal quantities (mass, force, acceleration, etc.).
- Handles ghost nodes if required by the basis.
"""
function setup_mesh(geom::Geometry{T1,T2,D},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    # Mesh & boundary condition setup
    L,nel,h = geom.L,geom.nel,geom.h
    xn,nel,nno = get_coords(geom)
    status,xB  = get_bc(xn,solver; ghosts=geom.ghost*h)
    # static array type for node coords, velocity, acceleration
    NN  = prod(length.(get_basis(solver.basis.which, T1, D).stencils))
    DIM = D
    SVD = SVector{DIM,T2}
    nn_total = nno[end]

    prop = MeshProperties{T1,T2,DIM}(
        T1.(nel                         ), # nel
        T1.(nno                         ), # nno
        T1(NN                           ), # nn
        T2.(L                           ), # L
        SVD(T2.(h)                      ), # h
        T2.(xB                          ), # xB
    )
    bcs = MeshBoundary(
        status
    )
    s = MeshSolidPhase{T1,T2,DIM}(
        prop,
        bcs,
        T2.(zeros(nn_total             )), # m
        T2.(zeros(nn_total,nn_total    )), # Mᵢⱼ
        T2.(zeros(DIM,nn_total        )), # oobf  (Matrix: atomic writes)
        T2.(zeros(nn_total             )), # oobp
        [zero(SVD) for _ in 1:nn_total]  , # a
        T2.(zeros(DIM,nn_total        )), # mv    (Matrix: atomic writes)
        [zero(SVD) for _ in 1:nn_total]  , # v
        T2.(zeros(nn_total             )), # p
    )
    mesh = Mesh{T1,T2,DIM}(
        prop,
        # nodal quantities
        SVD(T2.(vec(minimum(xn,dims=2))) )     , # x₀
        [SVD(T2.(xn[:,i])) for i in 1:nn_total], # x
        T2.(zeros(nn_total                  )), # ΔJ
        # solid phase
        s                                , # solid phase
        # thermal phase
        nothing                          , # thermal phase
    )
    return mesh
end
