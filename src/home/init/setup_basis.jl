
"""
    setup_basis(mesh::Mesh{T1,T2,D}, mpts::Point{T1,T2,D}, geom::Geometry{T1,T2,D}, solver::S) -> Basis

Construct the topology container linking `mesh` and `mpts`: element-to-node (`e2n`) and
element-to-element (`e2e`) connectivity from the mesh geometry and neighbor stencil, plus
zero-initialized point-to-node (`p2n`), point-to-element (`p2e`), element-to-points (`e2p`), and
point-to-point (`p2p`) connectivity — populated at runtime by the `p2e2n`/`nonlocal` kernels.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object (see `setup_mesh`).
- `mpts::Point{T1,T2,D}`: Material point object (see `setup_mpts`).
- `geom::Geometry{T1,T2,D}`: Geometry struct, provides the neighbor stencil (`geom.nn`) and element size (`geom.h`).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), selects the basis kind and nonlocal support length scale.

# Returns
- `Basis`: Topology container with fields `kind`, `e2n`, `e2e`, `p2n`, `p2e`, `e2p`, `p2p`.

# Example
```julia
basis = setup_basis(mesh, mpts, geom, solver)
```
"""
function setup_basis(mesh::Mesh{T1,T2,D}, mpts::Point{T1,T2,D}, geom::Geometry{T1,T2,D}, solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    nel, nno = mesh.prprt.nel, mesh.prprt.nno
    NN       = geom.nn.nn
    e2n_data = get_element_to_nodes(nel, nno, geom.nn)
    e2e_data = T1.(e2e(T1(D), nel, geom.h, solver))
    nmp      = mpts.nmp
    kind     = get_basis(solver.basis.which)
    return Basis{T1,D,NN,typeof(kind)}(
        kind,
        e2n_data,
        e2e_data,
        [zero(SVector{NN,T1}) for _ in 1:nmp]  , # p2n
        T1.(zeros(Int,nmp))                    , # p2e
        T1.(spzeros(Int,nmp,nel[end]))         , # e2p
        T1.(spzeros(Int,nmp,nmp))              , # p2p
    )
end
