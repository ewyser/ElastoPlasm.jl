
"""
    setup_basis(problem::MechanicalProblem{T1,T2,D}, solver::S) -> Basis

Construct the topology container linking `problem.mesh` and `problem.mpts`: element-to-node
(`e2n`) and element-to-element (`e2e`) connectivity from the mesh geometry and neighbor stencil,
plus zero-initialized point-to-node (`p2n`) and point-to-element (`p2e`) connectivity —
populated at runtime by the `p2e2n` kernel — and `type`, the per-axis node boundary-layer
classification used by `BSplineBasis`'s `eval_basis`.
`N`/`∂N` (per-particle cached shape values/gradients) are zero-initialized here and populated each
timestep by the `shpfun!` ignite kernel.

Takes `problem` rather than loose `mesh`/`mpts`/`geom`: all the geometric information
`setup_basis` needs (element size, counts) is already on `problem.mesh.prprt` — a
separate `geom::Geometry` argument was redundant (`geom.h` and `mesh.prprt.h` are the
same element-size data), and `mesh`+`mpts` were always the pair a `Problem` already
bundles.

# Arguments
- `problem::MechanicalProblem{T1,T2,D}`: Bundles `mesh`/`mpts` (see `setup_problem`); node
  count (`NN`) comes from the concrete basis kind's stencils.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), selects the basis kind and nonlocal support length scale.

# Returns
- `Basis`: Topology container with fields `kind`, `e2n`, `e2e`, `p2n`, `p2e`.

# Example
```julia
basis = setup_basis(problem, solver)
```
"""
function setup_basis(problem::MechanicalProblem{T1,T2,D}, solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    mesh, mpts = problem.mesh, problem.mpts
    nel, nno = mesh.prprt.nel, mesh.prprt.nno
    kind     = get_basis(solver.basis.which, T1, D)
    NN       = prod(length.(kind.stencils))
    e2n_data = get_element_to_nodes(nel, nno, kind.stencils)
    e2e_data = T1.(e2e(T1(D), nel, collect(mesh.prprt.h), solver))
    nmp      = mpts.nmp
    transfer = get_transfer(solver.basis.trsfr, T2, D, nmp)
    return Basis{T1,T2,D,NN,typeof(kind),typeof(transfer)}(
        kind,
        transfer,
        e2n_data,
        e2e_data,
        [zero(SVector{NN,T1}) for _ in 1:nmp]  , # p2n
        T1.(zeros(Int,nmp))                    , # p2e
        T1.(get_node_type(T1(D),nno))          , # type
        zeros(T2,NN,nmp)                       , # N
        zeros(T2,NN,D,nmp)                     , # ∂N
    )
end
