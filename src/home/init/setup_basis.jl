
"""
    setup_basis(problem, solver) -> Basis

The connectivity and shape-function state linking `problem.mesh` and `problem.mpts`: the basis kind
(`solver.basis.which`) and transfer scheme (`solver.basis.trsfr`), element-to-node (`e2n`) and
element-to-element (`e2e`) connectivity, the nodes' boundary layers (`type`), and the per-point
node and element indices (`p2n`, `p2e`) and shape values (`N`, `∂N`), zero here and filled every
step by the `p2e2n` and `shpfun!` kernels.
"""
function setup_basis(problem::MechanicalProblem{T1,T2,D}, solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    mesh, mpts = problem.mesh, problem.mpts
    nel, nno = mesh.prprt.nel, mesh.prprt.nno
    kind     = get_basis(solver.basis.which, T1, D)
    NN       = prod(length.(kind.stencils))
    e2n_data = get_element_to_nodes(nel, nno, kind.stencils)
    e2e_data = e2e(nel, collect(mesh.prprt.h), solver)
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
