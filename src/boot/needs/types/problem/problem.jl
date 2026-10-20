# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Problem: everything a run needs, built from a Case
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Problem

"""
    Problem(mpts, mesh, basis, time, solver)
    Problem(case::Case) -> Problem

Everything a run needs: the material points `mpts`, the background `mesh`, the `basis` linking them
(connectivity, shape functions, transfer scheme), the `time` stages and the `solver` (configuration
and kernels). Workflows run on it (`elastodynamic!(problem)`, `elastoplasm(problem; workflows)`), and
the `.jld2` of a simulation stores it (`save_problem`/`load_problem`).

`Problem(case)` builds it from a `Case`: the solver (`case_solver_options`), the mesh (`setup_mesh`),
the material points (`case_points`, `setup_mpts`), the time stages (`case.run`) and the basis
(`setup_basis`). The case itself is not kept: everything it describes is in the built objects, and
`setup_simulation` writes it as TOML beside the `.jld2`.
"""
struct Problem{T1,T2,D,P<:Point{T1,T2,D},B<:Basis{T1,T2,D},S<:AbstractSolver{T1,T2,D}}
    mpts   :: P
    mesh   :: Mesh{T1,T2,D}
    basis  :: B
    time   :: Time{T1,T2}
    solver :: S
end
@adapt_struct Problem

function Problem(case::Case)
    L, nel = domain(case)
    solver = get_solver(; dim=length(L), case_solver_options(case)...)
    mesh   = setup_mesh(L, nel, solver; thermal=case.thermal)
    mpts   = setup_mpts(mesh, solver, case_points(case, mesh, solver); thermal=case.thermal)
    time   = setup_time(solver; te=case.run.te, tg=case.run.tg, tep=case.run.tep)
    return Problem(mpts, mesh, setup_basis(mesh, mpts, solver), time, solver)
end
