# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Problem: everything a run needs, built by a Solver from its Case
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Problem

"""
    Problem(mpts, mesh, basis, time, solver)
    Problem(solver::Solver) -> Problem
    Problem(case::Case) -> Problem

Everything a run needs: the material points `mpts`, the background `mesh`, the `basis` linking them
(connectivity, shape functions, transfer scheme), the `time` stages and the `solver` (configuration,
run settings and kernels). Workflows run on it (`elastodynamic!(problem)`), and a `Simulation` holds it,
ready to run (`execute!`).

`Problem(solver)` builds it, for a solver built on a case (`Solver(case)`), from both: the mesh (`setup_mesh`), the
material points (`case_points`, `setup_mpts`), the time stages (`solver.run`) and the basis (`setup_basis`).
`Problem(case)` is `Problem(Solver(case))`. A saved problem's solver drops its case (`stored_solver`): everything it
describes is in the built objects, and `setup_simulation` writes it as TOML beside the `.jld2`.
"""
struct Problem{T1,T2,D,P<:Point{T1,T2,D},B<:Basis{T1,T2,D},S<:AbstractSolver{T1,T2,D},M<:Mesh{T1,T2,D}}
    mpts   :: P
    mesh   :: M
    basis  :: B
    time   :: Time{T1,T2}
    solver :: S
end

function Problem(solver::AbstractSolver)
    case = solver.case
    isnothing(case) && throw(ArgumentError("this solver solves no case (get_solver): build one with Solver(case)"))
    L, nel = case.layout.L, case.layout.nel
    mesh   = setup_mesh(L, nel, solver; thermal=case.thermal)
    mpts   = setup_mpts(mesh, solver, case_points(case, mesh, solver); thermal=case.thermal)
    time   = setup_time(solver; te=solver.run.te, tg=solver.run.tg, tep=solver.run.tep)
    return Problem(mpts, mesh, setup_basis(mesh, mpts, solver), time, solver)
end
function Problem(case::Case)
    return Problem(Solver(case))
end
