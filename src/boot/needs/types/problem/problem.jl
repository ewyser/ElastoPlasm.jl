# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Problem: what is solved (a Case) and how (a Solver), built
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Problem

"""
    Problem(mpts, mesh, basis, time, solver)
    Problem(case, solver) -> Problem
    Problem(case) -> Problem

What is solved and how, built: the material points `mpts`, the background `mesh`, the `basis` linking them
(connectivity, shape functions, transfer scheme), the `time` stages and the `solver` (configuration, run settings and
kernels). Workflows run on it (`elastodynamic!(problem)`), and a `Simulation` holds it, ready to run (`execute!`).

`Problem(case, solver)` builds it from both: the mesh from the case's domain and boundary conditions (`setup_mesh`),
the material points from its zones, constants, laws and random field (`case_points`, `setup_mpts`), in the solver's
precision, the time stages from the solver's run settings and the basis from the solver's (`setup_basis`). The case
is not kept: the built objects hold all it describes, and `setup_simulation` writes it as TOML beside the `.jld2`.
`Problem(case)` is `Problem(case, Solver(case))`.

# Fields
- `mpts::P`: the material points, a `Point`.
- `mesh::M`: the background grid, a `Mesh`.
- `basis::B`: the connectivity, shape functions and transfer scheme, a `Basis`.
- `time::Time{T1,T2}`: the time stages.
- `solver::S`: how it is solved, a `Solver`.

`T1`, `T2` and `D` are the index type, float type and dimension, shared by all its parts; `P`, `M`, `B` and `S` are
the parts' types, which `Problem(mpts, mesh, basis, time, solver)` takes from its arguments, requiring them to agree
on `T1`, `T2` and `D` (a `MethodError` otherwise).

# Throws
- `ArgumentError`: when the solver's dimension is not the case's.
"""
struct Problem{T1,T2,D,P,M,B,S}
    mpts   ::P
    mesh   ::M
    basis  ::B
    time   ::Time{T1,T2}
    solver ::S
end
function Problem(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D}, time::Time{T1,T2},
                 solver::Solver{T1,T2,D}) where {T1,T2,D}
    return Problem{T1,T2,D,typeof(mpts),typeof(mesh),typeof(basis),typeof(solver)}(mpts, mesh, basis, time, solver)
end

function Problem(case::Case, solver::Solver{T1,T2,D}) where {T1,T2,D}
    L, nel = case.layout.L, case.layout.nel
    length(L) == D || throw(ArgumentError("a $(D)D solver cannot solve a $(length(L))D case: build one for it, Solver(case)"))
    mesh   = setup_mesh(L, nel, case.bcs, solver; thermal=case.thermal)
    mpts   = setup_mpts(mesh, case_points(case, mesh), case.laws; thermal=case.thermal)
    time   = setup_time(solver; te=solver.run.te, tg=solver.run.tg, tep=solver.run.tep)
    return Problem(mpts, mesh, setup_basis(mesh, mpts, solver), time, solver)
end
function Problem(case::Case)
    return Problem(case, Solver(case))
end
