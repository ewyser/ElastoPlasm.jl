export setup_problem

"""
    setup_problem(geom, solver, ic_points; te=0.0, tg=0.0, tep=0.0) -> MechanicalProblem

Thin wrapper replacing the `setup_mesh` → `setup_material_constants` → `setup_mpts` →
`setup_time` sequence with a single call returning a `MechanicalProblem` bundling
`mesh`, `mpts` and `time` — the expensive, IC-defining part of a simulation, kept
independent of `Basis`/`Solver` (see `MechanicalProblem`'s docstring). Material
constants (`mat`, see `setup_material_constants`) are built internally rather than taken as an
argument, since `setup_material_constants` only ever needs `solver` (for its `T1,T2,D` type
parameters) — callers no longer build/pass a separate `mat`. `ic_points` is the
problem-specific material-point generator (e.g. `get_slump`), called as
`ic_points(mesh, mat, solver)` to build the `geom` kwarg `setup_mpts` expects —
this mirrors what `slump_problem` already did by hand, just bundled into one call so
a caller no longer has to thread `mesh`/`mpts`/`time` through three separate
`setup_*` calls itself.

# Arguments
- `geom::Geometry`: Domain geometry (element counts, spacing, ...).
- `solver::S`: Solver instance, used to infer numeric types and dispatch.
- `ic_points::Function`: Problem-specific point generator, `(mesh,mat,solver) -> NamedTuple`.
- `te`,`tg`,`tep`: (Optional) time-phase durations, forwarded to `setup_time`.

# Returns
- `MechanicalProblem`: Bundles `mesh::Mesh`, `mpts::Point`, `time::Time`.

# Example
```julia
problem = setup_problem(geom, solver, get_slump; te=10.0, tg=10.0, tep=5.0)
basis   = setup_basis(problem, solver)
```
"""
function setup_problem(geom::Geometry{T1,T2,D}, solver::S, ic_points::Function; te=0.0,tg=0.0,tep=0.0) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    mesh = setup_mesh(geom, solver)
    mat  = setup_material_constants(solver)
    mpts = setup_mpts(mesh, solver, mat; geom = ic_points(mesh, mat, solver))
    time = setup_time(solver; te=te, tg=tg, tep=tep)
    return MechanicalProblem(mesh, mpts, time)
end
