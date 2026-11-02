# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Solver: how a case is solved
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Solver, ExplicitSolver, ImplicitSolver

# how the equations are solved, `Solver`'s `solution` (the option `solution`, `get_solutions()`)
abstract type AbstractSolution end
"Explicit time integration: `elastodynamic!`, `elastoplastic!`, `thermodynamic!`."
struct ExplicitSolution <: AbstractSolution end
"Quasi-static, by dynamic relaxation: `elastoquasistatic!`."
struct ImplicitSolution <: AbstractSolution end
"""
    get_solutions() -> NamedTuple

The values of the option `solution`, each with its marker: `explicit = ExplicitSolution()`,
`implicit = ImplicitSolution()`: a solver's marker is `get_solutions().implicit`, or `table_value(get_solutions(),
name, "solution")` for a name given as text; `solution_name` reads it back, and `get_option` lists its names, in this
order.
"""
function get_solutions()
    return (; explicit = ExplicitSolution(), implicit = ImplicitSolution())
end

"""
    Solver{T1,T2,D,SL,…}

How a case is solved: the numerical options, the run settings and the kernels built on them. It holds nothing of a
case, so one solver solves several cases of its dimension (`Problem(case, solver)`).

# Fields
- `solution::SL`: the method, `ExplicitSolution()` or `ImplicitSolution()` (the option `solution`, `get_solutions()`);
  the workflows dispatch on it (`ExplicitSolver`, `ImplicitSolver`).
- `dtype::NamedTuple`: the precision, `(; precision, T0, bits)`: `T0 = (T1, T2)`.
- `basis::NamedTuple`: the basis and transfer scheme, `(; which, trsfr, C_pf, domain)`.
- `stab::NamedTuple`: stabilisation, `(; locking, damping, musl, sfrac)`.
- `nonloc::NamedTuple`: non-local regularisation of the plastic strain, `(; status, ls)`.
- `plot::NamedTuple`: live plots, `(; status, freq, dpi, what)`.
- `output::NamedTuple`: the HDF5 output, `(; status, every, what)`.
- `backend::NamedTuple`: where the kernels run, `(; select, exec)`: `exec` is the UnifiedBackend `Platform` selected.
- `run::NamedTuple`: the run settings, `(; te, tg, tep, workflow, lstps)` (`solver_run`).
- `cairn::NamedTuple`: the kernels, `(; ignite, mapsto, update, implicit)`; `nothing` in a saved solver.

`T1` and `T2` are the index and float types (`dtype.T0`), `D` the dimension, `SL` the solution's type; the other
parameters are the fields' types.

# Constructors
    Solver(case; run=(;), overrides…)
    Solver(; dim=2, sections…)

`Solver(case)` reads a case's settings; `Solver(; dim, sections…)` takes `get_default()` with `sections` merged over
it. Both end with `build_solver`.

# Throws
- `ArgumentError`: for physics among the options (`material`, `grf`, `bcs`: a case's), an unknown or invalid option,
  an unknown solution.

See also: `ExplicitSolver`, `ImplicitSolver`, `get_default`, `Problem`.
"""
struct Solver{T1,T2,D,SL<:AbstractSolution,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    solution ::SL
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    nonloc   ::NL
    plot     ::PT
    output   ::OT
    backend  ::BK
    run      ::RN
    cairn    ::CN
end
function Solver{T1,T2,D}(solution::SL,dtype::DT,basis::BA,stab::SB,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,cairn::CN) where {T1,T2,D,SL,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    return Solver{T1,T2,D,SL,DT,BA,SB,NL,PT,OT,BK,RN,CN}(solution,dtype,basis,stab,nonloc,plot,output,backend,run,cairn)
end

"A `Solver` of the explicit solution, what `elastodynamic!`, `elastoplastic!` and `thermodynamic!` take."
const ExplicitSolver{T1,T2,D} = Solver{T1,T2,D,ExplicitSolution}
"A `Solver` of the implicit solution, dynamic relaxation."
const ImplicitSolver{T1,T2,D} = Solver{T1,T2,D,ImplicitSolution}

"""
    solution_name(solver) -> String

A solver's solution as its option value, `"explicit"` or `"implicit"`.
"""
function solution_name(s::Solver)
    return string(findfirst(==(s.solution), get_solutions()))
end
