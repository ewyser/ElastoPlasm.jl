# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Solver: how a case is solved
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Solver

"""
    Solver{T1,T2,D,…}

How a case is solved: the numerical options, the run settings and the kernels built on them. It holds nothing of a
case, so one solver solves several cases of its dimension (`Problem(case, solver)`).

# Fields
- `solution::String`: the method, `"explicit"` (time integration: `elastodynamic!`, `elastoplastic!`,
  `thermodynamic!`, which refuse another, `explicit_only`) or `"implicit"` (quasi-static, by dynamic relaxation:
  `elastoquasistatic!`).
- `dtype::NamedTuple`: the precision, `(; precision, T0, bits)`: `T0 = (T1, T2)`.
- `basis::NamedTuple`: the basis and transfer scheme, `(; which, trsfr, C_pf, domain)`.
- `stab::NamedTuple`: stabilisation, `(; locking, damping, musl, sfrac)`.
- `nonloc::NamedTuple`: non-local regularisation of the plastic strain, `(; status, ls)`.
- `plot::NamedTuple`: live plots, `(; status, freq, dpi, what)`.
- `output::NamedTuple`: the HDF5 output, `(; status, every, what)`.
- `backend::NamedTuple`: where the kernels run, `(; select, exec)`: `exec` is the UnifiedBackend `Platform` selected.
- `run::NamedTuple`: the run settings, `(; te, tg, tep, workflow, lstps)` (`solver_run`).
- `cairn::NamedTuple`: the kernels, `(; ignite, mapsto, update, implicit)`; `nothing` in a saved solver.

`T1` and `T2` are the index and float types (`dtype.T0`), `D` the dimension; the other parameters are the fields'
types.

# Constructors
    Solver(case; run=(;), overrides…)
    Solver(; dim=2, sections…)

`Solver(case)` reads a case's settings; `Solver(; dim, sections…)` takes `get_default()` with `sections` merged over
it. Both end with `build_solver`.

# Throws
- `ArgumentError`: for physics among the options (`material`, `grf`, `bcs`: a case's), an unknown or invalid option,
  an unknown solution.

See also: `get_default`, `get_option`, `Problem`.
"""
struct Solver{T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    solution ::String
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
function Solver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,cairn::CN) where {T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    return Solver{T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}(solution,dtype,basis,stab,nonloc,plot,output,backend,run,cairn)
end
