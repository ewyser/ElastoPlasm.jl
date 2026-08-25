# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Problem Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export MechanicalProblem

"""
    MechanicalProblem{T1,T2,D,E,CM,TM,TV,TS} <: AbstractProblem{T1,T2,D,SP}

Bundles the IC-defining, expensive-to-regenerate part of a simulation — `Mesh`,
`Point` and `Time` — as one object, distinct from `Basis` (cheap, swappable,
still geometry-derived, built *against* a `Problem`) and `Solver` (config +
kernels). `SP` (the phase type parameter on `AbstractProblem`) is tied to
`PointSolidPhase{T1,T2,D,E,CM,TM,TV,TS}`, i.e. the type of `mpts.s` — reusing
the existing `MaterialPointPhase` hierarchy rather than introducing a second,
disconnected "phase" concept at the problem level.

`Basis`/`Solver` are deliberately *not* part of `Problem` — `basis.which` and
`transfer.trsfr` vary independently in practice, so bundling them would model a
relationship that doesn't exist. `Problem` is also, for now, only an outer-API/
persistence convenience: kernel/workflow signatures still take `mesh`/`mpts`/
`time` as separate positional args, not `Problem` itself — threading `Problem`
through every kernel call site is a much larger, currently out-of-scope lift.
"""
struct MechanicalProblem{T1,T2,D,E<:AbstractElasticity{T1,T2},CM<:AbstractConstitutiveModel{T2,D},TM,TV,TS} <: AbstractProblem{T1,T2,D,PointSolidPhase{T1,T2,D,E,CM,TM,TV,TS}}
    mesh ::Mesh{T1,T2,D}
    mpts ::Point{T1,T2,D,E,CM,TM,TV,TS}
    time ::Time{T1,T2}
end
@adapt_struct MechanicalProblem
