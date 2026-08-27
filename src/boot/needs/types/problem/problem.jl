# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Problem Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractProblem{T1,T2,D,SP<:AbstractMaterialPointPhase{T1,T2}} end

export MechanicalProblem

"""
    MechanicalProblem{T1,T2,D,CM,TM,TV,TS,ST,SC,SK} <: AbstractProblem{T1,T2,D,SP}

Bundles the IC-defining, expensive-to-regenerate part of a simulation — `Mesh`,
`Point` and `Time` — as one object, distinct from `Basis` (cheap, swappable,
still geometry-derived, built *against* a `Problem`) and `Solver` (config +
kernels). `SP` (the phase type parameter on `AbstractProblem`) is tied to
`PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}`, i.e. the type of `mpts.s` — reusing
the existing `MaterialPointPhase` hierarchy rather than introducing a second,
disconnected "phase" concept at the problem level.

`Basis`/`Solver` are deliberately *not* part of `Problem` — `basis.which` and
`basis.trsfr` (which basis kind, which P2G/G2P transfer scheme — sibling config
knobs, but genuinely independent axes: any basis kind pairs with any transfer
scheme) vary independently in practice, so bundling them into `Problem` would
model a relationship that doesn't exist. `Problem` is also, for now, only an outer-API/
persistence convenience: kernel/workflow signatures still take `mesh`/`mpts`/
`time` as separate positional args, not `Problem` itself — threading `Problem`
through every kernel call site is a much larger, currently out-of-scope lift.
"""
struct MechanicalProblem{T1,T2,D,CM<:AbstractConstitutiveModel{T2,D},TM,TV,TS,ST<:AbstractStrain,SC<:AbstractStress,SK<:AbstractStress} <: AbstractProblem{T1,T2,D,PointSolidPhase{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}}
    mesh ::Mesh{T1,T2,D}
    mpts ::Point{T1,T2,D,CM,TM,TV,TS,ST,SC,SK}
    time ::Time{T1,T2}
end
@adapt_struct MechanicalProblem
