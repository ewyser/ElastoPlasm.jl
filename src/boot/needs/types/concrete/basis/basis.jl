# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Neighbs

struct Neighbs{T,D,NN} <: AbstractNeighbs
    stencils::NTuple{D,UnitRange{T}}
    nn      ::T
end

# Generic constructor: infers dimension and dispatches
function Neighbs(stencils::NTuple{D,UnitRange{T}}) where {T,D}
    NN = prod(length.(stencils))
    Neighbs{T,D,NN}(stencils, T(NN))
end