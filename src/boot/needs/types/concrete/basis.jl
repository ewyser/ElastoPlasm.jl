# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export neighbs

struct neighbs{T<:Integer,D<:AbstractDimension} <: AbstractNeighbs
    stencils::NTuple{N,UnitRange{T}} where N
    nn::T
end

# Generic constructor: infers dimension and dispatches
function neighbs(stencils::NTuple{N,UnitRange{T}}) where {T,N}
    D = N == 1 ? OneDimension : N == 2 ? TwoDimension : N == 3 ? ThreeDimension : error("Unsupported dimension")
    nn = T(prod(length.(stencils)))
    neighbs{T,D}(stencils, nn)
end