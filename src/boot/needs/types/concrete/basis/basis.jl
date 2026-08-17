# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Neighbs, get_basis, Basis

struct Neighbs{T,D,NN} <: AbstractNeighbs
    stencils::NTuple{D,UnitRange{T}}
    nn      ::T
end

# Generic constructor: infers dimension and dispatches
function Neighbs(stencils::NTuple{D,UnitRange{T}}) where {T,D}
    NN = prod(length.(stencils))
    Neighbs{T,D,NN}(stencils, T(NN))
end

# Fuse the neighbor stencil into the basis type: each concrete basis defines its own
# stencil_range(::ConcreteBasis, ::Type{T1}), and this constructor builds the Neighbs from it.
function Neighbs(basis::AbstractBasis, dim::Integer, ::Type{T1}) where {T1}
    Neighbs(ntuple(_ -> stencil_range(basis, T1), dim))
end

"""
    get_basis(which::String) -> AbstractBasis

Map a solver's `basis.which` string to its concrete `AbstractBasis` instance.
"""
function get_basis(which::String)
    if which == "bsmpm"
        BSplineBasis()
    elseif which == "gimpm"
        GimpBasis()
    elseif which == "smpm"
        LinearBasis()
    else
        error("Unsupported basis type: $which")
    end
end

"""
    Basis{T1,D,NN,K<:AbstractBasis}

Owns the topology that links a `Mesh` and a `Point` container: which nodes neighbor each element
(`e2n`) and each material point (`p2n`), which elements neighbor each other (`e2e`) and each
material point (`e2p`), and which points neighbor each other (`p2p`, nonlocal regularization only).
Threaded through the solver as a third argument alongside `mpts`/`mesh`, e.g. `p2e2n(mpts,mesh,basis)`.
"""
struct Basis{T1,D,NN,K<:AbstractBasis}
    kind::K                              # BSplineBasis()/GimpBasis()/LinearBasis()
    e2n  ::Vector{SVector{NN,T1}}
    e2e  ::SparseMatrixCSC{T1,T1}
    p2n  ::Vector{SVector{NN,T1}}
    p2e  ::Vector{T1}
    e2p  ::Matrix{T1}
    p2p  ::Matrix{T1}
end
@adapt_struct Basis