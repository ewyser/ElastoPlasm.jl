# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export get_basis, Basis

"""
    get_basis(which::String, ::Type{T1}, dim::Integer) -> AbstractBasis

Map a solver's `basis.which` string to its concrete `AbstractBasis` instance, parametrized
by index type `T1` and problem dimension `dim`.
"""
function get_basis(which::String, ::Type{T1}, D::Integer) where {T1}
    if which == "bsmpm"
        return BSplineBasis{T1,D}()
    elseif which == "gimpm"
        return GimpBasis{T1,D}()
    elseif which == "smpm"
        return LinearBasis{T1,D}()
    else
        return error("Unsupported basis type: $which")
    end
end

"""
    Basis{T1,D,NN,K<:AbstractBasis}

Owns the topology that links a `Mesh` and a `Point` container: which nodes neighbor each element
(`e2n`) and each material point (`p2n`), which elements neighbor each other (`e2e`) and each
material point (`e2p`), and which points neighbor each other (`p2p`, nonlocal regularization only).
Also owns `type`, the per-axis node boundary-layer classification consumed by `BSplineBasis`'s
`eval_basis` — basis-kind-specific data, not mesh/point state.
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
    type ::Matrix{T1}
end
@adapt_struct Basis