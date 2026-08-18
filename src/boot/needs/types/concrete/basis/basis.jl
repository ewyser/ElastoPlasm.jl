# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export get_basis, Basis

"""
    get_basis(which::String, ::Type{T1}, dim::Integer) -> AbstractBasis

Map a solver's `basis.which` string to its concrete `AbstractBasis` instance, parametrized
by index type `T1` and problem dimension `dim`.
"""
function get_basis(which::String, ::Type{T1}, dim::Integer) where {T1}
    if which == "bsmpm"
        BSplineBasis{T1,dim}()
    elseif which == "gimpm"
        GimpBasis{T1,dim}()
    elseif which == "smpm"
        LinearBasis{T1,dim}()
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