# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractBasis end
abstract type AbstractTransfer end

export get_basis, Basis, ∂Nrow

"""
    ∂Nrow(∂N::AbstractArray{T2,3}, nn, p, ::Val{D}) -> SVector{D,T2}

Read neighbor `nn`'s cached shape gradient for particle `p` out of `basis.∂N` (a plain
`(NN,D,nmp)` array — see `Basis` docstring for why it's a plain `Array`, not a
`Vector{SMatrix}`) as an `SVector{D,T2}`, via a `Val`-unrolled `ntuple` so the small
per-read construction stays a zero-allocation stack value.
"""
@inline ∂Nrow(∂N::AbstractArray{T2,3}, nn, p, ::Val{D}) where {T2,D} = SVector{D,T2}(ntuple(d -> ∂N[nn,d,p], Val(D)))

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
    elseif which == "mlsmpm"
        return MLSBasis{T1,D}()
    else
        return error("Unsupported basis type: $which")
    end
end

"""
    Basis{T1,T2,D,NN,K<:AbstractBasis,TR<:AbstractTransfer}

Owns the topology that links a `Mesh` and a `Point` container: which nodes neighbor each element
(`e2n`) and each material point (`p2n`), and which elements neighbor each other (`e2e`, also used
by nonlocal plastic-strain regularization to bound its neighbor search — see `nonlocal.jl`).
Also owns `type`, the per-axis node boundary-layer classification consumed by `BSplineBasis`'s
`eval_basis` — basis-kind-specific data, not mesh/point state.
`N`/`∂N` cache each particle's shape values/gradients for all `NN` neighbor nodes, computed once
per particle per timestep by the `shpfun!` ignite kernel (via the kind-specific `eval_basis`) and
then just read by every P2G/G2P/update kernel, instead of recomputing per call site. They're plain
`(NN,nmp)`/`(NN,D,nmp)` arrays, not `Vector{SVector}`/`Vector{SMatrix}` — bundling all `NN`
per-particle values into one `SVector`/`SMatrix` object defeats Julia's escape analysis for
`NN` much above a handful, causing real per-particle heap allocation; scalar writes into a
pre-allocated plain `Array` do not.
`kind::K` dispatches `eval_basis`/`shpfun!` onto the concrete basis kind
(`BSplineBasis`/`GimpBasis`/`LinearBasis`/`MLSBasis`); `transfer::TR` does the same for
the P2G/G2P transfer scheme (`StdTransfer`/`TpicTransfer`/`ApicTransfer`) — both fields
follow the identical pattern: a config string (`solver.basis.which`/`solver.basis.trsfr`
— both live under the same `basis` config section, mirroring this struct owning both
`kind` and `transfer`) mapped once at setup time (`get_basis`/`get_transfer`) to a
concrete marker instance, then dispatched on via ordinary Julia multiple dispatch
instead of a runtime string branch.
`ApicTransfer` additionally carries `Bᵢⱼ`/`Dᵢⱼ` (`Vector{SMatrix}` per-particle
accumulators, 100% APIC-exclusive — see its own docstring), the same kind-specific-data-
lives-on-the-kind-struct pattern `GimpBasis`'s `stencils` already uses.
Threaded through the solver as a third argument alongside `mpts`/`mesh`, e.g. `p2e2n(mpts,mesh,basis)`.
"""
struct Basis{T1,T2,D,NN,K<:AbstractBasis,TR<:AbstractTransfer}
    kind    ::K                              # BSplineBasis()/GimpBasis()/LinearBasis()/MLSBasis()
    transfer::TR                             # StdTransfer()/TpicTransfer()/ApicTransfer()
    e2n  ::Vector{SVector{NN,T1}}
    e2e  ::SparseMatrixCSC{T1,T1}
    p2n  ::Vector{SVector{NN,T1}}
    p2e  ::Vector{T1}
    type ::Matrix{T1}
    N    ::Matrix{T2}
    ∂N   ::Array{T2,3}
end
@adapt_struct Basis