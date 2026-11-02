# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete Basis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractBasisKind end
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
    get_basis(which::String, ::Type{T1}, dim::Integer; domain="stretch") -> AbstractBasisKind

Map a solver's `basis.which` string to its concrete `AbstractBasisKind` instance, parametrized
by index type `T1` and problem dimension `dim`; GIMP's carries its domain update (`basis.domain`,
`get_domain`).
"""
function get_basis(which::String, ::Type{T1}, D::Integer; domain::String="stretch") where {T1}
    if which == "bsmpm"
        return BSplineBasis{T1,D}()
    elseif which == "gimpm"
        return GimpBasis{T1,D}(get_domain(domain))
    elseif which == "smpm"
        return LinearBasis{T1,D}()
    elseif which == "mlsmpm"
        return MLSBasis{T1,D}()
    else
        return error("Unsupported basis type: $which")
    end
end

"""
    Basis{T1,T2,D,NN,K<:AbstractBasisKind,TR<:AbstractTransfer}

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
struct Basis{T1,T2,D,NN,K<:AbstractBasisKind,TR<:AbstractTransfer,VN,A,VI,MI,M,A3}
    kind    ::K                              # BSplineBasis()/GimpBasis()/LinearBasis()/MLSBasis()
    transfer::TR                             # StdTransfer()/TpicTransfer()/ApicTransfer()
    e2n  ::VN
    e2e  ::A                                 # an Adjacency{T1}
    p2n  ::VN
    p2e  ::VI
    type ::MI
    N    ::M
    ∂N   ::A3
end
function Basis{T1,T2,D,NN}(kind::K, transfer::TR, e2n::VN, e2e::A, p2n::VN, p2e::VI, type::MI, N::M,
                           ∂N::A3) where {T1,T2,D,NN,K,TR,VN,A,VI,MI,M,A3}
    return Basis{T1,T2,D,NN,K,TR,VN,A,VI,MI,M,A3}(kind, transfer, e2n, e2e, p2n, p2e, type, N, ∂N)
end
function Adapt.adapt_structure(to, basis::Basis{T1,T2,D,NN}) where {T1,T2,D,NN}
    return adapt_fields(to, basis, Basis{T1,T2,D,NN})
end

"""
    Adjacency{T1}(ptr, idx)

A sparse adjacency list in compressed-column form, as a `SparseMatrixCSC`'s pattern without its values: the entries
of column `j` are `idx[ptr[j]:ptr[j+1]-1]`, read with `nzrange(a, j)` and `rowvals(a)` like a sparse matrix's. Two
plain index vectors, so it lives on a device as well as on the host. `basis.e2e` is one (the elements near each
element).
"""
struct Adjacency{T1,V}
    ptr ::V
    idx ::V
end
function Adjacency(ptr::V, idx::V) where {T1,V<:AbstractVector{T1}}
    return Adjacency{T1,V}(ptr, idx)
end
@adapt_struct Adjacency
@inline function SparseArrays.nzrange(a::Adjacency, j::Integer)
    return a.ptr[j]:(a.ptr[j+1]-1)
end
@inline function SparseArrays.rowvals(a::Adjacency)
    return a.idx
end

"""
    stencil_size(basis) -> Int

The number of nodes a material point's shape functions reach (`NN`, the length of each
`basis.p2n[p]`): 4/9/16 in 2D for linear/GIMP/cubic B-spline bases. A type parameter, so loops
`for nn ∈ 1:stencil_size(basis)` have a compile-time bound.
"""
stencil_size(::Basis{T1,T2,D,NN}) where {T1,T2,D,NN} = NN