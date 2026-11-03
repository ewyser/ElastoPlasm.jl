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
    Basis{T1,T2,D,NN,K,TR,…}

What links a `Mesh` and a `Point`: all the connectivity, the shape functions' kind and cached values, and the transfer
scheme. Kernels take it after the points and the mesh, `(mpts, mesh, basis, …)`.

# Fields
- `kind::K`: the basis kind, `BSplineBasis`, `GimpBasis`, `LinearBasis` or `MLSBasis` (`basis.which`, `get_basis`):
  `eval_basis`, hence `shpfun!`, dispatch on it.
- `transfer::TR`: the transfer scheme, `StdTransfer`, `TpicTransfer` or `ApicTransfer` (`basis.trsfr`,
  `get_transfer`): `p2n!` and `Bij` dispatch on it; APIC's carries its per-point matrices.
- `e2n::VN`: each element's nodes, an `SVector{NN,T1}` per element.
- `e2e::A`: the elements near each element, an `Adjacency` (also the non-local regularisation's neighbour search).
- `p2n::VN`: each point's `NN` nodes (0 where a node is outside the grid), updated every step (`p2e2n`).
- `p2e::VI`: each point's element.
- `type::MI`: the nodes' boundary layer per axis (`D × nno`), which `BSplineBasis`'s functions depend on.
- `N::M`: the shape values, `NN × nmp`, computed once per step (`shpfun!`) and read by every transfer.
- `∂N::A3`: the shape gradients, `NN × D × nmp`; read a row with `∂Nrow`.

`T1`, `T2` and `D` are the index type, float type and dimension; `NN` is the number of nodes a point reaches
(`stencil_size`), known once the mesh and the kind are; `K` and `TR` are the kind's and the scheme's types, what code
dispatches on. The rest are the arrays' types, host or device alike.

`N` and `∂N` are plain arrays, not vectors of `SVector`s or `SMatrix`es: bundling a point's `NN` values into one static
object defeats Julia's escape analysis for `NN` above a handful and allocates per point (see `gotchas.md`).

# Constructors
    Basis{T1,T2,D,NN}(kind, transfer, e2n, e2e, p2n, p2e, type, N, ∂N)

`setup_basis` builds the basis of a mesh and its points, from the solver's `basis` options.

See also: `Mesh`, `Point`, `get_basis`, `get_transfer`.
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
    Adjacency{T1,V}

A sparse adjacency list in compressed-column form, as a `SparseMatrixCSC`'s pattern without its values: two plain
index vectors, so it lives on a device as well as on the host. `basis.e2e` is one (the elements near each element).

# Fields
- `ptr::V`: where each column starts in `idx`: the entries of column `j` are `idx[ptr[j]:ptr[j+1]-1]`.
- `idx::V`: the entries, column after column.

`T1` is the index type, `V` the vectors' type. Read it as a sparse matrix's pattern, with `nzrange(a, j)` and
`rowvals(a)`.

# Constructors
    Adjacency(ptr, idx)
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