# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mesh Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Mesh

# The mesh, point and basis containers hold their arrays through type parameters (`V`, `VV`, `M`, …, appended after
# the parameters code dispatches on), so the same struct holds host `Array`s or a device's arrays (`adapt`, see
# `adapt_fields` in utils.jl); `Mesh{T1,T2,D}` etc. still name every version. They are left unbounded: bounds such
# as `V<:AbstractVector{T2}` refer to the other parameters, which made every subtype check on these types, hence
# compilation, about three times slower (a first problem built 10 → 34 s); the partial constructors (`Mesh{T1,T2,D}(…)`)
# and the setup functions give them their types.

"""
    MeshProperties{T1,T2,D,N,M}

The grid's sizes, a plain value of static vectors: read as it is on the host (launch sizes, plots) and in kernels, so
it needs no adapt rule.

# Fields
- `nel::SVector{N,T1}`: the elements per axis, then their total (`N = D+1`).
- `nno::SVector{N,T1}`: the nodes per axis, then their total.
- `L::SVector{D,T2}`: the domain's size [m].
- `h::SVector{D,T2}`: the element size per axis [m].
- `xB::SMatrix{D,2,T2,M}`: the domain's bounds, a `[lower upper]` row per axis (`M = 2D`).

# Constructors
    MeshProperties{T1,T2,D}(nel, nno, L, h, xB)

`N` and `M` come from the arguments' lengths.
"""
struct MeshProperties{T1,T2,D,N,M}
    # general information
    nel  ::SVector{N,T1}
    nno  ::SVector{N,T1}
    L    ::SVector{D,T2}
    h    ::SVector{D,T2}
    xB   ::SMatrix{D,2,T2,M}
end
function MeshProperties{T1,T2,D}(nel, nno, L, h, xB) where {T1,T2,D}
    N = length(nel)
    return MeshProperties{T1,T2,D,N,2D}(SVector{N,T1}(nel), SVector{N,T1}(nno), SVector{D,T2}(L), SVector{D,T2}(h),
                                         SMatrix{D,2,T2,2D}(xB))
end

"""
    MeshBoundary{B}

The Dirichlet conditions of the nodes (`get_bc`).

# Fields
- `status::B`: `D × nno`, `true` where a node's motion along that axis is blocked; `B` is the matrix's type (host or
  device).
"""
struct MeshBoundary{B<:AbstractMatrix{Bool}}
    status::B
end
@adapt_struct MeshBoundary

"""
    MeshSolidPhase{T1,T2,D,P,B,V,M,VV}

The solid phase of the mesh: the nodal fields of a step, one entry per node, reset and assembled every step.

# Fields
- `prprt::P`: the grid's sizes, the `MeshProperties`.
- `bcs::B`: the nodes' Dirichlet conditions, a `MeshBoundary`.
- `m::V`: the lumped mass [kg].
- `oobf::M`: the out-of-balance force, `D × nno` (a matrix: atomic writes).
- `a::VV`: the acceleration, an `SVector{D,T2}` per node.
- `mv::M`: the momentum, `D × nno` (a matrix: atomic writes).
- `v::VV`: the velocity (dynamic relaxation: the step's displacement `Δu`).

`T1`, `T2` and `D` are the index type, float type and dimension; `P` and `B` the fields' types; `V`, `M` and `VV` the
arrays' (of `T2`, a matrix of `T2`, of `SVector{D,T2}`; host or device).

# Constructors
    MeshSolidPhase{T1,T2,D}(prprt, bcs, m, oobf, a, mv, v)
"""
struct MeshSolidPhase{T1,T2,D,P,B,V,M,VV}
    prprt ::P
    bcs   ::B
    m     ::V  # consistent lumped mass matrix
    oobf  ::M  # out-of-balance mechanical load (atomic writes: keep Matrix)
    a     ::VV # acceleration per node : SVector{ndim,T2}
    mv    ::M  # momentum (atomic writes: keep Matrix)
    v     ::VV # velocity per node    : SVector{ndim,T2}
end
function MeshSolidPhase{T1,T2,D}(prprt::P, bcs::B, m::V, oobf::M, a::VV, mv::M, v::VV) where {T1,T2,D,P,B,V,M,VV}
    return MeshSolidPhase{T1,T2,D,P,B,V,M,VV}(prprt, bcs, m, oobf, a, mv, v)
end
function Adapt.adapt_structure(to, s::MeshSolidPhase{T1,T2,D}) where {T1,T2,D}
    return adapt_fields(to, s, MeshSolidPhase{T1,T2,D})
end

"""
    MeshThermalPhase{T1,T2,D,P,B,V}

The thermal phase of the mesh, one entry per node (a mesh's `t`, when the case is thermal).

# Fields
- `prprt::P`: the grid's sizes, the `MeshProperties`.
- `bcs::B`: the nodes' Dirichlet conditions, a `MeshBoundary`.
- `cᵢ::V`: the lumped heat capacity.
- `oobq::V`: the out-of-balance heat load.
- `dT::V`: the temperature's rate of change.
- `mcT::V`: the heat content, `m·c·T`.
- `T::V`: the temperature [K].

`T1`, `T2` and `D` are the index type, float type and dimension; `P`, `B` and `V` the fields' types.

# Constructors
    MeshThermalPhase{T1,T2,D}(prprt, bcs, cᵢ, oobq, dT, mcT, T)
"""
struct MeshThermalPhase{T1,T2,D,P,B,V}
    prprt ::P
    bcs   ::B
    cᵢ    ::V # consistent lumped heat capacity matrix
    oobq  ::V # out-of-balance heat load
    dT    ::V # temperature rate of change
    mcT   ::V # heat capacity
    T     ::V # temperature
end
function MeshThermalPhase{T1,T2,D}(prprt::P, bcs::B, cᵢ::V, oobq::V, dT::V, mcT::V, T::V) where {T1,T2,D,P,B,V}
    return MeshThermalPhase{T1,T2,D,P,B,V}(prprt, bcs, cᵢ, oobq, dT, mcT, T)
end
function Adapt.adapt_structure(to, t::MeshThermalPhase{T1,T2,D}) where {T1,T2,D}
    return adapt_fields(to, t, MeshThermalPhase{T1,T2,D})
end

"""
    Mesh{T1,T2,D,P,VV,V,S,TH}

The background grid (Eulerian): its nodes and their fields. It holds no connectivity, which is the `Basis`'s.

# Fields
- `prprt::P`: the grid's sizes, the `MeshProperties`.
- `x₀::SVector{D,T2}`: the lower corner [m].
- `x::VV`: the node coordinates, an `SVector{D,T2}` per node (x fastest, then y, then z).
- `ΔJ::V`: the mass-weighted volume ratio of the step per node, for F-bar (`stab.locking`).
- `s::S`: the solid phase, a `MeshSolidPhase`.
- `t::TH`: the thermal phase, a `MeshThermalPhase` for a thermal case, else `nothing`.

`T1`, `T2` and `D` are the index type, float type and dimension; `P`, `S` and `TH` the fields' types; `VV` and `V` the
arrays' (of `SVector{D,T2}`, of `T2`; host or device).

# Constructors
    Mesh{T1,T2,D}(prprt, x₀, x, ΔJ, s, t)

`setup_mesh` builds the grid of a case's domain.

See also: `Point`, `Basis`, `Problem`.
"""
struct Mesh{T1,T2,D,P,VV,V,S,TH}
    prprt ::P
    # nodal quantities
    x₀    ::SVector{D,T2}
    x     ::VV # node coordinates : SVector{ndim,T2}
    ΔJ    ::V
    # solid phase
    s     ::S
    # thermal phase (nothing without one)
    t     ::TH
end
function Mesh{T1,T2,D}(prprt::P, x₀::SVector{D,T2}, x::VV, ΔJ::V, s::S, t::TH) where {T1,T2,D,P,VV,V,S,TH}
    return Mesh{T1,T2,D,P,VV,V,S,TH}(prprt, x₀, x, ΔJ, s, t)
end
function Adapt.adapt_structure(to, mesh::Mesh{T1,T2,D}) where {T1,T2,D}
    return adapt_fields(to, mesh, Mesh{T1,T2,D})
end
