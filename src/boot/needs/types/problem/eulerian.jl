# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mesh Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

abstract type AbstractEulerian end
abstract type AbstractCartesianMesh{T1, T2}  <: AbstractEulerian end
abstract type AbstractMeshPhase{T1, T2}      <: AbstractCartesianMesh{T1,T2} end

export Mesh

# The mesh, point and basis containers hold their arrays through type parameters (`V`, `VV`, `M`, …, appended after
# the parameters code dispatches on), so the same struct holds host `Array`s or a device's arrays (`adapt`, see
# `adapt_fields` in utils.jl); `Mesh{T1,T2,D}` etc. still name every version. They are left unbounded: bounds such
# as `V<:AbstractVector{T2}` refer to the other parameters, which made every subtype check on these types, hence
# compilation, about three times slower (a first `Problem(case)` 10 → 34 s); the partial constructors (`Mesh{T1,T2,D}(…)`)
# and the setup functions give them their types.

# The grid's sizes: static vectors (`N = D+1`: per axis, then the total; `M = 2D`), so it is a plain value, read as it
# is on the host (launch sizes, plots) and inside kernels; it needs no adapt rule.
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

struct MeshBoundary{B<:AbstractMatrix{Bool}}
    status::B
end
@adapt_struct MeshBoundary

struct MeshSolidPhase{T1,T2,D,P,B,V,M,VV} <: AbstractMeshPhase{T1,T2}
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

struct MeshThermalPhase{T1,T2,D,P,B,V} <: AbstractMeshPhase{T1,T2}
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

struct Mesh{T1,T2,D,P,VV,V,S,TH} <: AbstractCartesianMesh{T1, T2}
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
