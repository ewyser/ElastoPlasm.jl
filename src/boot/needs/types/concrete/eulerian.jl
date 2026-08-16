# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mesh Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Mesh

struct MeshProperties{T1,T2,D,NN}
    # general information
    nel  ::Vector{T1}
    nno  ::Vector{T1}
    nn   ::T1
    L    ::Vector{T2}
    h    ::Vector{T2}
    xB   ::Matrix{T2}
end
@adapt_struct MeshProperties

struct MeshBoundary
    status::Matrix{Bool}
end
@adapt_struct MeshBoundary

struct MeshSolidPhase{T1,T2,D,NN} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2,D,NN}
    bcs   ::MeshBoundary
    m     ::Vector{T2} # consistent lumped mass matrix
    Mᵢⱼ   ::Matrix{T2}
    oobf  ::Matrix{T2} # out-of-balance mechanical load (atomic writes: keep Matrix)
    oobp  ::Vector{T2} # out-of-balance continuity residual (u-P)
    a     ::Vector{SVector{D,T2}} # acceleration per node : SVector{ndim,T2}
    mv    ::Matrix{T2} # momentum (atomic writes: keep Matrix)
    v     ::Vector{SVector{D,T2}} # velocity per node    : SVector{ndim,T2}
    p     ::Vector{T2} # nodal pressure (u-P)
end
@adapt_struct MeshSolidPhase

struct MeshThermalPhase{T1,T2,D,NN} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2,D,NN}
    bcs   ::MeshBoundary
    cᵢ    ::Vector{T2} # consistent lumped heat capacity matrix
    oobq  ::Vector{T2} # out-of-balance heat load
    dT    ::Vector{T2} # temperature rate of change
    mcT   ::Vector{T2} # heat capacity
    T     ::Vector{T2} # temperature
end
@adapt_struct MeshThermalPhase

struct Mesh{T1,T2,D,NN} <: UniformMesh{T1, T2}
    prprt ::MeshProperties{T1,T2,D,NN}
    # nodal quantities
    x₀    ::Vector{T2}
    x     ::Vector{SVector{D,T2}} # node coordinates : SVector{ndim,T2}
    type  ::Matrix{T1}
    ΔJ    ::Vector{T2}
    # solid phase
    s     ::MeshSolidPhase{T1,T2,D,NN}
    # thermal phase
    t     ::Union{Nothing, MeshThermalPhase{T1,T2,D,NN}}
    # connectivity
    e2n   ::Vector{SVector{NN,T1}}
    e2e   ::SparseMatrixCSC{T1,T1}
end
@adapt_struct Mesh
