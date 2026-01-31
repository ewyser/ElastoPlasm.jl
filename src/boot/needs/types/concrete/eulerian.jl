# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mesh Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Mesh

struct MeshProperties{T1,T2,D}
    # general information
    dim  ::T1
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

struct MeshSolidPhase{T1,T2,D} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2,D}
    bcs   ::MeshBoundary
    mᵢ    ::Vector{T2} # consistent lumped mass matrix
    Mᵢⱼ   ::Matrix{T2}
    oobf  ::Matrix{T2} # out-of-balance mechanical load
    a     ::Matrix{T2} # acceleration
    mv    ::Matrix{T2} # momentum
    v     ::Matrix{T2} # velocity
end
@adapt_struct MeshSolidPhase

struct MeshThermalPhase{T1,T2,D} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2,D}
    bcs   ::MeshBoundary
    cᵢ    ::Vector{T2} # consistent lumped heat capacity matrix
    oobq  ::Vector{T2} # out-of-balance heat load
    dT    ::Vector{T2} # temperature rate of change
    mcT   ::Vector{T2} # heat capacity
    T     ::Vector{T2} # temperature
end
@adapt_struct MeshThermalPhase

struct Mesh{T1,T2,D} <: UniformMesh{T1, T2}
    prprt ::MeshProperties{T1,T2,D}
    # nodal quantities
    x₀    ::Vector{T2}
    x     ::Matrix{T2}
    ΔJ    ::Matrix{T2}
    # solid phase
    s     ::MeshSolidPhase{T1,T2,D} # phase ::Vector{MeshPhase{T1,T2}}
    # thermal phase
    t     ::MeshThermalPhase{T1,T2,D} # phase ::Vector{MeshPhase{T1,T2}}
    # connectivity
    e2n   ::Matrix{T1}
    e2e   ::SparseMatrixCSC{T1,T1}
end
@adapt_struct Mesh
