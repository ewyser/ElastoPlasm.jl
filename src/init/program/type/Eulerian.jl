export AbstractEulerian
export UniformCartesian
export Mesh

abstract type AbstractEulerian end
abstract type UniformCartesian{T1, T2} <: AbstractEulerian end
abstract type NonUniformCartesian{T1, T2} <: AbstractEulerian end

struct Mesh{T1,T2,
    T3 <: AbstractVector{T1},
    T5 <: AbstractMatrix{T1},
    T7 <: AbstractArray{T1}, 
    T4 <: AbstractVector{T2}, 
    T6 <: AbstractMatrix{T2}, 
    T8 <: AbstractArray{T2}, 
} <: UniformCartesian{T1,T2}
    ndim   ::T1
    nel  ::T3
    nno  ::T3
    nn   ::T1
    L    ::T4
    h    ::T4
    minC ::T4 
    # nodal quantities
    xn   ::T6
    mn   ::T4
    Mn   ::T4
    oobf ::T6
    Dn   ::T6
    fn   ::T6
    an   ::T6
    pn   ::T6
    vn   ::T6
    Δun  ::T6
    ΔJn  ::T6
    bn   ::T8
    # mesh-to-node topology
    e2n  ::T5
    e2e  ::T5
    xB   ::T4
    # mesh boundary conditions
    bc   ::T6
end