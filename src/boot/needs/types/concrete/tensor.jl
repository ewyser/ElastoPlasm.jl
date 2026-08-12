abstract type AbstractTensor{T} end

abstract type AbstractStrain{S,T,L} <: AbstractTensor{T} end

@inline function get_tensor(strain::AbstractStrain{S, T, L}) where {S, T, L}
    return strain.dev + strain.vol * SMatrix{S,S,T,L}(I)
end

function Base.getindex(strain::AbstractStrain{S, T, L}, i::Int, j::Int) where {S, T, L}
    strain.dev[i,j] + (i == j ? strain.vol : zero(T))
end

struct InfinitesimalStrain{S, T, L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
end

struct LogarithmicStrain{S, T, L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
end

function LinearAlgebra.eigen(strain::LogarithmicStrain{S, T, L}) where {S, T, L}
    return eigen(Symmetric(_get_tensor(strain)))
end




abstract type AbstractStress{S,T,L} <: AbstractTensor{T} end

@inline function _get_tensor(stress::AbstractStress{S, T, L}) where {S, T, L}
    return stress.dev - stress.p * SMatrix{S,S,T,L}(I) # positive in compression
end

struct CauchyStress{S, T, L} <: AbstractStress{S, T, L}
    p  ::T
    dev::SMatrix{S,S,T,L}
end
struct KirchhoffStress{S, T, L} <: AbstractStress{S, T, L}
    p  ::T
    dev::SMatrix{S,S,T,L}
end






