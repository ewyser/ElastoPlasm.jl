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




abstract type AbstractStress{T} <: AbstractTensor{T} end

struct CauchyStress{S, T, L} <: AbstractStress{T}
    p  ::T
    dev::SMatrix{S,S,T,L}
end
struct KirchhoffStress{S, T, L} <: AbstractStress{T}
    p  ::T
    dev::SMatrix{S,S,T,L}
end



@inline function _get_tensor(stress::AbstractStress{S1, S2, T, L}) where {S1, S2, T, L}
    return stress.dev - stress.p * SMatrix{L,L,T}(I) # positive in compression
end


