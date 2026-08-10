struct InfinitesimalStrain{T,D} <: AbstractStrain{T,D}
    vol::T
    dev::SMatrix{D,D,T}
end
struct LogarithmicStrain{T,D} <: AbstractStrain{T,D}
    vol::T
    dev::SMatrix{D,D,T}
end


struct CauchyStress{T,D} <: AbstractStress{T,D}
    p  ::T
    dev::SMatrix{D,D,T}
end
struct KirchhoffStress{T,D} <: AbstractStress{T,D}
    p  ::T
    dev::SMatrix{D,D,T}
end


@inline function _get_tensor(strain::AbstractStrain{T,D}) where {T,D}
    return strain.dev + strain.vol * SMatrix{D,D,T}(I)
end
@inline function _get_tensor(stress::AbstractStress{T,D}) where {T,D}
    return stress.dev - stress.p * SMatrix{D,D,T}(I) # positive in compression
end


