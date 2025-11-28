# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dimension types for Mesh and MaterialPoint types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export TwoDimension, ThreeDimension
abstract type AbstractDimension end
struct TwoDimension{T1<:Integer} <: AbstractDimension 
    ndim::T1
    nstr::T1
end
function TwoDimension(::Type{T1}) where {T1<:Integer}
    ndim = T1(2) 
    nstr = T1(3) 
    return TwoDimension{T1}(ndim,nstr)
end
@adapt_struct TwoDimension

struct ThreeDimension{T1<:Integer} <: AbstractDimension 
    ndim::T1
    nstr::T1
end
function ThreeDimension(::Type{T1}) where {T1<:Integer}
    ndim = T1(3) 
    nstr = T1(6) 
    return ThreeDimension{T1}(ndim,nstr)
end
@adapt_struct ThreeDimension
#= definition example
mpts = Point{T1,T2,typeof(elast),typeof(rheo),TwoDimension}(
    ...all required fields..., 
)
=#