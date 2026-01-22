# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Abstract and concrete Dimension types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export Dimension, TwoDimension, ThreeDimension

struct TwoDimension{T1<:Integer} <: AbstractDimension 
    ndim::T1
    nstr::T1
end

struct ThreeDimension{T1<:Integer} <: AbstractDimension 
    ndim::T1
    nstr::T1
end

function Dimension(ndim::T1) where {T1<:Integer}
    if ndim == T1(2)
        return TwoDimension{T1}(2, 3)
    elseif ndim == T1(3)
        return ThreeDimension{T1}(3, 6)
    else
        error("Unsupported dimension type")
    end
end
#= definition example
mpts = Point{T1,T2,typeof(elast),typeof(rheo),TwoDimension}(
    ...all required fields..., 
)
=#