# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Abstract and concrete Dimension types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export get_dimension, OneDimension, TwoDimension, ThreeDimension

struct OneDimension <: AbstractDimension end

struct TwoDimension <: AbstractDimension end

struct ThreeDimension <: AbstractDimension end

#= definition example
mpts = Point{T1,T2,typeof(elast),typeof(rheo),TwoDimension}(
    ...all required fields..., 
)
=#