# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Abstract and concrete Dimension types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export get_dimension, OneDimension, TwoDimension, ThreeDimension

struct OneDimension <: AbstractDimension end

struct TwoDimension <: AbstractDimension end

struct ThreeDimension <: AbstractDimension end