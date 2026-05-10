# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Point,PointSolidPhase,PointFluidPhase,PointThermalPhase

struct LinearElasticity{T1,T2,D} <: AbstractElasticity{T1,T2}
    # tensor in voigt notation
    σᵢ   ::Matrix{T2}
    # tensor in matrix notation
    ϵᵢⱼ  ::Array{T2,3}
    ωᵢⱼ  ::Array{T2,3}
    σJᵢⱼ ::Array{T2,3}
end
function LinearElasticity(::Type{T1},::Type{T2},nmp::Integer,ndim::Integer) where {T1,T2}
    if ndim == 2 
        nstr = 3 
    elseif ndim == 3 
        nstr = 6 
    end
    σᵢ   = zeros(T2, nstr, nmp)
    ϵᵢⱼ  = zeros(T2, ndim, ndim, nmp)
    ωᵢⱼ  = zeros(T2, ndim, ndim, nmp)
    σJᵢⱼ = zeros(T2, ndim, ndim, nmp)
    return LinearElasticity{T1,T2,ndim}(σᵢ, ϵᵢⱼ, ωᵢⱼ, σJᵢⱼ)
end
@adapt_struct LinearElasticity

struct FiniteElasticity{T1,T2,D} <: AbstractElasticity{T1,T2}
    # tensor in voigt notation
    σᵢ   ::Matrix{T2}
    τᵢ   ::Matrix{T2}
    # tensor in matrix notation
    bᵢⱼ  ::Array{T2,3}
    ϵᵢⱼ  ::Array{T2,3}
end
function FiniteElasticity(::Type{T1}, ::Type{T2}, nmp::Integer, ndim::Integer) where {T1,T2}
    if ndim == 2 
        nstr = 3 
    elseif ndim == 3 
        nstr = 6 
    end
    σᵢ  = zeros(T2, nstr, nmp)
    τᵢ  = zeros(T2, nstr, nmp)
    bᵢⱼ = zeros(T2, ndim, ndim, nmp)
    ϵᵢⱼ = zeros(T2, ndim, ndim, nmp)
    return FiniteElasticity{T1,T2,ndim}(σᵢ, τᵢ, bᵢⱼ, ϵᵢⱼ)
end
@adapt_struct FiniteElasticity

struct DruckerPragerRheology{T1,T2,D} <: AbstractRheology{T1,T2}
    c₀   ::Vector{T2}
    cᵣ   ::Vector{T2}
    ϕ    ::Vector{T2}
    Δλ   ::Vector{T2}
    ϵpII ::Matrix{T2}
    ϵpV  ::Vector{T2}
end
@adapt_struct DruckerPragerRheology

struct PointSolidPhase{T1,T2,D,E<:AbstractElasticity,R<:AbstractRheology} <: MaterialPointPhase{T1,T2}
    u    ::Matrix{T2}
    v    ::Matrix{T2}
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    c₀   ::Vector{T2}
    cᵣ   ::Vector{T2}
    ϕ₀   ::Vector{T2}
    Δλ   ::Vector{T2}
    ϵpII ::Matrix{T2}
    ϵpV  ::Vector{T2}
    # tensor in voigt notation
    σᵢ   ::Matrix{T2}
    τᵢ   ::Matrix{T2}
    # tensor in matrix notation
    ∇vᵢⱼ ::Array{T2,3}
    ∇uᵢⱼ ::Array{T2,3}
    ΔFᵢⱼ ::Array{T2,3}
    Fᵢⱼ  ::Array{T2,3}
    bᵢⱼ  ::Array{T2,3}
    ϵᵢⱼ  ::Array{T2,3}
    ωᵢⱼ  ::Array{T2,3}
    σJᵢⱼ ::Array{T2,3}
    # new component-based fields
    elast::E
    rheo ::R
end
@adapt_struct PointSolidPhase

struct PointFluidPhase{T1,T2,D} <: MaterialPointPhase{T1,T2}
    # Add concrete fields as needed, e.g.:
    # v    ::Matrix{T2}
end
@adapt_struct PointFluidPhase

struct PointThermalPhase{T1,T2,D} <: MaterialPointPhase{T1,T2}
    c   ::Vector{T2} # specific heat capacity vector
    k   ::Vector{T2} # thermal conductivity vector
    q   ::Matrix{T2} # heat flux array
    T   ::Vector{T2} # temperature vector
end
@adapt_struct PointThermalPhase

struct Point{T1,T2,D<:AbstractDimension,B<:AbstractBasis,E<:AbstractElasticity,R<:AbstractRheology} <: MaterialPoint{T1,T2}
    # basis
    basis::B
    # general information
    ndim ::T1
    nmp  ::T1
    # CFL-related quantity
    vmax ::Vector{T2}
    # basis-related quantities
    Δnp  ::Array{T2,3}
    # APIC-related
    Bᵢⱼ  ::Array{T2,3}
    Dᵢⱼ  ::Array{T2,3}
    # connectivity
    nn   ::T1
    e2p  ::Matrix{T1}
    p2p  ::Matrix{T1}
    p2e  ::Vector{T1}
    p2n  ::Matrix{T1}
    # utils
    δᵢⱼ  ::Matrix{T2}
    # material point properties
    x    ::Matrix{T2}
    ℓ₀   ::Matrix{T2}
    ℓ    ::Matrix{T2}
    n₀   ::Vector{T2}
    n    ::Vector{T2}    
    Ω₀   ::Vector{T2}
    Ω    ::Vector{T2}
    ΔJ   ::Vector{T2}
    J    ::Vector{T2}
    # solid phase
    s    ::PointSolidPhase{T1,T2,D,E,R}
    # fluid phase
    f    ::Union{Nothing, PointFluidPhase{T1,T2,D}}
    # thermal phase
    t    ::Union{Nothing, PointThermalPhase{T1,T2,D}}
end
@adapt_struct Point
