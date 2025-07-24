export Self
export AbstractEulerian,UniformCartesian,Mesh
export AbstractLagrangian,MaterialPoint,Point

Base.@kwdef mutable struct Path
	root::String = ROOT
	out ::String = joinpath(dirname(ROOT),"dump")
	test::String = joinpath(dirname(ROOT),"test")
	lib ::Dict   = Dict()
    ncell::Int   = 0
end

Base.@kwdef mutable struct Distributed
    is_active::Bool       = false
    root     ::NamedTuple = (;name = "root", val = 0,)
    glob     ::String     = "glob.jld2"
end

Base.@kwdef mutable struct Execution
    functional::Vector{String} = ["Functional execution plateform(s):"]
	cpu::Dict = Dict()
	gpu::Dict = Dict()
end

Base.@kwdef mutable struct UI
	ui    ::Bool       = true
    plot  ::Bool       = false
    bckd  ::Bool       = true 
    logs  ::NamedTuple = (; log_info = true, log_warn = true, log_error = true) 
end

Base.@kwdef mutable struct Self
	sys ::Path
	ui  ::UI
	bckd::Execution
    mpi ::Distributed
end

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



abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end

struct Point{T1,T2,
    T3 <: AbstractVector{T1},
    T5 <: AbstractMatrix{T1},
    T7 <: AbstractArray{T1} , 
    T4 <: AbstractVector{T2}, 
    T6 <: AbstractMatrix{T2}, 
    T8 <: AbstractArray{T2} ,} <: MaterialPoint{T1, T2}
    # definition
    nmp  ::T1
    ℓ₀   ::T6
    ℓ    ::T6
    Ω₀   ::T4
    Ω    ::T4
    m    ::T4
    c₀   ::T4
    cᵣ   ::T4
    ϕ    ::T4
    Δλ   ::T4
    ϵpII ::T4
    ϵpV  ::T4
    ΔJ   ::T4
    J    ::T4
    x    ::T6
    u    ::T6
    v    ::T6
    p    ::T6
    # plot quantity
    z₀   ::T4
    # tensor in voigt notation
    σᵢ   ::T6
    τᵢ   ::T6
    # additional quantities
    I    ::T6
    ϕ∂ϕ  ::T8
    δnp  ::T8
    # tensor in matrix notation
    ∇vᵢⱼ ::T8
    ∇uᵢⱼ ::T8  
    ΔFᵢⱼ ::T8
    Fᵢⱼ  ::T8 
    bᵢⱼ  ::T8 
    ϵᵢⱼ  ::T8 
    ωᵢⱼ  ::T8 
    σJᵢⱼ ::T8 
    # topology & connectivity
    p2e  ::T4
    e2p  ::T6
    p2p  ::T6
    p2n  ::T6
end