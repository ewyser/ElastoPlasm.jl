export Self

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
    functional::Vector{String} = ["Functional execution platform(s):"]
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

# ============================================================================================================= # 
# Concrete Parametric AbstractArray Types
export Mesh,Point,Liquid,Solid

abstract type AbstractEulerian end
abstract type UniformCartesian{T1, T2} <: AbstractEulerian end
abstract type NonUniformCartesian{T1, T2} <: AbstractEulerian end


struct Boundary{B}
    status::Matrix{B}
end
@adapt_struct Boundary
struct Topology{T1,T2}
    e2n  ::Matrix{T1}
    e2e  ::SparseMatrixCSC{T1,T1}
    xB   ::Vector{T2}
end
@adapt_struct Topology
struct Field{T1,T2}
    x₀   ::Vector{T2}
    x    ::Matrix{T2}
    mᵢ   ::Vector{T2}
    Mᵢⱼ  ::Matrix{T2}
    oobf ::Matrix{T2}
    D    ::Matrix{T2}
    f    ::Matrix{T2}
    a    ::Matrix{T2}
    p    ::Matrix{T2}
    v    ::Matrix{T2}
    Δu   ::Matrix{T2}
    ΔJ   ::Matrix{T2}
    bij  ::Array{T2,3}
end
@adapt_struct Field


struct Mesh{T1,T2,B,NT}
    dim  ::T1
    nel  ::Vector{T1}
    nno  ::Vector{T1}
    nn   ::T1
    L    ::Vector{T2}
    h    ::Vector{T2}
    # nodal quantities
    x₀   ::Vector{T2}
    x    ::Matrix{T2}
    mᵢ   ::Vector{T2}
    Mᵢⱼ  ::Matrix{T2}
    oobf ::Matrix{T2}
    D    ::Matrix{T2}
    f    ::Matrix{T2}
    a    ::Matrix{T2}
    p    ::Matrix{T2}
    v    ::Matrix{T2}
    Δu   ::Matrix{T2}
    ΔJ   ::Matrix{T2}
    bij  ::Array{T2,3}
    # mesh-to-node topology
    e2n  ::Matrix{T1}
    e2e  ::SparseMatrixCSC{T1,T1}
    xB   ::Matrix{T2}
    # mesh boundary conditions
    bcs  ::Boundary{B}
end
@adapt_struct Mesh







abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end


struct Solid{T1,T2}
    u    ::Matrix{T2}
    v    ::Matrix{T2}
    p    ::Matrix{T2}
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    m    ::Vector{T2}
    c₀   ::Vector{T2}
    cᵣ   ::Vector{T2}
    ϕ    ::Vector{T2}
    Δλ   ::Vector{T2}
    ϵpII ::Matrix{T2}
    ϵpV  ::Vector{T2}
    ΔJ   ::Vector{T2}
    J    ::Vector{T2}
    # tensor in voigt notation
    σᵢ   ::Matrix{T2}
    τᵢ   ::Matrix{T2}
    # tensor in matrix notation
    δᵢⱼ  ::Matrix{T2}
    ∇vᵢⱼ ::Array{T2,3}
    ∇uᵢⱼ ::Array{T2,3}
    ΔFᵢⱼ ::Array{T2,3}
    Fᵢⱼ  ::Array{T2,3}
    Bᵢⱼ  ::Array{T2,3}
    ϵᵢⱼ  ::Array{T2,3}
    ωᵢⱼ  ::Array{T2,3}
    σJᵢⱼ ::Array{T2,3}
end
@adapt_struct Solid

struct Liquid{T1,T2}
    # Add concrete fields as needed, e.g.:
    # p    ::Vector{T2}
    # v    ::Matrix{T2}
end
@adapt_struct Liquid

struct Point{T1,T2}
    # general information
    ndim ::T1
    nmp  ::T1
    vmax ::Vector{T2}
    # basis-related quantities
    ϕ∂ϕ  ::Array{T2,3}
    δnp  ::Array{T2,3}
    # connectivity
    e2p  ::Matrix{T1}
    p2p  ::Matrix{T1}
    p2e  ::Vector{T1}
    p2n  ::Matrix{T1}
    # material point properties
    x    ::Matrix{T2}
    ℓ₀   ::Matrix{T2}
    ℓ    ::Matrix{T2}
    Ω₀   ::Vector{T2}
    Ω    ::Vector{T2}
    # solid phase
    s    ::Solid{T1,T2}
    # liquid phase
    l    ::Liquid{T1,T2}
end
@adapt_struct Point