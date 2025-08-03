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
struct Mesh{T1,T2,
    A3 <: AbstractArray{T1,1},
    A5 <: AbstractArray{T1,2},
    A7 <: AbstractArray{T1,3},
    A4 <: AbstractArray{T2,1}, 
    A6 <: AbstractArray{T2,2}, 
    A8 <: AbstractArray{T2,3}} <: UniformCartesian{T1,T2}
    dim  ::T1
    nel  ::A3
    nno  ::A3
    nn   ::T1
    L    ::A4
    h    ::A4
    # nodal quantities
    x₀   ::A4 
    x    ::A6
    mᵢ   ::A4
    Mᵢⱼ  ::A6
    oobf ::A6
    D    ::A6
    f    ::A6
    a    ::A6
    p    ::A6
    v    ::A6
    Δu   ::A6
    ΔJ   ::A6
    bij  ::A8
    # mesh-to-node topology
    e2n  ::A5
    e2e  ::SparseMatrixCSC{T1,T1}
    xB   ::A4
    # mesh boundary conditions
    bc   ::A6
end

abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end


struct Solid{T1,T2,
    A3 <: AbstractArray{T1,1},
    A5 <: AbstractArray{T1,2},
    A7 <: AbstractArray{T1,3},
    A4 <: AbstractArray{T2,1}, 
    A6 <: AbstractArray{T2,2}, 
    A8 <: AbstractArray{T2,3}} <: MaterialPoint{T1, T2}
    u    ::A6 
    v    ::A6
    p    ::A6
    # mechanical properties
    m    ::A4
    c₀   ::A4
    cᵣ   ::A4
    ϕ    ::A4
    Δλ   ::A4
    ϵpII ::A6
    ϵpV  ::A4
    ΔJ   ::A4
    J    ::A4
    # tensor in voigt notation
    σᵢ   ::A6
    τᵢ   ::A6 
    # tensor in matrix notation
    δᵢⱼ  ::A6
    ∇vᵢⱼ ::A8
    ∇uᵢⱼ ::A8
    ΔFᵢⱼ ::A8
    Fᵢⱼ  ::A8
    Bᵢⱼ  ::A8
    ϵᵢⱼ  ::A8
    ωᵢⱼ  ::A8
    σJᵢⱼ ::A8
end
@adapt_struct Solid

struct Liquid{T1,T2,
    A3 <: AbstractArray{T1,1},
    A5 <: AbstractArray{T1,2},
    A7 <: AbstractArray{T1,3},
    A4 <: AbstractArray{T2,1}, 
    A6 <: AbstractArray{T2,2}, 
    A8 <: AbstractArray{T2,3}} <: MaterialPoint{T1, T2}
end
@adapt_struct Liquid

struct Point{T1,T2,
    A3 <: AbstractArray{T1,1},
    A5 <: AbstractArray{T1,2},
    A7 <: AbstractArray{T1,3},
    A4 <: AbstractArray{T2,1}, 
    A6 <: AbstractArray{T2,2}, 
    A8 <: AbstractArray{T2,3}} <: MaterialPoint{T1, T2}
    # general information
    ndim ::T1
    nmp  ::T1
    vmax ::A4
    # basis-related quantities
    ϕ∂ϕ  ::A8
    δnp  ::A8
    # connectivity
    e2p  ::A5
    p2p  ::A5
    p2e  ::A3
    p2n  ::A5
    # material point properties
    x    ::A6
    ℓ₀   ::A6
    ℓ    ::A6
    Ω₀   ::A4
    Ω    ::A4
    # solid phase
    s    ::Solid{T1,T2,A3,A5,A7,A4,A6,A8}
    # liquid phase
    l    ::Liquid{T1,T2,A3,A5,A7,A4,A6,A8}
end
@adapt_struct Point