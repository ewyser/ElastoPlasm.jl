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

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Material Point Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Point,SolidMaterialPointPhase,FluidMaterialPointPhase,ThermalMaterialPointPhase

abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end
abstract type MaterialPointPhase{T1, T2} <: MaterialPoint{T1,T2} end

struct MaterialPointSolidPhase{T1,T2} <: MaterialPointPhase{T1,T2}
    u    ::Matrix{T2}
    v    ::Matrix{T2}
    # mechanical properties
    ρ₀   ::Vector{T2}
    ρ    ::Vector{T2}
    c₀   ::Vector{T2}
    cᵣ   ::Vector{T2}
    ϕ    ::Vector{T2}
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
end
@adapt_struct MaterialPointSolidPhase

struct MaterialPointFluidPhase{T1,T2} <: MaterialPointPhase{T1,T2}
    # Add concrete fields as needed, e.g.:
    # v    ::Matrix{T2}
end
@adapt_struct MaterialPointFluidPhase

struct MaterialPointThermalPhase{T1,T2} <: MaterialPointPhase{T1,T2}
    c   ::Vector{T2} # specific heat capacity vector
    k   ::Vector{T2} # thermal conductivity vector
    q   ::Matrix{T2} # heat flux array
    T   ::Vector{T2} # temperature vector
end
@adapt_struct MaterialPointThermalPhase

struct Point{T1,T2} <: MaterialPoint{T1,T2}
    # general information
    ndim ::T1
    nmp  ::T1
    # CFL-related quantity
    vmax ::Vector{T2}
    # basis-related quantities
    ϕ∂ϕ  ::Array{T2,3}
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
    s    ::MaterialPointSolidPhase{T1,T2}
    # fluid phase
    f    ::MaterialPointFluidPhase{T1,T2}
    # thermal phase
    t    ::MaterialPointThermalPhase{T1,T2}
end
@adapt_struct Point

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mesh Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Mesh

abstract type AbstractEulerian end
abstract type CartesianMesh{T1, T2}  <: AbstractEulerian end
abstract type UniformMesh{T1, T2}    <: CartesianMesh{T1, T2} end
abstract type NonUniformMesh{T1, T2} <: CartesianMesh{T1, T2} end
abstract type MeshPhase{T1, T2}      <: CartesianMesh{T1,T2} end

struct MeshProperties{T1,T2}
    # general information
    dim  ::T1
    nel  ::Vector{T1}
    nno  ::Vector{T1}
    nn   ::T1
    L    ::Vector{T2}
    h    ::Vector{T2}
    xB   ::Matrix{T2}
end
@adapt_struct MeshProperties

struct MeshBoundary{B}
    status::Matrix{B}
end
@adapt_struct MeshBoundary

struct MeshSolidPhase{T1,T2,B} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2}
    bcs   ::MeshBoundary{B}
    mᵢ    ::Vector{T2} # consistent lumped mass matrix
    Mᵢⱼ   ::Matrix{T2}
    oobf  ::Matrix{T2} # out-of-balance mechanical load
    a     ::Matrix{T2} # acceleration
    mv    ::Matrix{T2} # momentum
    v     ::Matrix{T2} # velocity
end
@adapt_struct MeshSolidPhase

struct MeshThermalPhase{T1,T2,B} <: MeshPhase{T1,T2}
    prprt ::MeshProperties{T1,T2}
    bcs   ::MeshBoundary{B}
    cᵢ    ::Vector{T2} # consistent lumped heat capacity matrix
    oobq  ::Vector{T2} # out-of-balance heat load
    Q     ::Vector{T2} # heat flux
    mcT   ::Vector{T2} # heat capacity
    T     ::Vector{T2} # temperature
end
@adapt_struct MeshThermalPhase

struct Mesh{T1,T2,B,NT} <: UniformMesh{T1, T2}
    prprt ::MeshProperties{T1,T2}
    # nodal quantities
    x₀    ::Vector{T2}
    x     ::Matrix{T2}
    ΔJ    ::Matrix{T2}
    # solid phase
    s     ::MeshSolidPhase{T1,T2,B} # phase ::Vector{MeshPhase{T1,T2}}
    # thermal phase
    t     ::MeshThermalPhase{T1,T2,B} # phase ::Vector{MeshPhase{T1,T2}}
    # connectivity
    e2n   ::Matrix{T1}
    e2e   ::SparseMatrixCSC{T1,T1}
end
@adapt_struct Mesh