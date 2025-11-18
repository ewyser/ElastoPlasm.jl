export Self

Base.@kwdef mutable struct Path
	root::String = " "
	out ::String = " "
	test::String = " "
	lib ::Dict   = Dict()
    ncell::Int   = 0
end
Base.@kwdef mutable struct Distributed
    is_active::Bool       = false
    root     ::NamedTuple = (;name = "root", val = 0,)
    glob     ::String     = "glob.jld2"
end
Base.@kwdef mutable struct Execution
    functional::Vector{String} = ["Available execution platform(s):"]
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
# Instruction, Time Types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export Instruction, Time

struct Instruction{T1<:Integer,T2<:Real,D<:AbstractDimension}
    dtype  ::NamedTuple
    basis  ::NamedTuple
    fwrk   ::NamedTuple
    bcs    ::NamedTuple
    grf    ::NamedTuple
    plast  ::NamedTuple
    nonloc ::NamedTuple
    plot   ::NamedTuple
    perf   ::NamedTuple
    backend::NamedTuple
    cairn  ::NamedTuple
end
@adapt_struct Instruction



struct Time{T1<:Integer,T2<:Real}
    t  ::Vector{T2}
    te ::T2
    tg ::T2
    tep::T2
end
@adapt_struct Time
