export Self

Base.@kwdef mutable struct Path
	root::String = " "
	dump::String = " "
	test::String = " "
	lib ::Dict   = Dict()
    ncell::Int   = 0
end
Base.@kwdef mutable struct Distributed
    is_active::Bool       = false
    root     ::NamedTuple = (;name = "root", val = 0,)
    glob     ::String     = "glob.jld2"
end
Base.@kwdef mutable struct UI
	ui    ::Bool       = true
    plot  ::Bool       = true
    bckd  ::Bool       = true 
    logs  ::NamedTuple = (; log_info = true, log_warn = true, log_error = true) 
end
Base.@kwdef mutable struct Self
	sys ::Path
	ui  ::UI
    mpi ::Distributed
end