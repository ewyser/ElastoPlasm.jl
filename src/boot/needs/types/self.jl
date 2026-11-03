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
"""
    Self

The package's state in a session, `ElastoPlasm.self`, set when it loads.

# Fields
- `sys::Path`: where things are: `root` (the package), `dump` (simulations' default parent folder), `test`, and `lib`
  (the source tree `superInc` loaded, for the documentation's pages); `ncell` is not read.
- `ui::UI`: the session: `plot` (plots can be drawn; `solver_instr` turns plotting off otherwise); `ui`, `bckd` and
  `logs` are not read.
- `mpi::Distributed`: multi-process settings, not read.
"""
Base.@kwdef mutable struct Self
	sys ::Path
	ui  ::UI
    mpi ::Distributed
end