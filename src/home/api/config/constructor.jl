export kwargser

"""
    kwargser(kwargs::Any; dim::Number=2, constructor::Type=Instruction) -> Instruction

Generate a configuration object by merging user-supplied keyword arguments with default values.

# Arguments
- `kwargs::Any`: Dictionary or named tuple of user-supplied keyword arguments to override defaults
- `::Type{Instruction}`: The type to construct
- `dim::Number=2`: Spatial dimension, used for kernel initialization

# Returns
- `Instruction`: Fully configured instruction object with kernels and precision settings

# Example
```julia
# Top-level overrides
instr = kwargser(Dict(:dtype => 64, :locking => true); dim=3)
println(instr.cairn.shpfun)  # prints the initialized shape function kernel

# To override nested fields, use a nested NamedTuple or Dict:
instr = kwargser((;fwrk = (deform="infinitesimal", trsfr="tpic", locking=false)))
println(instr.fwrk.deform)  # prints "infinitesimal"
```

# Notes
- Merges `kwargs` with default configuration
- Sets precision and attaches kernel functions for shape functions, mapping, and constitutive updates
- Warns about any unused keyword arguments
"""
function kwargser(kwargs::Any; dim::Number=2, constructor::Type=Instruction)
    ref    = get_default(constructor)
    # Only keep keys in kwargs that are also in ref
    valids = intersect(keys(kwargs), keys(ref))
    user   = NamedTuple{Tuple(valids)}(getindex.(Ref(kwargs), valids))
    instr  = merge(ref, user)
    # Warn about unused keys
    unused = setdiff(keys(kwargs), keys(ref))
    if !isempty(unused)
        @warn join(vcat("miscellaneous kwargs:", "\n\t- ".*String.(unused)))
    end
    # Set precision
    if instr.dtype == 64
        instr = merge(instr, (; dtype = (;T0=(Int64,Float64),bits=Int64(64),precision="64-bit precision (or double precision)"),))
    elseif instr.dtype == 32
        instr = merge(instr, (;dtype = (;T0=(Int32,Float32),bits=Int32(32),precision="32-bit precision (or single precision)"),))
    end
    # Get dimension type
    if dim == 1
        dim = OneDimension()
    elseif dim == 2
        dim = TwoDimension()
    elseif dim == 3
        dim = ThreeDimension()
    else
        error("Unsupported dimension type: $ndim")
    end
    # Verify plotting status from ui
    if !ElastoPlasm.info.ui.plot
        instr = merge(instr, (; plot = merge(instr.plot, (; status = false,)),))
    end    
    # Set execution backend
    if !haskey(instr[:backend],:exec)
        exec  = select_execution_backend(instr[:backend][:select]; mpi_status=instr[:backend][:distributed])
        instr = merge(instr, (; backend = merge(instr[:backend], (; exec = exec,)),))
    end
    # Add cairns to instr     
    cairn = (;
        shpfun = init_shpfun(instr),
        mapsto = init_mapsto(instr),
        update = init_update(instr),
    )
    instr = merge(instr, (; cairn = cairn,))
    # Create Instruction type
    return constructor{instr[:dtype][:T0]...,typeof(dim)}(
        instr[:dtype],
        instr[:basis],
        instr[:fwrk],
        instr[:bcs],
        instr[:grf],
        instr[:plast],
        instr[:nonloc],
        instr[:plot],
        instr[:perf],
        instr[:backend],
        instr[:cairn],
    )
end
