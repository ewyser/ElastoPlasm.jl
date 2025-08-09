export kwargser
"""
    kwargser(type::Symbol, kwargs::Any; dim::Number=2) -> NamedTuple

Generate a configuration named tuple for simulation or API routines by merging user-supplied keyword arguments with default reference values.

# Arguments
- `type::Symbol`: The configuration type (e.g., `:instr`) to load defaults for.
- `kwargs::Any`: Dictionary or named tuple of user-supplied keyword arguments to override defaults.
- `dim::Number=2`: (Optional) Spatial dimension, used for kernel initialization (only relevant for `:instr`).

# Returns
- `NamedTuple`: Merged configuration. For `:instr`, includes kernels and precision settings.

# Example
```julia
# Top-level overrides
cfg = kwargser(:instr, Dict(:dtype => 64, :locking => true); dim=3)
println(cfg.cairn.shpfun)  # prints the initialized shape function kernel

# To override nested fields, use a nested NamedTuple or Dict:
cfg = kwargser(:instr, (;fwrk = (deform="infinitesimal", trsfr="tpic", locking=false)))
println(cfg.fwrk.deform)  # prints "infinitesimal"
```

# Notes
- For `type == :instr`, merges `kwargs` with required defaults, sets precision, and attaches kernel functions for shape functions, mapping, and constitutive updates.
- For other types, merges `kwargs` with defaults and returns the result.
- Warns about any unused keyword arguments.
"""
function kwargser(type::Symbol, kwargs::Any; dim::Number=2)
    ref    = require(type)
    # Only keep keys in kwargs that are also in ref
    valids = intersect(keys(kwargs), keys(ref))
    user   = NamedTuple{Tuple(valids)}(getindex.(Ref(kwargs), valids))
    instr  = merge(ref, user)
    # Warn about unused keys
    unused = setdiff(keys(kwargs), keys(ref))
    if !isempty(unused)
        @warn join(vcat("miscellaneous kwargs for require():", "\n\t- ".*String.(unused)))
    end
    # Set precision
    if instr.dtype == 64
        instr = merge(instr, (dtype = (;T0=(Int64,Float64),bits=Int64(64),precision="FP64 precision"),))
    elseif instr.dtype == 32
        instr = merge(instr, (dtype = (;T0=(Int32,Float32),bits=Int32(32),precision="FP32 precision"),))
    end
    # Add cairns (abstract kernels) to instr set
    if type == :instr
        cairn = (
            shpfun = init_shpfun(dim, instr),
            mapsto = init_mapsto(dim, instr),
            elastoplast = (
                update = init_update(instr),
                elast  = init_elast(instr),
                plast  = init_plast(instr),
            )
        )
        instr = merge(instr, (cairn = cairn,))
    end
    return instr
end