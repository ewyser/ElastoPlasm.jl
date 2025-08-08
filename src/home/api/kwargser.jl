export kwargser
"""
    kwargser(type::Symbol, kwargs::Any; dim::Number=2) -> Dict

Generate a configuration dictionary for simulation or API routines by merging user-supplied keyword arguments with default reference values.

# Arguments
- `type::Symbol`: The configuration type (e.g., `:instr`) to load defaults for.
- `kwargs::Any`: Dictionary of user-supplied keyword arguments to override defaults.
- `dim::Number=2`: (Optional) Spatial dimension, used for kernel initialization (only relevant for `:instr`).

# Returns
- `Dict`: A dictionary containing the merged configuration. For `:instr`, this includes kernels and precision settings.

# Behavior
- For `type == :instr`, merges `kwargs` with required defaults, sets precision, and attaches kernel functions for shape functions, mapping, and constitutive updates.
- For other types, merges `kwargs` with defaults and returns the result.
- Warns about any unused keyword arguments.

# Example
```julia
cfg = kwargser(:instr, Dict(:dtype => 64, :locking => true); dim=3)
println(cfg[:cairn][:shpfun])  # prints the initialized shape function kernel
```
"""
function kwargser(type::Symbol,kwargs::Any;dim::Number = 2)
    ref  = require(type)
    key0 = collect(keys(kwargs))
    
    k,v  = [],[]
    if type == :instr
        for (i,key) ∈ enumerate(collect(keys(ref)))
            push!(k,key)
            if haskey(kwargs,key)
                push!(v,kwargs[key])
                filter!(x->x≠key,key0)
            else
                push!(v,ref[key])
            end
        end
        # display @warn message
        if !isempty(key0)
            @warn join(vcat("miscellaneous kwargs for require():","\n\t- ".*String.(key0)))
        end
        # zip & set precision
        instr = Dict(zip(k,v))
        if instr[:dtype] == 64
            instr[:dtype] = (;T0=(Int64,Float64),bits=Int64(64),precision="FP64")
        elseif instr[:dtype] == 32
            instr[:dtype] = (;T0=(Int32,Float32),bits=Int32(32),precision="FP32")
        end
        # add cairns (abstract kernels) to instr set
        instr[:cairn] = (;
            shpfun = init_shpfun(dim,instr),
            mapsto = init_mapsto(dim,instr),
            elastoplast = (;
                update = init_update(instr),
                elast  = init_elast(instr),
                plast  = init_plast(instr),
            )
        )
        return instr
    else
        for (i,var) ∈ enumerate(key)
            if haskey(kwargs,var)
                push!(values,kwargs[var])
                if var == :select
                    selection = kwargs[var]
                end
            else
                push!(values,ref[var])
                if var == :select
                    selection = ref[var]
                end
            end
        end
        return Dict(zip(key,value)) #return (;zip(Tuple(Symbol(x) for x in key),value)...)
    end
end