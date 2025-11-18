export get_default, get_option, cli, kwargser
"""
    get_default(::Type{Instruction}) -> NamedTuple

Return default configuration values for the specified type.

# Arguments
- `::Type{Instruction}`: The type to get defaults for

# Returns
- `NamedTuple`: Default configuration values for simulation, including precision, basis, deformation framework, mapping scheme, locking mitigation, random field generator, plasticity, non-local regularization, plotting options, and performance mode.

# Example
```julia
cfg = get_default(Instruction)
println(cfg.basis.which)  # prints the default basis type
```

# Configuration Keys
- `:dtype`   — Arithmetic precision (e.g., 64 for Float64)
- `:basis`   — Shape function type and options
- `:fwrk`    — Deformation framework and transfer scheme
- `:bcs`     — Boundary condition settings
- `:grf`     — Gaussian Random Field generator options
- `:plast`   — Plasticity onset and flow law
- `:nonloc`  — Non-local regularization options
- `:plot`    — Plotting options
- `:perf`    — Performance mode options
"""
function get_default(::Type{Instruction})
    return (;
        dtype = 64,
        basis = (;
            which = "bsmpm",
            how = nothing,
            ghost = false,
        ),
        fwrk  = (;
            deform = "finite",
            trsfr = "std",
            C_pf = 1.0, 
            musl = true,
            locking = true,
            damping = 0.1
        ),
        bcs   = (;
            dirichlet = [
                :roller :roller;
                :roller :roller;
                :roller :roller], # for 2d, this translates to [lower_x upper_x;lower_y upper_y]
        ),
        grf   = (;
            status = false,
            covariance = "gaussian",
            param = (; 
                Iₓ= [2.5,2.5,2.5], 
                Nₕ = 5000, 
                kₘ = 100,
            ),
        ),
        plast = (;
            status = false,
            constitutive = "DP",
        ),
        nonloc = (;
            status=true,
            ls=0.5,
        ),
        plot  = (;
            status = true,
            freq   = 1.0,
            dpi    = 500,
            what   = [(;mpts=(name="epII",cblim=(0.0,1.5)),),],
        ),
        perf  = (;
            status=false,
        ),
        backend = (;
            select="host",
            distributed=false
        ),
    )
end

"""
    get_option(::Type{Instruction}) -> NamedTuple

Return all possible values for each configuration field as vectors.

# Arguments
- `::Type{Instruction}`: The type to get options for

# Returns
- `NamedTuple`: Structure matching `get_default(Instruction)` but with vectors of possible values instead of single values.

# Example
```julia
opts = get_option(Instruction)
println(opts.basis.which)  # ["bsmpm", "gimpm", "smpm"]
println(opts.fwrk.deform)  # ["finite", "infinitesimal"]
```

# Notes
- Useful for validation, UI generation, and documentation
- Fields with continuous values (like damping) show reasonable ranges
- Boolean fields are represented as `[true, false]`
"""
function get_option(::Type{Instruction})
    return (
            dtype = ("Select arithmetic precision",[32, 64]),
            basis = (
                which = ("Select basis type",["bsmpm", "gimpm", "smpm"]),
                how = ("Select material point domain update",[nothing]),
                ghost = ("Add ghost nodes ?",[true, false]),
            ),
            fwrk  = (
                deform = ("Select the deformation framework",["finite", "infinitesimal"]),
                trsfr = ("Select the mapping scheme",["std", "tpic", "apic"]),
                C_pf = ("Select picflip ratio",[1.0, 0.99, 0.95]),
                musl = ("Enable musl update",[true, false]),
                locking = ("Enable volumetric locking mitigation",[true, false]),
                damping = ("Select damping coefficient",[0.0, 0.1, 0.2, 0.4])
            ),
            grf   = (
                status = ("Enable Gaussian Random Field ",[true, false]),
                covariance = ("Select covariance function",["gaussian", "exponential"]),
                param = ( 
                    Iₓ = ("Select correlation length",[[1.0,1.0,1.0], [2.5,2.5,2.5], [5.0,5.0,5.0]]), 
                    Nₕ = ("Select number of points",[1000, 5000, 10000]), 
                    kₘ = ("Select maximum wavenumber",[50, 100, 200]),
                ),
            ),
            plast = (
                status = ("Enable plasticity",[true, false]),
                constitutive = ("Select constitutive model",["DP", "MC", "VM"]),
            ),
            nonloc = (
                status = ("Enable nonlocal effects",[true, false]),
                ls = ("Select length scale",[0.125 ,0.25, 0.5]),
            ),
            plot  = (
                status = ("Enable plotting",[true, false]),
                freq   = ("Select plot frequency",[0.1, 0.5, 1.0, 5.0]),
                dpi    = ("Select plot resolution",[100, 300, 500, 1000]),
                what   = ("Select plot variable(s)",[
                    (;mpts=(name="P"       ,cblim=nothing  ,),),
                    (;mpts=(name="epII"    ,cblim=(0.0,1.5),),),
                    (;mpts=(name="n"       ,cblim=(0.0,1.0),),),
                    (;mpts=(name="J"       ,cblim=(0.5,1.5),),),
                    (;mpts=(name="vol0"    ,cblim=nothing  ,),),
                ]),
            ),
            perf  = (
                status = ("Enable performance monitoring",[true, false]),
            ),
            backend = (
                select = ("Select backend",["host", "cuda", "rocm"]),
                distributed = ("Enable distributed computing",[true, false])
            ),
        )
end


function process_cli_option(value, default_val, key_path=())
    if isa(value, NamedTuple)
        # Check if this NamedTuple has a status field
        if haskey(value, :status)
            # Process status first
            status_result = process_cli_option(value.status, default_val.status, (key_path..., :status))
            if status_result == false
                # Return defaults for this entire section
                return default_val
            end
            # Process remaining fields (excluding status since it's already processed)
            other_keys = filter(k -> k != :status, keys(value))
            processed = Dict{Symbol,Any}(:status => status_result)
            for k in other_keys
                processed[k] = process_cli_option(getfield(value, k), getfield(default_val, k), (key_path..., k))
            end
            return NamedTuple{Tuple(keys(value))}(getindex(processed, k) for k in keys(value))
        end
        return NamedTuple{keys(value)}(process_cli_option(getfield(value, k), getfield(default_val, k), (key_path..., k)) for k in keys(value))
    elseif isa(value, Tuple) && length(value) == 2
        prompt, vals = value
        # Special handling for plot.what: allow multi-selection of NamedTuples, keep original definition
        if key_path[end] == :what && occursin("plot", string(key_path))
            menu = MultiSelectMenu([string(v) for v in vals], pagesize=length(vals))
            idxs = request(prompt * " (Space to select, Enter to confirm):", menu)
            selected = [vals[i] for i in idxs]
            return selected
        else
            menu = RadioMenu(string.(vals), pagesize=length(vals))
            return vals[request(prompt, menu)]
        end
    else
        return value
    end
end

"""
    cli(::Type{Instruction}, dim::Number=2) -> Instruction

Interactively build a configuration using an automated CLI menu system.
The function recursively processes all configuration options and presents 
RadioMenu selections for each option.

# Arguments
- `::Type{Instruction}`: The type to build interactively
- `dim::Number=2`: Spatial dimension (2D or 3D)

# Configuration Format
Options are defined as tuples `("prompt text", [values])` or nested NamedTuple structures.
The function automatically:
- Creates RadioMenu for tuple-formatted options
- Recursively processes nested NamedTuple configurations
- Handles arbitrary nesting depth (e.g., `grf.param.Iₓ`)
- Skips sub-options when `status = false` is selected

# Returns
- `Instruction`: Validated configuration processed through `kwargser(kwargs, Instruction; dim=dimension)`

# Example
```julia
instr = cli(Instruction; dimension::Number=2)
# Presents interactive menus for:
# - dtype (Float32/Float64)
# - basis functions
# - GRF parameters (if enabled)
# - plotting options (if enabled)
# - performance monitoring
# - backend selection
```
"""
function cli(; ui::Bool=false)    
    default_config = get_default(Instruction)
    if !ui
        kwargs = Dict{Any,Any}()
        for (K, V) ∈ pairs(default_config)
            kwargs[K] = process_cli_option(V, getfield(default_config, K), (K,))
        end
    else
        # Process all options recursively
        @info """
        Starting interactive ϵlastσPlasm 👻 v$(get_version()) command-line interface (cli):
        - use arrow keys to navigate and Enter to select.
        """
        kwargs    = Dict{Any,Any}()
        top_level = collect(keys(get_option(Instruction)))
        menu      = MultiSelectMenu(string.(top_level), pagesize=length(top_level))
        selection = Set(top_level[i] for i ∈ request("Select simulation configuration option(s):", menu))
        for (K, V) ∈ pairs(get_option(Instruction))
            if K ∈ selection
                kwargs[K] = process_cli_option(V, getfield(default_config, K), (K,))
            else
                kwargs[K] = getfield(default_config, K)
            end
        end
    end
    return kwargs
end
"""
    kwargser(kwargs::Any, ::Type{Instruction}; dim::Number=2) -> Instruction

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
instr = kwargser(Dict(:dtype => 64, :locking => true), Instruction; dim=3)
println(instr.cairn.shpfun)  # prints the initialized shape function kernel

# To override nested fields, use a nested NamedTuple or Dict:
instr = kwargser((;fwrk = (deform="infinitesimal", trsfr="tpic", locking=false)), Instruction)
println(instr.fwrk.deform)  # prints "infinitesimal"
```

# Notes
- Merges `kwargs` with default configuration
- Sets precision and attaches kernel functions for shape functions, mapping, and constitutive updates
- Warns about any unused keyword arguments
"""
function kwargser(kwargs::Any, ::Type{Instruction}; dim::Number=2)
    ref    = get_default(Instruction)
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
    return Instruction{instr[:dtype][:T0]...,typeof(dim)}(
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