export get_option, process_cli_option, cli

"""
    get_option() -> NamedTuple

Return all possible values for each configuration field as vectors.

# Returns
- `NamedTuple`: Structure matching `get_default()` but with vectors of possible values instead of single values.

# Example
```julia
opts = get_option()
println(opts.basis.which)  # ["bsmpm", "gimpm", "smpm", "mlsmpm"]
println(opts.strain.deform)  # ["finite", "infinitesimal"]
```

# Notes
- Useful for validation, UI generation, and documentation
- Fields with continuous values (like damping) show reasonable ranges
- Boolean fields are represented as `[true, false]`
"""
function get_option()
    return (

            solution = ("Select solver formulation",["explicit","implicit","dynamic relaxation"]),
            dtype = (
                T0 = ("Select arithmetic types",[(Int64,Float64),(Int32,Float32)]),
                bits = ("Select arithmetic precision",[Int64(64),Int32(32)]),
                precision = ("Select precision",["64-bit precision (or double precision)","32-bit precision (or single precision)"]),
            ),
            basis = (
                which = ("Select basis type",["bsmpm", "gimpm", "smpm", "mlsmpm"]),
                how = ("Select material point domain update",[nothing]),
                trsfr = ("Select the mapping scheme",["std", "tpic", "apic"]),
                C_pf = ("Select picflip ratio",[1.0, 0.99, 0.95]),
            ),
            strain   = (
                deform = ("Select the deformation framework",["finite", "infinitesimal"]),
            ),
            stab     = (
                locking = ("Enable volumetric locking mitigation",[true, false]),
                damping = ("Select damping coefficient",[0.0, 0.1, 0.2, 0.4]),
                musl = ("Enable musl update",[true, false]),
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
                what   = ("Select plot variable(s)", get_variable_plot_options()),
            ),
            perf  = (
                status = ("Enable optimized implementation",[true, false]),
            ),
            backend = (
                select = ("Select backend",["host", "cuda", "rocm"]),
                distributed = ("Enable distributed computing",[true, false])
            ),
        )
end


"""
    process_cli_option(value, default_val, key_path=())

Recursively process configuration options to build an interactive CLI menu system.

# Arguments
- `value`: Current option value (NamedTuple, Tuple, or scalar)
- `default_val`: Default value for this option
- `key_path`: Tuple tracking the current path in the configuration tree

# Returns
- Processed value based on user selection or default

# Notes
- Handles nested NamedTuples recursively
- Creates RadioMenu for tuple-formatted options
- Special handling for plot.what to allow multi-selection
- Respects status fields to skip disabled sections
"""
function process_cli_option(value, default_val, key_path=())
    if isa(value, NamedTuple)
        # Check if this NamedTuple has a status field
        if haskey(value, :status)
            # Process status first
            status_result = process_cli_option(value.status, default_val.status, (key_path..., :status))
            if status_result == false
                # Return defaults for this entire section
                return merge(default_val,(;status=false))
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
    elseif isa(value, Tuple) && length(value) == 2 && isa(value[1], AbstractString)
        prompt, vals = value
        # Special handling for plot.what: allow multi-selection of NamedTuples, keep original definition
        if key_path[end] == :what && occursin("plot", string(key_path))
            # Extract names from either mpts or mesh variables
            names = String[] # sorter is names = [string(haskey(v, :mpts) ? v.mpts.name : v.mesh.name) for v ∈ vals]
            for v ∈ vals
                if haskey(v, :mpts)
                    push!(names, string(v.mpts.name))
                elseif haskey(v, :mesh)
                    push!(names, string(v.mesh.name))
                else
                    push!(names, string(v))
                end
            end
            menu  = MultiSelectMenu(names, pagesize=length(names))
            idxs  = request(prompt * " (Space to select, Enter to confirm):", menu)
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
    cli(; ui::Bool=false) -> Dict

Interactively build a configuration using an automated CLI menu system.
The function recursively processes all configuration options and presents 
RadioMenu selections for each option.

# Arguments
- `ui::Bool=false`: If true, allows user to select which top-level options to configure; if false, uses defaults for all unselected options

# Configuration Format
Options are defined as tuples `("prompt text", [values])` or nested NamedTuple structures.
The function automatically:
- Creates RadioMenu for tuple-formatted options
- Recursively processes nested NamedTuple configurations
- Handles arbitrary nesting depth (e.g., `grf.param.Iₓ`)
- Skips sub-options when `status = false` is selected

# Returns
- `Dict`: Configuration dictionary ready to be splatted as kwargs into `slump_problem`/`get_solver`.

# Example
```julia
jld2 = slump_problem(L, nel; cli(ui=true)...)
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
    default_config = get_default()
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
        top_level = collect(keys(get_option()))
        menu      = MultiSelectMenu(string.(top_level), pagesize=length(top_level))
        selection = Set(top_level[i] for i ∈ request("Select simulation configuration option(s):", menu))
        for (K, V) ∈ pairs(get_option())
            if K ∈ selection
                kwargs[K] = process_cli_option(V, getfield(default_config, K), (K,))
            else
                kwargs[K] = getfield(default_config, K)
            end
        end
    end
    return kwargs::Dict{Any,Any}
end
