export cli, process_cli_option

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
    elseif isa(value, Tuple) && length(value) == 2
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
- `Dict`: Configuration dictionary ready to be passed to `kwargser()`

# Example
```julia
kwargs = cli(ui=true)
instr = kwargser(kwargs; dim=2)
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
