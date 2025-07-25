"""
    superList(lists::Vector{String}; root::String=info.sys.root)

Description:
---
A description

Example:
---
An example
```julia

julia> 
```

Note:
---
A note
"""
function superList(lists::Vector{String}; root::String=info.sys.root, tree::Bool=false,)
	sucess = ["Welcome to ϵlastσPlasm 👻 \nsuperInc() jls parser:"]
	for (k,dir) ∈ enumerate(lists)
		incs = superInc(joinpath(root,dir))
		if tree
			relative = Vector{String}()
			for path ∈ incs.paths
				push!(relative,last(split(path,"/ElastoPlasm.jl")))
			end
			push!(sucess,join(tree(relative)))	
		end
		push!(info.sys.lib,("$(dir)" => incs.jls))
		push!(sucess      ,"✓ $(dir)"  )
	end
	return sucess
end
"""
    superInc(DIR::String)

Description:
---
A description

Example:
---
An example
```julia

julia> 
```

Note:
---
A note
"""
function superInc(DIR::String)
	jls,paths = Vector{String}(),Vector{String}()
	for (root, dirs, files) ∈ walkdir(DIR)
		for file ∈ files
			f = joinpath(root,file)
			if last(splitext(f)) == ".jl" 
				include(f); push!(jls,file); push!(paths,f)
			end
		end
	end
	return incs = (; jls = jls, paths = paths,)
end
"""
    tree(sucess, prefix="\n\t", level=0, max_level=1)

Description:
---
A description

Example:
---
An example
```julia

julia> 
```

Note:
---
A note
"""
function tree(sucess, prefix="\n\t", level=0, max_level=1)
    if level > max_level
        return nothing
    end
    n,printout = length(sucess),[]
    for (i, name) ∈ enumerate(sucess)
        connector = i == n ? "└── " : "├── "
		push!(printout,prefix*connector*name)
    end
	return printout
end
"""
    rootflush(info)

Description:
---
A description

Example:
---
An example
```julia

julia> 
```

Note:
---
A note
"""
function rootflush(info)
	if !isdir(info.sys.out)
		msg = ["Creating:\n+ $(trunc_path(info.sys.out))"]
		mkdir(info.sys.out) 
	else
		msg,files = ["Nothing to flush at /dump"],readdir(info.sys.out;join=true)
		if !isempty(files)
			msg = ["Flushing:"]
			for file ∈ files
				if !occursin(info.mpi.glob,file) 
					rm(file,recursive=true)  
					push!(msg,"✗ $(trunc_path(file))")
				end
			end
		end
	end
	return msg
end

"""
    trunc_path(full_path::AbstractString; anchor::AbstractString="ElastoPlasm.jl") -> String

Returns the subpath of `full_path` starting from the directory name `anchor`.

# Arguments
- `full_path`: The full absolute or relative path.
- `anchor`: The folder name from which you want to keep the rest of the path (e.g., "ElastoPlasm.jl").

# Returns
- A truncated path string like `"ElastoPlasm.jl/dump/slump"`.

# Example
```julia
truncate_path_from("C:/Users/lili8/Documents/GitHub/ElastoPlasm.jl/dump/slump", "ElastoPlasm.jl")
# => "ElastoPlasm.jl/dump/slump"
"""
function trunc_path(full_path::AbstractString; anchor::AbstractString="ElastoPlasm.jl")
	parts = splitpath(full_path)
	idx = findfirst(==(anchor), parts)
	return isnothing(idx) ? full_path : joinpath(parts[idx:end]...)
end

"""
    plasming_logs(instr::NamedTuple) -> String

Generates a summary log string describing the current simulation configuration for ElastoPlasm.

# Arguments
- `instr::NamedTuple`: Instruction/configuration named tuple containing simulation options.

# Returns
- `String`: A multi-line string summarizing the simulation setup, including active threads, strain formulation, calculation cycle, and optional features such as F-bar locking mitigation and non-local plastic regularization.

# Example
```julia
logstr = plasming_logs(instr)
println(logstr)
```
"""
function plasming_logs(instr)
    # Build the list of log lines
    logs = [
        "Launching ϵlastσPlasm 👻 v$(get_version()):",
        "- $(nthreads()) active thread(s)",
        "- $(instr[:fwrk][:deform]) strain formulation",
        "- $(instr[:basis][:which]) calculation cycle",
    ]
    # Add optional lines only if the corresponding flags are true
    if instr[:fwrk][:locking]
        push!(logs, "- F-bar locking mitigation")
    end

    if instr[:nonloc][:status]
        push!(logs, "- non-local plastic regularization")
    end
    return join(logs,"\n")
end