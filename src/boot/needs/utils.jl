"""
    superList(lists::Vector{String}; root::String=info.sys.root, tree::Bool=false) -> Vector{String}

Loads and parses a list of module directories, optionally displaying a tree structure of included files.

# Arguments
- `lists::Vector{String}`: List of directory names to include.
- `root::String=info.sys.root`: Root directory for the modules.
- `tree::Bool=false`: If true, displays a tree structure of included files.

# Returns
- `Vector{String}`: Messages summarizing the inclusion status for each directory.

# Example
```julia
msgs = superList(["src/boot", "src/home"])
println.(msgs)
```
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
    superInc(DIR::String) -> NamedTuple

Recursively includes all `.jl` files in the given directory and its subdirectories.

# Arguments
- `DIR::String`: Directory to search for Julia files.

# Returns
- `NamedTuple`: Contains `jls` (filenames) and `paths` (absolute paths).

# Example
```julia
incs = superInc("src/boot")
println(incs.jls)
```
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
    tree(sucess, prefix="\n\t", level=0, max_level=1) -> Vector{String}

Formats a list of strings into a tree-like structure for display.

# Arguments
- `sucess`: List of strings to format.
- `prefix="\n\t"`: String prefix for each line.
- `level=0`: Current tree depth.
- `max_level=1`: Maximum tree depth.

# Returns
- `Vector{String}`: Tree-formatted strings.

# Example
```julia
tree(["boot", "home"])
```
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
    rootflush(info) -> Vector{String}

Creates or flushes the output directory, removing files that do not match a given pattern.

# Arguments
- `info`: Struct containing system and MPI information.

# Returns
- `Vector{String}`: Messages about created or deleted files.

# Example
```julia
msgs = rootflush(info)
println.(msgs)
```
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
					push!(msg,"\e[31m- $(trunc_path(file))\e[0m")
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
- `anchor`: The folder name from which you want to keep the rest of the path.

# Returns
- `String`: Truncated path string.

# Example
```julia
trunc_path("C:/Users/lili8/Documents/GitHub/ElastoPlasm.jl/dump/slump", "ElastoPlasm.jl")
# => "ElastoPlasm.jl/dump/slump"
```
"""
function trunc_path(full_path::AbstractString; anchor::AbstractString="ElastoPlasm.jl")
	parts = splitpath(full_path)
	idx = findfirst(==(anchor), parts)
	return isnothing(idx) ? full_path : joinpath(parts[idx:end]...)
end


function ic_log(mesh,mp,time)
	# build the list of constant log lines
    logs = [
		"Summary:",
		"- elements: $(mesh.nel[end])",
		"- material points: $(mp.nmp)", 
		"- simulation time ∈ $(time.t) s:",
    ]
    # add optional lines
    if isa(time.tg,AbstractFloat)
        push!(logs, "   - gravity ramp-up: $(time.tg ) s")
    end
	if isa(time.te,AbstractFloat)
        push!(logs, "   - elastodynamic  : $(time.te ) s")
    end
	if isa(time.tep,AbstractFloat)
        push!(logs, "   - elastoplastic  : $(time.tep) s")
    end
	return join(logs,"\n")::String
end

"""
    elastoplasm_log(instr::NamedTuple) -> String

Generates a summary log string describing the current simulation configuration for ElastoPlasm.

# Arguments
- `instr::NamedTuple`: Instruction/configuration named tuple containing simulation options.

# Returns
- `String`: Multi-line string summarizing the simulation setup.

# Example
```julia
logstr = plasming_logs(instr)
println(logstr)
```
"""
function elastoplasm_log(instr; msg::String="elastodynamic")
    # build the list of log lines
    logs = [
        "Launching ϵlastσPlasm 👻 v$(get_version()):",
        "- $(nthreads()) active thread(s)",
		"- $msg workflow",
        "- $(instr[:fwrk][:deform]) strain formulation",
        "- $(instr[:basis][:which]) calculation cycle",
    ]
    # add optional lines only if the corresponding flags are true
    if instr[:fwrk][:locking]
        push!(logs, "- F-bar locking mitigation")
    end

    if instr[:nonloc][:status] && msg == "elastoplastic"
        push!(logs, "- non-local plastic regularization")
    end
    return join(logs,"\n")
end