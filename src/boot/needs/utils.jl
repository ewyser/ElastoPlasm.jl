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
	root = basename(info.sys.out)
	if !isdir(info.sys.out)
		msg = ["Creating /dump at:\n+ $(root)"]
		mkdir(info.sys.out) 
	else
		msg,files = ["Nothing to flush at /dump"],readdir(info.sys.out;join=true)
		if !isempty(files)
			msg = ["Flushing /dump:"]
			for file ∈ files
				if !occursin(info.mpi.glob,file) 
					rm(file,recursive=true)  
					push!(msg,"✗ $(basename(file))")
				end
			end
		end
	end
	return msg
end