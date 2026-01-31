
"""
    superInc(lists::Vector{String}; root::String=info.sys.root, tree::Bool=false) -> Vector{String}

Recursively includes all `.jl` files in the given directories and their subdirectories, optionally displaying a tree structure of included files.

# Arguments
- `lists::Vector{String}`: List of directory names to include.
- `root::String=info.sys.root`: Root directory for the modules.
- `tree::Bool=false`: If true, displays a tree structure of included files.

# Returns
- `Vector{String}`: Messages summarizing the inclusion status for each directory.

# Example
```julia
msgs = superInc(["src/boot", "src/home"])
println.(msgs)
```
```
"""
function superInc(lists::Vector{String}; root::String="",lib::Dict=Dict(), tree::Bool=false)
    sucess = ["superInc() jls parser:"]
    all_errors = Dict{String, Vector{Tuple{String, Exception, Vector{Base.StackTraces.StackFrame}}}}()
    
    for (k,dir) ∈ enumerate(lists)
        dict = superDir(joinpath(root,dir))
        dir_errors = Tuple{String, Exception, Vector{Base.StackTraces.StackFrame}}[]
        
        # Collect and include all .jl files in this subtree
        function collect_and_include_jls(d, path="")
            files = String[]
            for (k,v) ∈ d
                if isa(v, Dict)
                    append!(files, collect_and_include_jls(v, joinpath(path, k)))
                elseif endswith(k, ".jl")
                    try
                        include(v)
                        push!(files, k)
                    catch e
                        bt = Base.catch_backtrace()
                        stack = Base.stacktrace(bt)
                        rel_path = joinpath(path, k)
                        push!(dir_errors, (rel_path, e, stack))
                        @warn "Failed to include $(rel_path)" exception=(e, bt)
                    end
                end
            end
            return files
        end
        
        jls_files = collect_and_include_jls(dict)
        
        if !isempty(dir_errors)
            all_errors[dir] = dir_errors
            push!(sucess, "⚠ $(dir) ($(length(jls_files)) loaded, $(length(dir_errors)) failed)")
        else
            push!(sucess, "✓ $(dir)")
        end
        
        if tree
            push!(sucess, join(tree(collect(keys(dict)))))
        end
        
        # Store the nested dictionary and errors for each directory
        push!(lib, (dir => Dict("files" => dict, "errors" => dir_errors)))
    end
    
    # Store consolidated errors in lib
    if !isempty(all_errors)
        push!(lib, ("_errors" => all_errors))
    end
    
    return sucess
end

"""
    superDir(DIR::String) -> Dict

Recursively builds a nested Dict representing the directory structure and `.jl` files.

# Arguments
- `DIR::String`: Directory to search for Julia files and subdirectories.

# Returns
- `Dict`: Nested dictionary where keys are directory or file names. `.jl` files map to their absolute paths.

# Example
```julia
tree = superDict("src/boot")
println(tree)
```
"""
function superDir(DIR::String)
    d = Dict{String,Any}()
    for entry ∈ readdir(DIR; join=true)
        name = splitpath(entry)[end]
        if isdir(entry)
            d[name] = superDir(entry)
        elseif endswith(name, ".jl")
            d[name] = entry
        end
    end
    return d
end