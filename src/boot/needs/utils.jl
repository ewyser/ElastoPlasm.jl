
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
"""
function superInc(lists::Vector{String}; root::String=info.sys.root, tree::Bool=false,)
    sucess = ["superInc() jls parser:"]
    for (k,dir) ∈ enumerate(lists)
        dict = superDir(joinpath(root,dir))
        # Collect and include all .jl files in this subtree
        function collect_and_include_jls(d)
            files = String[]
            for (k,v) ∈ d
                if isa(v, Dict)
                    append!(files, collect_and_include_jls(v))
                elseif endswith(k, ".jl")
                    include(v)
                    push!(files, k)
                end
            end
            return files
        end
        jls_files = collect_and_include_jls(dict)
        if tree
            push!(sucess, join(tree(collect(keys(dict)))))
        end
        # Store the nested dictionary for each directory for more granularity
        push!(info.sys.lib, ("$(dir)" => dict))
        push!(sucess, "✓ $(dir)")
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

"""
    get_version() -> String

Return the current project version as a string, as specified in the Julia project file.

# Returns
- `String`: The project version.

# Example
```julia
v = get_version()
println(v)
```
"""
function get_version()
    return string(Pkg.project().version)
end

"""
    welcome_log(; greeting::String="Welcome to ϵlastσPlasm 👻 v$(get_version())")

Prints a styled welcome message to the console, highlighting "Welcome" and vertical bars in green and bold.

# Arguments
- `greeting::String`: The greeting message to display at the top (default: "Welcome to ϵlastσPlasm 👻 v$(get_version())").

# Returns
- `Nothing`. Prints the welcome message to the console.

# Example
```julia
welcome_log()
welcome_log(greeting="Hello from ElastoPlasm!")
```
"""
function welcome_log(; greeting::String="Welcome to ϵlastσPlasm 👻 v$(get_version())", showcase::String = "slumping") 
    printstyled("┌ $greeting\n", color=:green, bold=true)
    printstyled("│", color=:green, bold=true); println(" New comer ? Try $(showcase) out")
    if showcase == "slumping"
        printstyled("│", color=:green, bold=true); println("   L,nel  = [64.1584,64.1584/4.0],[40,10];")
        printstyled("│", color=:green, bold=true); println("   ic,cfg = ic_slump(L,nel);")
        printstyled("└", color=:green, bold=true); println("   out    = slump(ic,cfg; workflow=\"all-in-one\");")
    elseif showcase == "collapsing"
        printstyled("│", color=:green, bold=true); println("   plot      = (; status=true, freq=1.0, what=[\"sigxx\"], dims=(500.0,250.0) );")
        printstyled("│", color=:green, bold=true); println("   fwrk      = (; deform = \"finite\",trsfr = \"musl\",locking = false,damping = 0.0);")
        printstyled("│", color=:green, bold=true); println("   nel       = [5,10];")
        printstyled("│", color=:green, bold=true); println("   ν,E,ρ0,l0 = 0.0,1.0e4,80.0,50.0;")
        printstyled("│", color=:green, bold=true); println("   ic, cfg   = ic_collapse(nel, ν, E, ρ0, l0; plot,fwrk);")
        printstyled("└", color=:green, bold=true); println("   out       = collapse(ic, cfg);")
    else
        printstyled("└", color=:green, bold=true); println("   ...$(showcase) ?!? \e[5m¯\\_(ツ)_/¯\e[0m")
    end

    return nothing
end

"""
    ic_log(mesh, mpts, time) -> String

Generates a summary log string for the initial simulation setup, including mesh, material points, and time parameters.

# Arguments
- `mesh`: Mesh object containing element information.
- `mpts`: Material point object containing number of points.
- `time`: Named tuple with time parameters (`t`, `te`, `tg`, `tep`).

# Returns
- `String`: Multi-line summary of the simulation setup.

# Example
```julia
summary = ic_log(mesh, mpts, time)
println(summary)
```
"""
function ic_log(mesh,mpts,time,instr)
    # build the list of constant log lines
    logs = [
        "Summary: ",
        "- $(instr.dtype.precision)",
        "- elements: $(mesh.nel[end])",
        "- material points: $(mpts.nmp)", 
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
    elastoplasm_log(instr::NamedTuple; msg::String="elastodynamic") -> String

Generates a summary log string describing the current simulation configuration for ElastoPlasm.

# Arguments
- `instr::NamedTuple`: Instruction/configuration named tuple containing simulation options.
- `msg::String`: (Optional) Workflow name, default is "elastodynamic".

# Returns
- `String`: Multi-line string summarizing the simulation setup.

# Example
```julia
logstr = elastoplasm_log(instr)
println(logstr)
```
"""
function elastoplasm_log(instr; msg::String="elastodynamic")
    # build the list of log lines
    logs = [
        "Launching ϵlastσPlasm 👻 v$(get_version()):",
        "└ $(nthreads()) active thread(s)",
        "- $msg workflow",
        "- $(instr[:fwrk][:deform]) strain formulation",
        "- $(instr[:basis][:which]) calculation cycle",
    ]
    # add optional lines only if the corresponding flags are true
    if instr[:fwrk][:locking]
        push!(logs, "- F-bar locking mitigation")
    end
    if instr[:nonloc][:status] && msg ∈ ["elastoplastic","all-in-one"]
        push!(logs, "- non-local plastic regularization")
    end
    return join(logs,"\n")
end

"""
    exit_log(message::String)

Print a styled exit message to the console, using green, bold, and blinking text if supported by the terminal.

# Arguments
- `message::String`: The message to display at program exit.

# Returns
- `Nothing`. Prints the message to the console.

# Example
```julia
exit_log("Simulation finished successfully.")
```
"""
function exit_log(message::String)
    try
        return printstyled("└ $message",color=:green,bold=true,blink=true)
    catch
        return printstyled("└ $message",color=:blink)
    end
end