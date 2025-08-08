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
    get_vals(mesh, mpts, it) -> Vector{Tuple{String, Any}}

Return a summary of mesh and simulation state for progress display or logging.

# Arguments
- `mesh`: Mesh data structure, must contain `nel`.
- `mpts`: Material point data structure, must contain `nmp`.
- `it`: Current iteration or step count.

# Returns
- `Vector{Tuple{String, Any}}`: List of (label, value) pairs summarizing the state.

# Example
```julia
vals = get_vals(mesh, mpts, 10)
println(vals)
```
"""
function get_vals(mesh,mpts,it)
    # save vals
    vals = [
            ("nel,np",(round(Int64,prod(mesh.nel[1:end-1])),mpts.nmp)),
            ("iteration(s)",it),
            ]
    return vals
end

"""
    msg(message)

Print a styled message to the console in red, bold, and blinking text if supported.

# Arguments
- `message`: The message string to display.

# Returns
- Nothing. Prints the message to the console.

# Example
```julia
msg("This is a warning message.")
```
"""
function msg(message)
    message = "└ "*message
    try
        return printstyled(message,color=:red,bold=true,blink=true)
    catch
        return printstyled(message,color=:blink)
    end
end
