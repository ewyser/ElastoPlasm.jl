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
            ("nel,np",(round(Int64,prod(mesh.prprt.nel[1:end-1])),mpts.nmp)),
            ("iteration(s)",it),
            ]
    return vals
end

"""
    make_progress(total::Int; desc::String="Processing", dt::Real=0.1, barlen::Int=10) -> Progress

Create a configured progress bar for tracking iterations.

# Arguments
- `total::Int`: Total number of steps/iterations
- `desc::String`: Description text to display (default: "Processing")
- `dt::Real`: Minimum time interval between updates in seconds (default: 0.1)
- `barlen::Int`: Length of the progress bar (default: 10)

# Returns
- `Progress`: Configured progress bar ready to use with `next!()`

# Example
```julia
prog = make_progress(100; desc="Simulating")
for i in 1:100
    # do work...
    next!(prog; desc="Simulating \$i/100...")
end
```
"""
function make_progress(total::Int; desc::String="Processing", dt::Real=0.1, barlen::Int=10)
    prog = Progress(total; dt=dt, desc=desc*" 0/$total...", barlen=barlen)
    sleep(1.01*dt)
    update!(prog, 0)
    return prog::Progress
end