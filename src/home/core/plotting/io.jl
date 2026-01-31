export get_plot_field, save_plot

"""
    get_plot_field(mpts, mesh, opts; P=Vector{Any}(), dpi=100) -> DisplayObject

Generate and display plots for all fields specified in `opts[:what]`.

# Arguments
- `mpts`: Material point data structure.
- `mesh`: Mesh data structure.
- `opts`: Named tuple or dictionary specifying what to plot and plot settings (e.g., `what`, `dims`).
- `P`: (Optional) Vector to collect plot objects.
- `dpi`: (Optional) DPI resolution for plots (default: 100).

# Returns
- `DisplayObject`: The displayed plot figure.

# Example
```julia
opts = (;what=["P", "epII"], dims=(500,500))
get_plot_field(mpts, mesh, opts)
```
"""
function get_plot_field(mpts,mesh,opts; P::Vector{Any}=[], dpi::Int=100) 
    # scale
    scale = length(opts[:what]) 
    sx,sy = opts[:dims][1],scale*opts[:dims][2]
    # plotting
    config_plot(); opts.backend
    for variable in opts[:what]
        k = first(keys(variable))
        entry = variable[k]
        if k == :mpts
            p0 = what_plot_field(mpts, (;opts..., what=entry.name, cblim=entry.cblim))
        elseif k == :mesh
            p0 = what_plot_field(mesh, (;opts..., what=entry.name, cblim=entry.cblim))
        else
            error("Unknown plot type: $k")
        end
        push!(P, p0)
    end
    # display fig
    fig   = display(plot(P...;layout=(scale,1),size=(sx,sy)))
    return fig
end

"""
    save_plot(opts) -> Nothing

Save the current plot to a file specified in `opts.file` and log the output path.

# Arguments
- `opts`: Named tuple or dictionary containing the `file` path for saving the plot.

# Returns
- `Nothing`. Logs the generated file path.

# Example
```julia
opts = (;file="output.png")
save_plot(opts)
```
"""
function save_plot(opts)
    savefig(opts.file)
    return @info "Generated fig: \n\e[32m+ $(trunc_path(opts.file))\e[0m"
end
