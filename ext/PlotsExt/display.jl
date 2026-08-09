export what_plot_field, get_plot_field, save_plot

"""
    what_plot_field(mpts, opts) -> Plot

Select and prepare a material point field for plotting based on the `opts` specification.

# Arguments
- `mpts`: Material point data structure.
- `opts`: Named tuple or dictionary specifying what to plot and plot settings (e.g., `what`, `dims`, `tit`, `backend`).

# Returns
- `Plot`: Plot object for the selected field.

# Example
```julia
opts = (;what="P", dims=(500,250), tit="", backend=gr())
p = what_plot_field(mpts, opts)
display(p)
```

# Notes
- Supports fields: pressure (`P`), plastic strain (`epII`), volumetric plastic strain (`epV`), displacement (`du`), initial vertical position (`z0`), initial cohesion (`coh0`), and initial friction angle (`phi0`).
- Throws an error if the requested field is not defined.
"""
function what_plot_field(mpts::Point{T1,T2,D,B,E,R},opts) where {T1,T2,D<:AbstractDimension,B<:AbstractBasis,E<:AbstractElasticity,R<:AbstractRheology}
    # Extract data using the data function from opts
    data = opts.data(mpts)
    if isnothing(opts.cblim)
        dmin, dmax = minimum(data), maximum(data)
        if dmin == dmax
            clim = (dmin - 0.1 * abs(dmin) - 1e-10, dmax + 0.1 * abs(dmax) + 1e-10)
        else
            clim = (dmin, dmax)
        end
    else
        clim = opts.cblim
    end
    # Plotting
    return plot(
        if length(mpts.x[1]) == 2
            getindex.(mpts.x, 1), getindex.(mpts.x, 2)
        elseif length(mpts.x[1]) == 3
            getindex.(mpts.x, 1), getindex.(mpts.x, 3)
        end,
        seriestype  = :scatter,
        marker_z    = data,
        xlabel      = L"$x-$direction"*" [m]",
        ylabel      = L"$z-$direction"*" [m]",
        label       = opts.label,
        color       = opts.cb,
        clim        = clim,
        xlim        = opts.xlim,
        ylim        = opts.ylim,
        title       = opts.tit,
        aspect_ratio= 1,
        size        = opts.dims,
    )
end

@views function what_plot_field(mesh::Mesh{T1,T2,D},opts) where {T1,T2,D<:AbstractDimension}
    if opts.what == "v"
        d     = [sqrt(v[1]^2 + v[2]^2) for v in mesh.s.v]
        lab   = L"$v(x_n)$"*" [m/s]"
        tit   = "nodal solid velocity"
        cb    = :viridis
    elseif opts.what == "vx"
        d     = getindex.(mesh.s.v, 1)
        lab   = L"$v_x(x_n)$"*" [m/s]"
        tit   = "nodal solid x-velocity"
        cb    = :vik
    elseif opts.what == "vz"
        d     = getindex.(mesh.s.v, 2)
        lab   = L"$v_z(x_n)$"*" [m/s]"
        tit   = "nodal solid z-velocity"
        cb    = :vik
    elseif opts.what == "m"
        d     = mesh.s.m
        lab   = L"$m(x_n)$"*" [kg]"
        tit   = "nocdal solid mass"
        cb    = :viridis
    elseif opts.what == "c"
        d     = mesh.t.cᵢ
        lab   = L"$c(x_n)$"*" [J/(kg·K)]" 
        tit   = "nocdal specific heat capacity"
        cb    = :viridis
    elseif opts.what == "T"
        d     = mesh.t.T
        lab   = L"$T(x_n)$"*" [K]"
        tit   = "nodal temperature"
        cb    = :thermal    
    elseif opts.what == "bcs"
        d     = mesh.t.bcs.status[1,:]
        lab   = L"$T(x_n)$"*" [K]"
        tit   = "nodal temperature"
        cb    = :viridis       
    else        
        throw(error("UndefinedPlotOption: $(opts.what)"))
    end
    
    if length(mesh.x[1]) == 2
        x = getindex.(mesh.x, 1)[1:mesh.prprt.nno[2]:end]
        z = getindex.(mesh.x, 2)[1:mesh.prprt.nno[2]    ]
    elseif length(mesh.x[1]) == 3
        x = getindex.(mesh.x, 1)[1:mesh.prprt.nno[3]:end]
        z = getindex.(mesh.x, 3)[1:mesh.prprt.nno[3]    ]
    end
    d = reshape(d, mesh.prprt.nno[2], mesh.prprt.nno[1])
    p = heatmap(
        x, z, d,
        xlabel       = L"$x-$direction [m]",
        ylabel       = L"$z-$direction [m]",
        label        = lab,
        color        = cb,
        clim         = opts.cblim,
        xlim         = opts.xlim,
        ylim         = opts.ylim,
        title        = "$tit, at $(opts.tit)",
        aspect_ratio = 1,
        size         = opts.dims
    )
    return p
end

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
    for variable ∈ opts[:what]
        if haskey(variable, :mpts)
            merged = merge(variable.mpts, (;tit="$(variable.mpts.name), at $(opts.tit))", xlim=opts.xlim, ylim=opts.ylim, dims=opts.dims))
            p0 = what_plot_field(mpts, merged)
        elseif haskey(variable, :mesh)
            merged = merge(variable.mesh, (;tit="$(variable.mesh.name), at $(opts.tit))", xlim=opts.xlim, ylim=opts.ylim, dims=opts.dims))
            p0 = what_plot_field(mesh, merged)
        else
            error("Unknown plot type: $variable")
        end
        push!(P, p0)
    end
    return display(plot(P...;layout=(scale,1),size=(sx,sy)))
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
