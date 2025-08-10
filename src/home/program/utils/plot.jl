export  save_plot,get_plot_field,what_plot_field

"""
    what_plot_field(mpts, mesh, opts)

Select and prepare a field for plotting from the material point and mesh data, based on the `opts` specification.

# Arguments
- `mpts`: Material point data structure.
- `mesh`: Mesh data structure.
- `opts`: Named tuple or dictionary specifying what to plot and plot settings (e.g., `what`, `dims`, `tit`, `backend`).

# Returns
- `Plot`: Plot object for the selected field.

# Example
```julia
opts = (;what="P", dims=(500,250), tit="", backend=gr())
p = what_plot_field(mpts, mesh, opts)
display(p)
```

# Notes
- Supports fields: pressure (`P`), plastic strain (`epII`), volumetric plastic strain (`epV`), displacement (`du`), initial vertical position (`z0`), initial cohesion (`coh0`), and initial friction angle (`phi0`).
- Throws an error if the requested field is not defined.
"""
@views function what_plot_field(mpts,mesh,opts)
    if opts.what == "P"
        if size(mpts.s.σᵢ,1) == 3
            d   = -(mpts.s.σᵢ[1,:]+mpts.s.σᵢ[2,:])/2/1e3
            lab = L"$p=-\sigma_{ii,p}/2$"*" [kPa]"
        elseif size(mpts.s.σᵢ,1) == 6
            d   = -(mpts.s.σᵢ[1,:]+mpts.s.σᵢ[2,:]+mpts.s.σᵢ[3,:])/3/1e3
            lab = L"$p=-\sigma_{ii,p}/3$"*" [kPa]"
        end            
        tit   = "pressure"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (minimum(d),maximum(d))
        end
    elseif opts.what == "epII"
        d     = mpts.s.ϵpII[1,:]
        lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"*" [-]"
        tit   = "plastic strain"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif opts.what == "epV"
        d     = mpts.s.ϵpV
        lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
        tit   = "volumetric plastic strain"
        cb    = :seismic
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (-maximum(abs.(d)),maximum(abs.(d)))
        end
    elseif opts.what == "du"
        d     = sqrt.(mpts.s.u[1,:].^2+mpts.s.u[2,:].^2)
        lab   = L"$\Delta u$"*" [m]"
        tit   = "displacement"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif opts.what == "z0"
        d     = mpts.z₀
        lab   = L"$z_p(t_0)$"*" [m]"
        tit   = "initial vertical position"
        cb    = palette(:grayC,5) 
        cblim = (0.0,maximum(d))
    elseif opts.what == "coh0"
        d     = mpts.s.c₀./1e3
        lab   = L"c_0(x_p)"*" [kPa]"
        tit   = "initial cohesion field"
        cb    = :vik
        coh0  = sum(d)/length(d)
        cblim = (coh0-coh0/2,coh0+coh0/2)
    elseif opts.what == "phi0"
        d     = mpts.s.ϕ
        lab   = L"$\phi_0(x_p)$"*" [rad]"
        tit   = "initial friction angle"
        cb    = :viridis
        cblim = (minimum(d),maximum(d)) 
    elseif opts.what == "n0"
        d     = mpts.n
        lab   = L"$\n_{0}(x_p)$"*" [-]"
        tit   = "initial porosity"
        cb    = :viridis
        cblim = (0.0,1.0) 
    elseif opts.what == "n"
        d     = mpts.n
        lab   = L"$\n(x_p)$"*" [-]"
        tit   = "porosity"
        cb    = :viridis
        cblim = (0.0,1.0) 
    else
        throw(error("UndefinedPlotOption: $(opts.what)"))
    end

    # plotting
    p = plot(
        if mesh.dim == 2
            mpts.x[1,:],mpts.x[2,:]
        elseif mesh.dim == 3
            mpts.x[1,:],mpts.x[3,:]
        end,
        seriestype  = :scatter,
        marker_z    = d,
        xlabel      = L"$x-$direction"*" [m]",
        ylabel      = L"$z-$direction"*" [m]",
        label       = lab,
        color       = cb,
        clim        = cblim,
        ylim        = (-10.0,20.0),
        title       = "$tit $(opts.tit)",
        aspect_ratio= 1,
        size        = opts.dims,
    )
    return p
end

"""
    get_plot_field(mpts, mesh, opts; P=Vector{{Any}}())

Generate and display plots for all fields specified in `opts[:what]`.

# Arguments
- `mpts`: Material point data structure.
- `mesh`: Mesh data structure.
- `opts`: Named tuple or dictionary specifying what to plot and plot settings (e.g., `what`, `dims`).
- `P`: (Optional) Vector to collect plot objects.

# Returns
- `Nothing`. Displays the plot(s).

# Example
```julia
opts = (;what=["P", "epII"], dims=(500,500))
get_plot_field(mpts, mesh, opts)
```
"""
function get_plot_field(mpts,mesh,opts; P::Vector{Any}=[]) 
    # plotting
    config_plot(); opts.backend
    for (k,variable) ∈ enumerate(opts[:what])
        opts = (;opts...,what=variable)
        p0   = what_plot_field(mpts,mesh,(;opts...,what=variable))
        push!(P,p0)
    end
    scale = length(P) 
    sx,sy = opts[:dims][1],scale*opts[:dims][2]
    fig   = display(plot(P...;layout=(scale,1),size=(sx,sy)))
    return fig
end

"""
    save_plot(opts)

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