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
@views function what_plot_field(mpts::Point{T1,T2,D,E,R},opts; dim::String="d") where {T1,T2,D<:AbstractDimension,E<:AbstractElasticity,R<:AbstractRheology}
    if D <: TwoDimension
        dim = "2"
    elseif D <: ThreeDimension
        dim = "3"
    end
    
    if opts.what == "P"
        if size(mpts.s.σᵢ,1) == 3
            d   = -(mpts.s.σᵢ[1,:]+mpts.s.σᵢ[2,:])/2/1e3
        elseif size(mpts.s.σᵢ,1) == 6
            d   = -(mpts.s.σᵢ[1,:]+mpts.s.σᵢ[2,:]+mpts.s.σᵢ[3,:])/3/1e3
        end 
        lab   = L"$p(x_p)$"*" [kPa]"           
        tit   = "pressure"
        cb    = :viridis
    elseif opts.what == "sigxx"
        d     = mpts.s.σᵢ[1,:]/1e3
        lab   = L"$\sigma_{xx}$"*" [kPa]"
        tit   = "stress tensor x-component"
        cb    = :viridis
    elseif opts.what == "epsxx"
        d     = mpts.s.ϵᵢⱼ[1,1,:]
        lab   = L"$\epsilon_{xx}$"*" [-]"
        tit   = "strain tensor x-component"
        cb    = :viridis      
    elseif opts.what == "epII"
        d     = mpts.s.ϵpII[1,:]
        lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"*" [-]"
        tit   = "plastic strain"
        cb    = :viridis
    elseif opts.what == "epV"
        d     = mpts.s.ϵpV
        lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
        tit   = "volumetric plastic strain"
        cb    = :seismic
    elseif opts.what == "du"
        d     = sqrt.(mpts.s.u[1,:].^2+mpts.s.u[2,:].^2)
        lab   = L"$\Delta u$"*" [m]"
        tit   = "displacement"
        cb    = :viridis
    elseif opts.what == "z0"
        d     = mpts.z₀
        lab   = L"$z_p(t_0)$"*" [m]"
        tit   = "initial vertical position"
        cb    = palette(:grayC,5) 
    elseif opts.what == "coh0"
        d     = mpts.s.c₀./1e3
        lab   = L"c_0(x_p)"*" [kPa]"
        tit   = "initial cohesion field"
        cb    = :vik
    elseif opts.what == "phi0"
        d     = mpts.s.ϕ
        lab   = L"$\phi_0(x_p)$"*" [rad]"
        tit   = "initial friction angle"
        cb    = :viridis
    elseif opts.what == "n0"
        d     = mpts.n
        lab   = L"$\n_{0}(x_p)$"*" [-]"
        tit   = "initial porosity"
        cb    = :viridis
    elseif opts.what == "n"
        d     = mpts.n
        lab   = L"n(x_p)"*" [-]"
        tit   = "porosity"
        cb    = :seismic
    elseif opts.what == "J"
        d     = mpts.J 
        lab   = L"J(x_p)"*" [-]"
        tit   = "deformation determinant"
        cb    = :seismic
    elseif opts.what == "ms"
        d     = mpts.s.ρ.*mpts.Ω
        lab   = L"m_{s}(x_p)"*" [-]"
        tit   = "solid mass"
        cb    = :viridis
    elseif opts.what == "vol0"
        d     = mpts.Ω₀
        lab   = L"\Omega_{s}(x_p)"*" [m"*L"^{%$dim}"*"]" 
        tit   = "initial solid volume"
        cb    = :viridis   
    elseif opts.what == "vol"
        d     = mpts.Ω
        ndim  = size(mpts.x, 1)
        lab   = L"\Omega_{s}(x_p)"*" [m"*L"^{%$dim}"*"]" 
        tit   = "solid volume"
        cb    = :viridis    
    elseif opts.what == "T"
        d     = mpts.t.T
        lab   = L"T(x_p)"*" [K]"
        tit   = "temperature"
        cb    = :thermal
    elseif opts.what == "q"
        d     = sqrt.(mpts.t.q[1,:].^2 .+ mpts.t.q[2,:].^2)
        lab   = L"$q(x_p)$"*" [K]"
        tit   = "heat flux"
        cb    = :viridis             
    else        
        throw(error("UndefinedPlotOption: $(opts.what)"))
    end

    if isnothing(opts.cblim)
        dmin, dmax = minimum(d), maximum(d)
        if dmin == dmax
            clim = (dmin - 0.1 * abs(dmin) - 1e-10, dmax + 0.1 * abs(dmax) + 1e-10)
        else
            clim = (dmin, dmax)
        end
    else
        clim = opts.cblim
    end
    # plotting
    p = plot(
        if size(mpts.x,1) == 2
            mpts.x[1,:],mpts.x[2,:]
        elseif size(mpts.x,1) == 3
            mpts.x[1,:],mpts.x[3,:]
        end,
        seriestype  = :scatter,
        marker_z    = d,
        xlabel      = L"$x-$direction"*" [m]",
        ylabel      = L"$z-$direction"*" [m]",
        label       = lab,
        color       = cb,
        clim        = clim,
        xlim        = opts.xlim,
        ylim        = opts.ylim,
        title       = "$tit, at $(opts.tit)",
        aspect_ratio= 1,
        size        = opts.dims,
    )
    return p
end
@views function what_plot_field(mesh::Mesh,opts)
    if opts.what == "v"
        d     = sqrt.(mesh.s.v[1,:].^2 .+ mesh.s.v[2,:].^2)
        lab   = L"$v(x_n)$"*" [m/s]"
        tit   = "nodal solid velocity"
        cb    = :viridis
    elseif opts.what == "vx"
        d     = mesh.s.v[1,:]
        lab   = L"$v_x(x_n)$"*" [m/s]"
        tit   = "nodal solid x-velocity"
        cb    = :vik
    elseif opts.what == "vz"
        d     = mesh.s.v[2,:]
        lab   = L"$v_z(x_n)$"*" [m/s]"
        tit   = "nodal solid z-velocity"
        cb    = :vik
    elseif opts.what == "m"
        d     = mesh.s.mᵢ
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
    
    if size(mesh.x,1) == 2
        x = mesh.x[1,:][1:mesh.prprt.nno[2]:end]
        z = mesh.x[2,:][1:mesh.prprt.nno[2]    ]
    elseif size(mesh.x,1) == 3
        x = mesh.x[1,:][1:mesh.prprt.nno[3]:end]
        z = mesh.x[3,:][1:mesh.prprt.nno[3]    ]
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