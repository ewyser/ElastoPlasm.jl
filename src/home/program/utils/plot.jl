@views function what_plot_field(mp,mesh,opts)
    if opts.what == "P"
        if size(mp.σᵢ,1) == 3
            d   = -(mp.σᵢ[1,:]+mp.σᵢ[2,:])/2/1e3
            lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}\right)/2$"
        elseif size(mp.σᵢ,1) == 6
            d   = -(mp.σᵢ[1,:]+mp.σᵢ[2,:]+mp.σᵢ[3,:])/3/1e3
            lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)/3$"
        end            
        tit   = "pressure"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (minimum(d),maximum(d))
        end
    elseif opts.what == "epII"
        d     = mp.ϵpII[1,:]
        lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"
        tit   = "plastic strain"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif opts.what == "epV"
        d     = mp.ϵpV
        lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
        tit   = "volumetric plastic strain"
        cb    = :seismic
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (-maximum(abs.(d)),maximum(abs.(d)))
        end
    elseif opts.what == "du"
        d     = sqrt.(mp.u[1,:].^2+mp.u[2,:].^2)
        lab   = L"$\Delta u$"
        tit   = "displacement"
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif opts.what == "z0"
        d     = mp.z₀
        lab   = L"$z_p(t_0)$"
        tit   = "initial vertical position"
        cb    = palette(:grayC,5) 
        cblim = (0.0,maximum(d))
    elseif opts.what == "coh0"
        d     = mp.c₀./1e3
        lab   = L"c_0(x_p) [kPa]"
        tit   = "initial cohesion field"
        cb    = :vik
        coh0  = sum(d)/length(d)
        cblim = (coh0-coh0/2,coh0+coh0/2)
    elseif opts.what == "phi0"
        d     = mp.ϕ
        lab   = L"$\phi_0(x_p)$"
        tit   = "initial friction angle"
        cb    = :viridis
        cblim = (minimum(d),maximum(d)) 
    else
        throw(error("UndefinedPlotOption: $(opts.what)"))
    end

    # plotting
    opts.backend
    p = plot(
        if mesh.dim == 2
            mp.x[1,:],mp.x[2,:]
        elseif mesh.dim == 3
            mp.x[1,:],mp.x[3,:]
        end,
        seriestype  = :scatter,
        marker_z    = d,
        xlabel      = L"$x-$direction [m]",
        ylabel      = L"$z-$direction [m]",
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

function get_plot_field(mp,mesh,opts; P::Vector{Any}=[]) 
    # plotting
    for (k,variable) ∈ enumerate(opts[:what])
        opts = (;opts...,what=variable)
        p0   = what_plot_field(mp,mesh,(;opts...,what=variable))
        push!(P,p0)
    end
    scale = length(P) 
    sx,sy = opts[:dims][1],scale*opts[:dims][2]
    P     = plot(P...;layout=(scale,1),size=(sx,sy))
    return display(P)
end

function save_plot(opts)
    savefig(opts.file)
    return @info "Generated fig: \n\e[32m+ $(trunc_path(opts.file))\e[0m"
end