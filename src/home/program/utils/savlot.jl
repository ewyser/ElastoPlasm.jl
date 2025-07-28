@views function whichplot(what,temp,mp,mesh,instr)
    type = what
    
    if type == "P"
        if size(mp.σᵢ,1) == 3
            d   = -(mp.σᵢ[1,:]+mp.σᵢ[2,:])/2/1e3
            lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}\right)/2$"
        elseif size(mp.σᵢ,1) == 6
            d   = -(mp.σᵢ[1,:]+mp.σᵢ[2,:]+mp.σᵢ[3,:])/3/1e3
            lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)/3$"
        end            
        tit   = "pressure, "*temp
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (minimum(d),maximum(d))
        end
    elseif type == "epII"
        d     = mp.ϵpII[1,:]
        lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"
        tit   = "plastic strain, "*temp
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif type == "epV"
        d     = mp.ϵpV
        lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
        tit   = "volumetric plastic strain, "*temp
        cb    = :seismic
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (-maximum(abs.(d)),maximum(abs.(d)))
        end
    elseif type == "du"
        d     = sqrt.(mp.u[1,:].^2+mp.u[2,:].^2)
        lab   = L"$\Delta u$"
        tit   = "displacement, "*temp
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif type == "z0"
        d     = mp.z₀
        lab   = L"$z_p(t_0)$"
        tit   = "initial vertical position, "*temp
        cb    = palette(:grayC,5) 
        cblim = (0.0,maximum(d))
    elseif type == "coh0"
        d     = mp.c₀./1e3
        lab   = L"c_0(x_p) [kPa]"
        tit   = "initial cohesion field, "*temp
        cb    = :vik
        coh0  = sum(d)/length(d)
        cblim = (coh0-coh0/2,coh0+coh0/2)
    else
        throw(error("UndefinedPlotOption: $(type)"))
    end
    # plot option
    ms = 2.25
    ms = 0.4*instr[:plot][:dims][1]/mesh.nel[1]
    # plotting
    gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,)
    p0 = plot(
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
        title       = tit,
        aspect_ratio= 1,
    )
    return p0
end

@views function savlot(mp,mesh,t,instr) 
    if !instr[:plast][:status]
        mp.z₀.= mp.x[end,:] 
    end
    if instr[:plot][:status]
        P,tit = [],L"$t = $"*string(round(t,digits=1))*" [s]"
        for (k,what) ∈ enumerate(instr[:plot][:what])
            p0 = whichplot(what,tit,mp,mesh,instr)
            push!(P,p0)
        end
        scale = length(P) 
        sx,sy = instr[:plot][:dims][1],scale*instr[:plot][:dims][2]
        display(plot(P...;layout=(scale,1),size=(sx,sy))) 
    else
        nothing
    end
    return nothing
end