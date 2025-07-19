@views function whichplot(what,temp,mpD,meD,instr)
    type = what
    
    if type == "P"
        if size(mpD.σᵢ,1) == 3
            d   = -(mpD.σᵢ[1,:]+mpD.σᵢ[2,:])/2/1e3
            lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}\right)/2$"
        elseif size(mpD.σᵢ,1) == 6
            d   = -(mpD.σᵢ[1,:]+mpD.σᵢ[2,:]+mpD.σᵢ[3,:])/3/1e3
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
        d     = mpD.ϵpII
        lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"
        tit   = "plastic strain, "*temp
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif type == "epV"
        d     = mpD.ϵpV
        lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
        tit   = "volumetric plastic strain, "*temp
        cb    = :seismic
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (-maximum(abs.(d)),maximum(abs.(d)))
        end
    elseif type == "du"
        d     = sqrt.(mpD.u[:,1].^2+mpD.u[:,2].^2)
        lab   = L"$\Delta u$"
        tit   = "displacement, "*temp
        cb    = :viridis
        if minimum(d) == maximum(d)
            cblim = (-1.0,1.0)
        else
            cblim = (0.0,maximum(d))
        end
    elseif type == "z0"
        d     = mpD.z₀
        lab   = L"$z_p(t_0)$"
        tit   = "initial vertical position, "*temp
        cb    = palette(:grayC,5) 
        cblim = (0.0,maximum(d))
    else
        err_msg = "$(type): plot option undefined"
        throw(error(err_msg))
    end
    # plot option
    ms = 2.25
    ms = 0.4*instr[:plot][:dims][1]/meD.nel[1]
    # plotting
    gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,)
    p0 = plot(
        if meD.nD == 2
            mpD.x[:,1],mpD.x[:,2]
        elseif meD.nD == 3
            mpD.x[:,1],mpD.x[:,3]
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

@views function savlot(mpD,meD,t,instr) 
    if !instr[:plast][:status]
        mpD.z₀ .= mpD.x[:,end] 
    end
    if instr[:plot][:status]
        P,tit = [],L"$t = $"*string(round(t,digits=1))*" [s]"
        for (k,what) ∈ enumerate(instr[:plot][:what])
            p0 = whichplot(what,tit,mpD,meD,instr)
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