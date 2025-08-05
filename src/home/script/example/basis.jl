export basisfunction

function geom_basis(mesh,cmp,instr; ni = 2)
    if mesh.dim == 2
        x          = collect(mesh.xB[1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[2])
        z          = collect(mesh.xB[3]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[4])
        nmp        = [length(x),length(z),length(x)*length(z)]
        xp         = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp         = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])
    elseif mesh.dim == 3
        x          = collect(mesh.xB[1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[2])
        y          = collect(mesh.xB[3]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[4])
        z          = collect(mesh.xB[5]+(0.5*mesh.h[3]/ni):mesh.h[3]/ni:mesh.xB[6])
        nmp        = [length(x),length(y),length(z),length(x)*length(y)*length(z)]
        xp         = repeat(reshape(x,1     ,nmp[1],1     ),nmp[3],1     ,nmp[2])
        yp         = repeat(reshape(y,1     ,1     ,nmp[2]),nmp[3],nmp[1],1     )
        zp         = repeat(reshape(z,nmp[3],1     ,1     ),1     ,nmp[1],nmp[2])
    end
    if mesh.dim == 2 
        xp = vcat(vec(xp)',vec(zp)') 
    elseif mesh.dim == 3 
        xp = vcat(vec(xp)',vec(yp)',vec(zp)') 
    end
    nmp  = size(xp,2)
    coh0 = ones(nmp).*cmp[:c0]
    cohr = ones(nmp).*cmp[:cr]
    phi  = ones(nmp).*cmp[:ϕ0]
    return ni,nmp,(;xp=xp,coh0=coh0,cohr=cohr,phi=phi,)
end

function basisfunction(; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Calculating shape functions for basis functions"
    nel,L = [4,4],[1.0,1.0]
    ni    = 20
    # init & kwargs
    instr = kwargser(:instr,kwargs;dim=length(L))
    paths = set_paths(fid,info.sys.out;interactive=false)  
    T0    = instr[:dtype].T0  
    T1,T2 = first(T0),last(T0)  
    # mesh & mp initial conditions
    mesh  = setup_mesh(nel,L,instr)
    cmpr  = setup_cmpr(mesh,instr)                       
    mp    = setup_mps(mesh,cmpr;define=geom_basis(mesh,cmpr,instr; ni = ni))
    # calculate tplgy and shpfun
    shpfun(mp,mesh,instr)

    # extract and store value of mp.ϕ∂ϕ
    N   = spzeros(T2,mesh.nno[end],mp.nmp)
    dNx = spzeros(T2,mesh.nno[end],mp.nmp)
    dNy = spzeros(T2,mesh.nno[end],mp.nmp)
    X   = spzeros(T2,mesh.nno[end],mp.nmp)
    Y   = spzeros(T2,mesh.nno[end],mp.nmp)
    for p = 1:mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if iszero(no) continue end
            N[no,p]   = mp.ϕ∂ϕ[nn,p,1]
            dNx[no,p] = mp.ϕ∂ϕ[nn,p,2]
            dNy[no,p] = mp.ϕ∂ϕ[nn,p,3]
            X[no,p]   = mp.x[1,p]
            Y[no,p]   = mp.x[2,p]
        end
    end

    x_coord = reshape(mesh.x[1, :], Int(mesh.nno[2]), Int(mesh.nno[1]))
    y_coord = reshape(mesh.x[2, :], Int(mesh.nno[2]), Int(mesh.nno[1]))

    config_plot()
    gr(size=(250,250),legend=false,markersize=2.0,markerstrokecolor=:auto)
    for no = 1:mesh.nno[end]
        p0 = plot( x_coord ,y_coord ,zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")
        p0 = plot!(x_coord',y_coord',zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")

        x = X[no,:]
        y = Y[no,:]
        n = N[no,:]

        x[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        y[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        n[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization

        p0 = plot!(
            x,y,n,
            seriestype = :scatter,
            xlims = (0.0, 1.0),
            ylims = (0.0, 1.0),
            zlims = (0.0, 1.0),
            framestyle = :box,
            grid = true,
            )
        display(plot(p0;layout=(1,1))) 
        path = mkpath(joinpath(paths[:plot],"$(instr[:basis][:which])","N_$(T2)"))
        savefig(joinpath(path,"N_node_$(no).png"))
    end
    for no = 1:mesh.nno[end]
        p0 = plot( x_coord ,y_coord ,zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")
        p0 = plot!(x_coord',y_coord',zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")

        x = X[no,:]
        y = Y[no,:]
        n = dNx[no,:]

        x[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        y[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        n[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization

        p0 = plot!(
            x,y,n,
            seriestype = :scatter,
            xlims = (0.0, 1.0),
            ylims = (0.0, 1.0),
            zlims = (-2.0, 2.0),
            framestyle = :box,
            grid = true,
            )
        display(plot(p0;layout=(1,1))) 
        path = mkpath(joinpath(paths[:plot],"$(instr[:basis][:which])","dNx_$(T2)"))
        savefig(joinpath(path,"dNx_node_$(no).png"))
    end
    for no = 1:mesh.nno[end]
        p0 = plot( x_coord ,y_coord ,zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")
        p0 = plot!(x_coord',y_coord',zeros(size(x_coord)),seriestype = :line, color = :black, label = "mesh nodes")

        x = X[no,:]
        y = Y[no,:]
        n = dNy[no,:]

        x[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        y[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization
        n[iszero.(n)] .= NaN  # replace zeros with NaN for better visualization

        p0 = plot!(
            x,y,n,
            seriestype = :scatter,
            xlims = (0.0, 1.0),
            ylims = (0.0, 1.0),
            zlims = (-2.0, 2.0),
            framestyle = :box,
            grid = true,
            )
        display(plot(p0;layout=(1,1))) 
        path = mkpath(joinpath(paths[:plot],"$(instr[:basis][:which])","dNy_$(T2)"))
        savefig(joinpath(path,"dNy_node_$(no).png"))
    end
    return nothing
end