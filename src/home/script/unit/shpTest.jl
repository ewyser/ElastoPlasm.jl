function shpfunCheck(shp,instr,paths)
    nel,L  = 2,[1.0]

    mesh,ni = meshSetup(nel,L,instr),10

    xp     = collect(mesh.xB[1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[2])
    nel    = mesh.nel[end]
    nmp    = length(xp)
    # constructor
    mp = (
        ndim   = mesh.dim,
        nmp  = nmp,
        x    = xp,
        ℓ    = ones(nmp).*3.0.*L./nmp,
        ϕ∂ϕ  = zeros(mesh.nn,nmp ,mesh.dim+1   ),
        δnp  = zeros(mesh.nn,mesh.dim,nmp      ),
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,mesh.nn,nmp),
    )
    instr[:cairn] = (;shpfun = init_shpfun(mesh.dim,instr[:basis]),)
    # calculate tplgy and shpfun
    shpfun(mp,mesh,instr)
    # extract and store value of mp.ϕ∂ϕ
    xp,ϕ,∂ϕ = zeros(nmp,nel+1),zeros(nmp,nel+1),zeros(nmp,nel+1)  
    PoU = zeros(Float64,nmp)
    for mp ∈ 1:nmp
        ϕ∂ϕ = mp.ϕ∂ϕ[:,mp,:]
        for (k,nn) ∈ enumerate(mp.p2n[:,mp]) if nn<0 continue end
            xp[mp,nn] = mp.x[mp]  
            ϕ[mp,nn]  = ϕ∂ϕ[k,1]  
            ∂ϕ[mp,nn] = ϕ∂ϕ[k,2]  
            PoU[mp]  +=ϕ∂ϕ[k,1]
        end
    end

    ϕmax  = maximum(abs.(ϕ))
    ϕ[ϕ.<1e-4].=NaN
    ∂ϕmax = maximum(abs.(∂ϕ))
    ∂ϕ[isnan.(ϕ)].=NaN

    colors = [:aqua,:fuchsia,:navy,:indigo,:tomato,:olive]


    configPlot()
    T = [L"\phi_n(x_p)",L"\partial_x\phi_n(x_p)",L"$\sum_n\phi_n(x_p)=1$",L"$\Delta = x_p-\sum_n\phi_n(x_p)x_n$"]
    gr(size=(2.0*250,3*125),legend=false,markersize=2.25,markerstrokecolor=:auto)
    
    
    

    xp0,ϕ0,colors = mesh.xn,zeros(length(mesh.xn)),[]
    p0 = plot(ylim= (0.0,1.1*ϕmax),title= T[1],)
    for i = 1:length(mesh.xn)
        clr = i % length(mesh.xn) + 1
        p0 = plot!(xp[:,i],ϕ[:,i],seriestype = :line   ,color=clr,)
        push!(colors,clr)
    end
    p0 = plot!(xp0,ϕ0,seriestype = :scatter,markershape= :square,markersize = 2*2.25, markercolor=colors,)

    xp0,ϕ0,colors = mesh.xn,zeros(length(mesh.xn)),[]
    p1 = plot(ylim = (-1.1*∂ϕmax,1.1*∂ϕmax) , title = T[2],)
    for i = 1:length(mesh.xn)
        clr = i % length(mesh.xn) + 1
        p1 = plot!(xp[:,i],∂ϕ[:,i],seriestype = :line   ,color=clr,)
        push!(colors,clr)
    end
    p1 = plot!(xp0,ϕ0,seriestype = :scatter,markershape= :square,markersize = 2*2.25, markercolor=colors,)

    p2 = plot(
        mp.x,PoU,
        seriestype = :line,
        color      = :black,
        xlabel     = L"$x-$direction [m]",
        ylim       = (1.0-0.1,1.0+0.1),
        yscale     = :log10,
        title      = T[3],
    )
    p2 = plot!(
        mesh.xn,ones(size(mesh.xn)),
        seriestype = :scatter,
        markershape= :square,
        color      = colors,
        markersize = 2*2.25,
        yscale     = :log10,
        title      = T[3],
    )

    display(plot(p0,p1,p2;layout=(3,1))) 
    savefig(joinpath(paths[:plot],"summary_$(shp)"))
    return PoU
end
function shpTest(;ξ::Real=0.90,ghost::Bool=false)
    fid   = splitext(basename(@__FILE__))
    instr = require(:instr)
    paths  = setPaths(first(fid), sys.out;interactive=false)
    for (k,ξ) ∈ enumerate([0.9])
        @testset "partition of unity (PoU) testset ξ = $(round(ξ,digits=2))" verbose = true begin 
            for shp ∈ ["bsmpm","smpm","gimpm"]
                if shp == "gimpm"
                    instr[:basis] = (;which=shp,how="undeformed",ghost=ghost)
                else
                    instr[:basis] = (;which=shp,how=nothing,ghost=ghost)
                end
                @testset "$(shp): $(round(ξ,digits=2)) < PoU < $(round(1.0+(1.0-ξ),digits=2))" verbose = true begin
                    PoU = shpfunCheck(shp,instr,paths)
                    @test ξ < minimum(PoU) < 1.0+(1.0-ξ)
                end
            end
        end
    end
    return nothing
end
export shpTest