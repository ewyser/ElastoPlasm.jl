function inicmp(meD,cmParam,instr,ni;ℓ₀=0.0)
    coh0,cohr,phi0= cmParam[:c0],cmParam[:cr],cmParam[:ϕ0]
    if meD.nD == 2
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:ℓ₀-0.5*meD.h[2]/ni
        npx,npz     = length(xL),length(zL)
        xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL )) 
        xp,zp       = vec(xp),vec(zp)
    elseif meD.nD == 3
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        yL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:meD.xB[4]
        zL          = meD.xB[5]+(0.5*meD.h[3]/ni):meD.h[3]/ni:ℓ₀-0.5*meD.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        xp,yp,zp    = vec(xp),vec(yp),vec(zp)
    end
    if meD.nD == 2 
        xp = hcat(xp,zp) 
    elseif meD.nD == 3 
        xp = hcat(xp,yp,zp) 
    end
    id   = shuffle(collect(1:size(xp,1)))
    cohr = ones(size(xp,1)).*coh0
    cohr = ones(size(xp,1)).*cohr
    phi  = ones(size(xp,1)).*phi0
    return ni,size(xp,1),(;xp=xp,coh0=coh0,cohr=cohr,phi=phi,)
end

@views function plotStuff(mpD,t,type,ctr,title)
    xlab,ylab = L"$x-$direction",L"$z-$direction"
    gr(size=(2*250,2*125),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.0,markerstrokecolor=:match,)
    temp = title
    if type == "P"
        p = -mpD.σ[end,:]/1e3
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=p,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$\sigma_{zz}$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(0.0,50.0),
            title=temp,
            show=true,
            )  
    end
    return ctr+=1
end
function compactTest(dim,nel,varPlot,ν,E,ρ0,l0; kwargs...)
    @info "execution of compactTest()"
    if dim == 2
        L       = [l0/nel,l0       ]                                           # domain geometry
    elseif dim == 3
        L       = [l0/nel,l0/nel,l0]                                            # domain geometry
    end
    # independant physical constant
    instr  = kwargser(:instr,kwargs)
    fid    = splitext(basename(@__FILE__))
    paths  = setPaths(first(fid), sys.out;interactive=false)
    # independant physical constant
    g       = 9.81                      
    ni      = 2                         
    # constitutive model
    cm(dim,instr; E=E,ν=ν,ρ0=ρ0)
    cmParam = cm(length(L),instr)
    tg      = ceil((1.0/cmParam.c)*(2.0*l0)*40.0)
    T,te    = 1.25*tg,1.25*tg   
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)    
    setgeom = inicmp(meD,cmParam,instr,ni;ℓ₀=l0) 
    mpD     = pointSetup(meD,cmParam,instr;define=setgeom)                                        
    z0      = maximum(mpD.x[:,end])
    # action
    out = ϵp23De!(mpD,meD,cmParam,g,T,te,tg,instr)    
    # analytics
    if meD.nD==2
        xN,yN = abs.(mpD.σᵢ[2,:]),z0
    elseif meD.nD==3
        xN,yN = abs.(mpD.σᵢ[3,:]),z0
    end
    xA,yA = abs.(cmParam.ρ0.*g.*(l0.-z0)),z0
    err   = sum(sqrt.((xN.-xA).^2).*mpD.Ω₀)/(abs(g[end])*cmParam.ρ0*l0*sum(mpD.Ω₀)) 
    return (xN,yN,xA,yA),meD.h,err 
end
@views function compacTest(dim,trsfrAp)
    ϕ∂ϕType    = "gimpm"
    fwrkDeform = "finite"
    store,H,error = [],[],[]
    try
        @testset "convergence using $(ϕ∂ϕType), $(fwrkDeform) deformation, $(trsfrAp) mapping" begin
            # geometry
            n         = [0,1,2,3]#[0,1,2,3,4,5,6]
            nel       = 2.0.^n
            # initial parameters 
            l0,ν,E,ρ0 = 50.0,0.0,1.0e4,80.0
            # init error
            ϵ         = 1.0
            for (it,nel) in enumerate(nel)
                #action
                DAT,h,err = compactTest(dim,nel,"P",ν,E,ρ0,l0;basis=ϕ∂ϕType,fwrk=fwrkDeform,trsfr=trsfrAp,vollock=false)
                push!(store,DAT )
                push!(H ,h[end])
                push!(error,err)
                # test
                @test (err < ϵ )
                ϵ = err
            end 
        end
    catch

    end
    xN,yN,xA,yA = store[end]
    return H,error
end
#=
if @isdefined(perf)
    if perf
        runCompacTest([2,3],[:mUSL,:mUSL])
    else
        runCompacTest([2,3],[:mUSL,:mUSL])
    end
else
    runCompacTest([2,3,2,3],[:mUSL,:mUSL,:tpicUSL,:tpicUSL])
end
=#
export compacTest

#compacTest(2,"mUSL")



































#=
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    p1 = plot(xN.*1e-3,yN,seriestype=:scatter, label="$(dim)D $(ϕ∂ϕType), $(trsfrAp) mapping")
    p1 = plot!(xA.*1e-3,yA,label=L"\sum_{p}\dfrac{||\sigma_{yy}^p-\sigma_{yy}^a(x_p)||V_0^p}{(g\rho_0l_0)V_0}",xlabel=L"$\sigma_{yy}$ [kPa]",ylabel=L"$y-$position [m]") 
    display(plot(p1; layout=(1,1), size=(450,250)))
    savefig(path_test*"$(dim)D_numericVsAnalytic_compacTest_$(ϕ∂ϕType)_$(fwrkDeform)_$(trsfrAp).png")
=#