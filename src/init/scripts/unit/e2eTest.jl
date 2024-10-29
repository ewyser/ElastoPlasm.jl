# L,nel = [64.1584,12.80],40
# e2eTest(L,nel;)
function e2eTest(L::Vector{Float64},nel::Int64; kwargs...)
    configPlot()
    # init & kwargs
    instr  = kwargser(:instr,kwargs;dim=length(L))
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    ni      = 2    
    # constitutive model
    cmParam = cm(length(L),instr)
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)
    setgeom = inislump(meD,cmParam,ni,instr)                       
    mpD     = pointSetup(meD,cmParam,instr;define=setgeom)


    instr[:cairn][:shpfun].tplgy!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())

    ls      = cmParam[:nonlocal][:ls]
    mpD.e2p.= Int(0)
    mpD.p2p.= Int(0)
    mpD.ϵpII[:,2].= 0.0
    W,w     = spzeros(mpD.nmp),spzeros(mpD.nmp,mpD.nmp)
    for proc ∈ ["tplgy"]
        instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mpD,meD,ls,proc; ndrange=mpD.nmp);sync(CPU())
    end





    fSize = (2.0*250,2*125)
    mSize = 0.25*fSize[1]/meD.nel[1]
    gr(size=fSize,legend=true,markersize=2.25,markerstrokecolor=:auto)
    xn = reshape(meD.xn[:,1],meD.nno[2],meD.nno[1])
    yn = reshape(meD.xn[:,2],meD.nno[2],meD.nno[1])
    xe = xn[1:end-1,1:end-1].+0.5*meD.h[1]
    ye = yn[1:end-1,1:end-1].+0.5*meD.h[2]
    for p ∈ 1:mpD.nmp
        ps = findall(!iszero,mpD.e2p[:,mpD.p2e[p]])
        
        plot(xn  ,yn ,seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)
        plot!(xn',yn',seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)

        plot!(xe  ,ye ,seriestype=:scatter,shape=:cross,color=:black,markersize=2.5)

        scatter!(mpD.x[:,1],mpD.x[:,2]  ,c=:black,alpha=0.05,markersize=mSize    ,)
        scatter!(mpD.x[ps,1],mpD.x[ps,2],c=:black,alpha=0.20,markersize=mSize    ,)
        scatter!((mpD.x[p,1],mpD.x[p,2]),c=:green,alpha=1.00,markersize=mSize,markershape=:cross,legend=false,aspect_ratio=1,display=true)
    end
    return msg("(✓) Done! exiting...")
end
export e2eTest
#e2eTest([64.1584,12.80],40;basis="bsmpm")