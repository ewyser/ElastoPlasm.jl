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
    mesh     = meshSetup(nel,L,instr)
    setgeom = inislump(mesh,cmParam,ni,instr)                       
    mp     = pointSetup(mesh,cmParam,instr;define=setgeom)


    instr[:cairn][:shpfun].tplgy!(mp,mesh; ndrange=(mp.nmp));sync(CPU())

    ls      = cmParam[:nonlocal][:ls]
    mp.e2p.= Int(0)
    mp.p2p.= Int(0)
    mp.ϵpII[:,2].= 0.0
    W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
    for proc ∈ ["tplgy"]
        instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mp,mesh,ls,proc; ndrange=mp.nmp);sync(CPU())
    end





    fSize = (2.0*250,2*125)
    mSize = 0.25*fSize[1]/mesh.nel[1]
    gr(size=fSize,legend=true,markersize=2.25,markerstrokecolor=:auto)
    xn = reshape(mesh.xn[:,1],mesh.nno[2],mesh.nno[1])
    yn = reshape(mesh.xn[:,2],mesh.nno[2],mesh.nno[1])
    xe = xn[1:end-1,1:end-1].+0.5*mesh.h[1]
    ye = yn[1:end-1,1:end-1].+0.5*mesh.h[2]
    for p ∈ 1:mp.nmp
        ps = findall(!iszero,mp.e2p[:,mp.p2e[p]])
        
        plot(xn  ,yn ,seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)
        plot!(xn',yn',seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)

        plot!(xe  ,ye ,seriestype=:scatter,shape=:cross,color=:black,markersize=2.5)

        scatter!(mp.x[:,1],mp.x[:,2]  ,c=:black,alpha=0.05,markersize=mSize    ,)
        scatter!(mp.x[ps,1],mp.x[ps,2],c=:black,alpha=0.20,markersize=mSize    ,)
        scatter!((mp.x[p,1],mp.x[p,2]),c=:green,alpha=1.00,markersize=mSize,markershape=:cross,legend=false,aspect_ratio=1,display=true)
    end
    return msg("(✓) Done! exiting...")
end
export e2eTest
#e2eTest([64.1584,12.80],40;basis="bsmpm")