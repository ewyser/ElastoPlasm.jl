# L,nel = [64.1584,12.80],40
# e2eTest(L,nel;)
function e2eTest(L::Vector{Float64},nel::Int64; kwargs...)
    config_plot()
    # init & kwargs
    instr  = kwargser(:instr,kwargs;dim=length(L))
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    ni      = 2    
    # constitutive model
    cmp = setup_cmpr(length(L),instr)
    # mesh & mpts setup
    mesh     = setup_mesh(nel,L,instr)
    setgeom = inislump(mesh,cmp,ni,instr)                       
    mpts     = setup_mps(mesh,cmp,instr;define=setgeom)


    instr[:cairn][:shpfun].tplgy!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())

    ls      = cmp[:nonlocal][:ls]
    mpts.e2p.= Int(0)
    mpts.p2p.= Int(0)
    mpts.ϵpII[:,2].= 0.0
    W,w     = spzeros(mpts.nmp),spzeros(mpts.nmp,mpts.nmp)
    for proc ∈ ["tplgy"]
        instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mpts,mesh,ls,proc; ndrange=mpts.nmp);sync(CPU())
    end





    fSize = (2.0*250,2*125)
    mSize = 0.25*fSize[1]/mesh.nel[1]
    gr(size=fSize,legend=true,markersize=2.25,markerstrokecolor=:auto)
    xn = reshape(mesh.x[:,1],mesh.nno[2],mesh.nno[1])
    yn = reshape(mesh.x[:,2],mesh.nno[2],mesh.nno[1])
    xe = xn[1:end-1,1:end-1].+0.5*mesh.h[1]
    ye = yn[1:end-1,1:end-1].+0.5*mesh.h[2]
    for p ∈ 1:mpts.nmp
        ps = findall(!iszero,mpts.e2p[:,mpts.p2e[p]])
        
        plot(xn  ,yn ,seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)
        plot!(xn',yn',seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25,markersize=0.01,)

        plot!(xe  ,ye ,seriestype=:scatter,shape=:cross,color=:black,markersize=2.5)

        scatter!(mpts.x[:,1],mpts.x[:,2]  ,c=:black,alpha=0.05,markersize=mSize    ,)
        scatter!(mpts.x[ps,1],mpts.x[ps,2],c=:black,alpha=0.20,markersize=mSize    ,)
        scatter!((mpts.x[p,1],mpts.x[p,2]),c=:green,alpha=1.00,markersize=mSize,markershape=:cross,legend=false,aspect_ratio=1,display=true)
    end
    return msg("(✓) Done! exiting...")
end
export e2eTest
#e2eTest([64.1584,12.80],40;basis="bsmpm")