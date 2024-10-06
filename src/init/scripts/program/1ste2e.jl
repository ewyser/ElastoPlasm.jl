function teste2e(L::Vector{Float64},nel::Int64; kwargs...)
    configPlot()
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init slump geometry"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup

    twoDtplgy!(mpD,meD)


    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    for p ∈ 1:mpD.nmp
        scatter(mpD.x[:,1],mpD.x[:,2], c=:black, alpha=0.1, legend=false)
        scatter!((mpD.x[p,1],mpD.x[p,2]), c=:green, alpha=1.0, legend=false)
        ps,active = [],meD.e2e[mpD.p2e[p]]
        for k ∈ 1:length(active)
            ps = vcat(ps,findall(x->x==active[k],mpD.p2e))
        end
        scatter!(mpD.x[ps,1],mpD.x[ps,2], c=:red, alpha=0.2, legend=false,aspect_ratio=1,display=true)
    end






    return msg("(✓) Done! exiting...")
end
export teste2e