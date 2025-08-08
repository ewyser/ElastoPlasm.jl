using BenchmarkTools
@views function allocTest(L::Vector{Float64},nel::Int64; kwargs...)
    @warn "Unit testing"
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "Init test geometry"
    # independant physical constant
    if length(L) == 2
        g = vec([0.0,9.81])
    elseif length(L) == 3
        g = vec([0.0,0.0,9.81])
    end                                                           # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmp = setup_cmpr(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mpts setup
    mesh     = setup_mesh(nel,L,instr)                                            # mesh geometry setup
    mpts     = setup_mps(mesh,L,cmp,instr[:grf],typeD)                      # material point geometry setup
    @info "Mesh & mpts feature(s):" instr[:shpfun] instr[:fwrk] instr[:trsfr] instr[:vollock] nel nthreads()
    # plot & time stepping parameters
    tw,Δt,it,ctr,toc,flag,ηmax,ηtot,dt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action
    suite = BenchmarkGroup()
    suite["shpfun"] = BenchmarkGroup(["string", "unicode"])
    suite["solve" ] = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"] = BenchmarkGroup(["string", "unicode"])
    suite["elasto"] = BenchmarkGroup(["string", "unicode"])
    suite["plast" ] = BenchmarkGroup(["string", "unicode"])
    

    suite["shpfun"]["shpfun!"] = @benchmarkable shpfun!($mpts,$mesh,$instr)
    suite["mapsto"]["p->n"   ] = @benchmarkable mapsto!($mpts,$mesh,$g,$dt,$instr,"p->n")
    suite["solve" ]["solve!" ] = @benchmarkable solve!($mesh,$dt)
    #suite["mapsto"]["p<-p"   ] = @benchmarkable mapsto!($mpts,$mesh,$g,$dt,$instr,"p<-n")
    suite["elasto"]["-------"] = @benchmarkable ηmax = elastoplast!($mpts,$mesh,$cmp,$dt,$instr)
    suite["elasto"]["strain!"] = @benchmarkable strain!($mpts,$mesh,$dt,$instr)
    suite["elasto"]["ΔFbar! "] = @benchmarkable ΔFbar!($mpts,$mesh)
    suite["elasto"]["stress!"] = @benchmarkable stress!($mpts,$cmp,$instr,:update)
    
#=
    mpts.Δλ[:] .= 1.0

    ls      = cmp[:nonlocal][:ls]
    mpts.e2p.= Int(0)
    mpts.p2p.= Int(0)
    W,w     = spzeros(mpts.nmp),spzeros(mpts.nmp,mpts.nmp)
    @isdefined(nonloc!) ? nothing : nonloc! = nonlocal(CPU())
    for proc ∈ ["p->q","p<-q"]
        nonloc!(W,w,mpts,mesh,ls,proc; ndrange=mpts.nmp);sync(CPU())
    end

    suite["plast" ]["ϵII0p2q"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mpts,$mesh,$cmp[:nonlocal][:ls],"p->q"; ndrange=$mpts.nmp);sync(CPU())
    suite["plast" ]["ϵII0q2p"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mpts,$mesh,$cmp[:nonlocal][:ls],"p<-q"; ndrange=$mpts.nmp);sync(CPU())
=#
    @info "run benchmarks..."
    benchmark = mean(run(suite))
    for (k,KEY) ∈ enumerate(keys(benchmark))
        for (l,key) ∈ enumerate(keys(benchmark[KEY]))
            benchmark[KEY][key]
        end
    end
    return benchmark
end
export allocTest
#allocCheck([64.1584,64.1584/4,12.80],80)