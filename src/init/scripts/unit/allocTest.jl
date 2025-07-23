using BenchmarkTools
@views function allocTest(L::Vector{Float64},nel::Int64; kwargs...)
    @warn "unit testing"
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init test geometry"
    # independant physical constant
    if length(L) == 2
        g = vec([0.0,9.81])
    elseif length(L) == 3
        g = vec([0.0,0.0,9.81])
    end                                                           # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    mesh     = meshSetup(nel,L,instr)                                            # mesh geometry setup
    mp     = pointSetup(mesh,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup
    @info "mesh & mp feature(s):" instr[:shpfun] instr[:fwrk] instr[:trsfr] instr[:vollock] nel nthreads()
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot,Δt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action
    suite = BenchmarkGroup()
    suite["shpfun"] = BenchmarkGroup(["string", "unicode"])
    suite["solve" ] = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"] = BenchmarkGroup(["string", "unicode"])
    suite["elasto"] = BenchmarkGroup(["string", "unicode"])
    suite["plast" ] = BenchmarkGroup(["string", "unicode"])
    

    suite["shpfun"]["shpfun!"] = @benchmarkable shpfun!($mp,$mesh,$instr)
    suite["mapsto"]["p->n"   ] = @benchmarkable mapsto!($mp,$mesh,$g,$Δt,$instr,"p->n")
    suite["solve" ]["solve!" ] = @benchmarkable solve!($mesh,$Δt)
    #suite["mapsto"]["p<-p"   ] = @benchmarkable mapsto!($mp,$mesh,$g,$Δt,$instr,"p<-n")
    suite["elasto"]["-------"] = @benchmarkable ηmax = elastoplast!($mp,$mesh,$cmParam,$Δt,$instr)
    suite["elasto"]["strain!"] = @benchmarkable strain!($mp,$mesh,$Δt,$instr)
    suite["elasto"]["ΔFbar! "] = @benchmarkable ΔFbar!($mp,$mesh)
    suite["elasto"]["stress!"] = @benchmarkable stress!($mp,$cmParam,$instr,:update)
    
#=
    mp.Δλ[:] .= 1.0

    ls      = cmParam[:nonlocal][:ls]
    mp.e2p.= Int(0)
    mp.p2p.= Int(0)
    W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
    @isdefined(nonloc!) ? nothing : nonloc! = nonlocal(CPU())
    for proc ∈ ["p->q","p<-q"]
        nonloc!(W,w,mp,mesh,ls,proc; ndrange=mp.nmp);sync(CPU())
    end

    suite["plast" ]["ϵII0p2q"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mp,$mesh,$cmParam[:nonlocal][:ls],"p->q"; ndrange=$mp.nmp);sync(CPU())
    suite["plast" ]["ϵII0q2p"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mp,$mesh,$cmParam[:nonlocal][:ls],"p<-q"; ndrange=$mp.nmp);sync(CPU())
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