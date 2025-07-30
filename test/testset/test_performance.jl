using BenchmarkTools, KernelAbstractions
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync

@testset "+ $(basename(@__FILE__))" verbose = true begin
   baseline = Dict(
        "shpfun" => (; mean_time=0.49, memory=missing, allocs=missing),
    )
    suite = BenchmarkGroup()
    suite["shpfun"] = BenchmarkGroup(["string", "unicode"])

    L,nel  = [64.1584,64.1584/4.0],[80,20];
    ic,cfg = ic_slump(L,nel; fid = "test/performance");
    # unpack mesh, mp, cmpr, instr, paths as aliases
    mesh,mp,cmpr = ic[:mesh]  , ic[:mp]    , (ic[:cmpr])
    instr,paths  = cfg[:instr], cfg[:paths]
    time         = ic[:time]

    
    suite["shpfun"]["tplgy!"] = @benchmarkable shpfun!($mp,$mesh,$instr)
    suite["shpfun"]["ϕ∂ϕ!"  ] = @benchmarkable shpfun!($mp,$mesh,$instr)
    # calculate/update topology
    suite["shpfun"]["tplgy!"] = @benchmarkable begin
       $cfg.instr[:cairn][:shpfun].tplgy!($mp,$mesh; ndrange=($mp.nmp));sync(CPU()) 
    end
    # calculate shape functions
    suite["shpfun"]["ϕ∂ϕ!"  ] = @benchmarkable begin
       $cfg.instr[:cairn][:shpfun].ϕ∂ϕ!($mp,$mesh; ndrange=($mp.nmp));sync(CPU()) 
    end

    @info "Run benchmarks..."
    benchmark = mean(run(suite))
    results   = Dict()
    # Print readable summary for each benchmark with statistics


    @info "Benchmark summary:"
    for (group, results) in benchmark
        println("Group: $group")
        tot_mean,tot_memory,tot_alloc = 0.0, 0, 0
        for (name, res) in results
            mean_time = round(res.time / 1e6,digits=2)
            memory    = res.memory
            allocs    = res.allocs
            println("  $name:")
            println("    mean time: $(mean_time) ms")
            println("    memory   : $(memory) bytes"                 )
            println("    allocs   : $(allocs)"                       )
            tot_mean  = tot_mean + mean_time
            tot_alloc = tot_alloc + allocs
            tot_memory = tot_memory + memory
        end
        println("  +-------------------------------------------+")
        println("    tot time  : $(tot_mean) ms"                 )
        println("    tot allocs: $(tot_alloc)"                   )
        println("    tot memory: $(tot_memory) bytes"            )
        # Store results in a dictionary for later use
        results[group] = (; mean_time=tot_mean, memory=tot_memory, allocs=tot_alloc)
        
        @testset "$(group): total execution time" begin 
            @test_broken results[group].mean_time ≈ baseline[group].mean_time 
        end
    end
end










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
    # mesh & mp setup
    mesh     = setup_mesh(nel,L,instr)                                            # mesh geometry setup
    mp     = setup_mps(mesh,L,cmp,instr[:grf],typeD)                      # material point geometry setup
    @info "Mesh & mp feature(s):" instr[:shpfun] instr[:fwrk] instr[:trsfr] instr[:vollock] nel nthreads()
    # plot & time stepping parameters
    tw,Δt,it,ctr,toc,flag,ηmax,ηtot,dt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action
    suite = BenchmarkGroup()
    suite["shpfun"] = BenchmarkGroup(["string", "unicode"])
    suite["solve" ] = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"] = BenchmarkGroup(["string", "unicode"])
    suite["elasto"] = BenchmarkGroup(["string", "unicode"])
    suite["plast" ] = BenchmarkGroup(["string", "unicode"])
    

    suite["shpfun"]["shpfun!"] = @benchmarkable shpfun!($mp,$mesh,$instr)
    suite["mapsto"]["p->n"   ] = @benchmarkable mapsto!($mp,$mesh,$g,$dt,$instr,"p->n")
    suite["solve" ]["solve!" ] = @benchmarkable solve!($mesh,$dt)
    #suite["mapsto"]["p<-p"   ] = @benchmarkable mapsto!($mp,$mesh,$g,$dt,$instr,"p<-n")
    suite["elasto"]["-------"] = @benchmarkable ηmax = elastoplast!($mp,$mesh,$cmp,$dt,$instr)
    suite["elasto"]["strain!"] = @benchmarkable strain!($mp,$mesh,$dt,$instr)
    suite["elasto"]["ΔFbar! "] = @benchmarkable ΔFbar!($mp,$mesh)
    suite["elasto"]["stress!"] = @benchmarkable stress!($mp,$cmp,$instr,:update)
    
#=
    mp.Δλ[:] .= 1.0

    ls      = cmp[:nonlocal][:ls]
    mp.e2p.= Int(0)
    mp.p2p.= Int(0)
    W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
    @isdefined(nonloc!) ? nothing : nonloc! = nonlocal(CPU())
    for proc ∈ ["p->q","p<-q"]
        nonloc!(W,w,mp,mesh,ls,proc; ndrange=mp.nmp);sync(CPU())
    end

    suite["plast" ]["ϵII0p2q"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mp,$mesh,$cmp[:nonlocal][:ls],"p->q"; ndrange=$mp.nmp);sync(CPU())
    suite["plast" ]["ϵII0q2p"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mp,$mesh,$cmp[:nonlocal][:ls],"p<-q"; ndrange=$mp.nmp);sync(CPU())
=#
    @info "run benchmarks..."
    benchmark = mean(run(suite))
    display(benchmark)
    # Print readable summary for each benchmark with statistics
    for (group, results) in benchmark
        println("Group: $group")
        for (name, res) in results
            println("  $name:")
            println("    mean time:   $(res.time) ns")
            println("    memory:     $(res.memory) bytes")
            println("    allocs:     $(res.allocs)")
            println("    samples:    $(res.samples)")
        end
    end
    return benchmark
end
export allocTest
#allocCheck([64.1584,64.1584/4,12.80],80)