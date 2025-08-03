function define_benchs(ic,cfg)
    suite = BenchmarkGroup()
    suite["shpfun"]      = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"]      = BenchmarkGroup(["string", "unicode"])
    suite["elastoplast"] = BenchmarkGroup(["string", "unicode"])
    # unpack mesh, mp, cmpr, instr, paths as aliases
    mesh,mp,cmpr = ic[:mesh]  , ic[:mp]    , (ic[:cmpr])
    instr,paths  = cfg[:instr], cfg[:paths]
    time         = ic[:time]
    g, dt, η     = [0.0, -9.81], 1e-3, 0.1
    # calculate/update topology
    suite["shpfun"]["tplgy!"] = @benchmarkable begin
       $cfg.instr[:cairn][:shpfun].tplgy!($mp,$mesh; ndrange=($mp.nmp));sync(CPU()) 
    end
    # calculate shape functions
    suite["shpfun"]["ϕ∂ϕ!"  ] = @benchmarkable begin
       $cfg.instr[:cairn][:shpfun].ϕ∂ϕ!($mp,$mesh; ndrange=($mp.nmp));sync(CPU()) 
    end
    # map material point to node
    suite["mapsto"]["p2n!"  ] = @benchmarkable begin
       $cfg.instr[:cairn][:mapsto][:map].p2n!(ndrange=$mp.nmp,$mp,$mesh,$g);sync(CPU())
    end
    # solve Eulerian momentum equation
    suite["mapsto"]["solve!"] = @benchmarkable begin
        $cfg.instr[:cairn][:mapsto][:map].solve!(ndrange=$mesh.nno[end],$mesh,$dt,$η);sync(CPU())
    end
    # map back solution to material point
    suite["mapsto"]["n2p!"  ] = @benchmarkable begin
        $cfg.instr[:cairn][:mapsto][:map].n2p!(ndrange=$mp.nmp,$mp,$mesh,$dt);sync(CPU())
    end
    # volumetric locking correction
    suite["mapsto"]["augm"] = @benchmarkable begin
        # accumulate material point contributions
        $cfg.instr[:cairn][:mapsto][:augm].p2n!(ndrange=$mp.nmp,$mp,$mesh);sync(CPU())
        # solve for nodal incremental displacement
        $cfg.instr[:cairn][:mapsto][:augm].solve!(ndrange=$mesh.nno[end],$mesh);sync(CPU())
        # update material point's displacement
        $cfg.instr[:cairn][:mapsto][:augm].Δu!(ndrange=$mp.nmp,$mp,$mesh,$dt);sync(CPU())
    end
    # get incremental deformation tensor
    suite["elastoplast"]["deform!"] = @benchmarkable begin
        $cfg.instr[:cairn][:elastoplast][:update].deform!(ndrange=$mp.nmp,$mp,$mesh,$dt);sync(CPU())
    end
    # volumetric locking correction
    suite["elastoplast"]["locking"] = @benchmarkable begin
        # mapping to mesh 
        $cfg.instr[:cairn][:elastoplast][:update].ΔJn!(ndrange=$mp.nmp,$mp,$mesh);sync(CPU())
        # compute nodal determinant of incremental deformation 
        $cfg.instr[:cairn][:elastoplast][:update].ΔJs!(ndrange=$mesh.nno[end],$mesh);sync(CPU())
        # compute determinant Jbar 
        $cfg.instr[:cairn][:elastoplast][:update].ΔJp!(ndrange=$mp.nmp,$mp,$mesh,1/$mesh.dim);sync(CPU())
    end
    # elastic predictor
    suite["elastoplast"]["elast"] = @benchmarkable begin
        $cfg.instr[:cairn][:elastoplast][:elast].elast!(ndrange=$mp.nmp,$mp,$cmpr.Del);sync(CPU())
    end
    return suite
end

function run_bench(L,nel)
    ic,cfg = ic_slump(L,nel; fid = "test/performance");
    suite  = define_benchs(ic,cfg)
    if length(L) == 2
        dim = "2d" 
    elseif length(L) == 3
        dim = "3d"
    end

    @info "Run $dim benchmarks..."
    current = Dict()
    benchmark = mean(run(suite))
    for (group, results) ∈ benchmark
        println("Group: kernel(s) in $group")
        tot_mean,tot_memory,tot_alloc = 0.0, 0, 0
        for (name, res) in results
            mean_time = round(res.time / 1e6,digits=2)
            memory    = res.memory
            allocs    = res.allocs
            println("  $name:")
            println("    mean time: $(mean_time) ms")
            println("    memory   : $(memory) bytes"                 )
            println("    allocs   : $(allocs)"                       )
            tot_mean   = round(tot_mean + mean_time,digits=2)
            tot_alloc  = tot_alloc + allocs
            tot_memory = tot_memory + memory
        end
        # Store results in a dictionary for later use
        current[group] = (; mean_time=tot_mean, memory=tot_memory, allocs=tot_alloc)
        #println("")
    end
    @info "Overall summary:"
    for (key, group) ∈ current
        println("Results for group: $key")
        println("  mean time: $(group.mean_time) ms")
        println("  memory   : $(group.memory) bytes")
        println("  allocs   : $(group.allocs)")
    end
    println("")
    return current
end


@testset "+ $(basename(@__FILE__))" verbose = true begin

    current = Dict(
        "2d" => Dict(),
        "3d" => Dict(),
    )

    L,nel   = [64.1584,64.1584/4.0],[80,20];
    current["2d"] = run_bench(L,nel)



    path = joinpath(DATASET, "performance_baseline.jld2")
    cpu_name = split(string(Sys.cpu_info()[1]), ":")[1]
    if "performance_baseline.jld2" ∉ readdir(DATASET) 
        @info "No baseline found, saving current performance results."
        # Save the current benchmark results as a baseline
        baseline = Dict(cpu_name => current); save(path, baseline)
    else
        @info "Baseline found, comparing current performance against it."
        # Load the baseline for comparison
        baseline = load(path)
        if !haskey(baseline,cpu_name)
            baseline[cpu_name] = current
            save(path, baseline)
        end
        # Compare current results with the baseline
        for (group, res) ∈ baseline[cpu_name]["2d"]
            @testset "- $group" verbose = true begin
                if !haskey(current["2d"], group)
                    @test_broken "Current performance results do not contain group $group"                    
                else
                    println("Testing $(group):")
                    println("current: $(current["2d"][group].mean_time) ms ≤ baseline: $(res.mean_time) ms")
                    @test current["2d"][group].mean_time ≤ res.mean_time * (1 + 0.05)
                end
            end
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