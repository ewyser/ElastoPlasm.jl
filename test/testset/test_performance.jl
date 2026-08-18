function define_benchs(ic,cfg)
    suite = BenchmarkGroup()
    suite["ignite"]      = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"]      = BenchmarkGroup(["string", "unicode"])
    suite["elastoplast"] = BenchmarkGroup(["string", "unicode"])
    # unpack mesh, mpts, basis, cmpr, instr, paths as aliases
    mesh,mpts,basis,cmpr = ic["mesh"], ic["mpts"], ic["basis"], ic["cmpr"]
    instr          = cfg["solver"]
    time           = ic["time"]
    g, dt, η, C_pf = [0.0, -9.81], 1e-3, 0.1, 0.99
    # calculate/update topology
    suite["ignite"]["tplgy!"] = @benchmarkable begin
       $instr.cairn.ignite.tplgy!($mpts,$mesh,$basis; ndrange=($mpts.nmp));sync(CPU())
    end
    # cache shape values/gradients for all NN neighbors of every particle
    instr.cairn.ignite.tplgy!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
    suite["ignite"]["shpfun!"] = @benchmarkable begin
       $instr.cairn.ignite.shpfun!($mpts,$mesh,$basis; ndrange=($mpts.nmp));sync(CPU())
    end
    # map material point to node
    suite["mapsto"]["p2n!"  ] = @benchmarkable begin
       $instr.cairn.mapsto.map.p2n!($mpts,$mesh,$basis,$g; ndrange=$mpts.nmp);sync(CPU())
    end
    # solve Eulerian momentum equation
    suite["mapsto"]["solve!"] = @benchmarkable begin
        $instr.cairn.mapsto.map.solve!($mesh,$dt,$η; ndrange=$mesh.prprt.nno[end]);sync(CPU())
    end
    # map back solution to material point
    suite["mapsto"]["n2p!"  ] = @benchmarkable begin
        $instr.cairn.mapsto.map.n2p!($mpts,$mesh,$basis,$dt,$C_pf; ndrange=$mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    suite["mapsto"]["augm"] = @benchmarkable begin
        # accumulate material point contributions
        $instr.cairn.mapsto.augm.p2n!($mpts,$mesh,$basis; ndrange=$mpts.nmp);sync(CPU())
        # solve for nodal incremental displacement
        $instr.cairn.mapsto.augm.solve!($mesh; ndrange=$mesh.prprt.nno[end]);sync(CPU())
    end
    # get incremental deformation tensor
    suite["elastoplast"]["deform!"] = @benchmarkable begin
        $instr.cairn.update.deform!($mpts,$mesh,$basis,$dt; ndrange=$mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    suite["elastoplast"]["ΔJn"] = @benchmarkable begin
        # mapping to mesh
        $instr.cairn.update.ΔJn!($mpts,$mesh,$basis; ndrange=$mpts.nmp);sync(CPU())
    end
    suite["elastoplast"]["ΔJp!"] = @benchmarkable begin
        # compute determinant Jbar
        $instr.cairn.update.ΔJp!($mpts,$mesh,$basis,1.0/(length($mesh.prprt.nel)-1); ndrange=$mpts.nmp);sync(CPU())
    end
    # elastic predictor
    suite["elastoplast"]["elast"] = @benchmarkable begin
        $instr.cairn.update.elast!($mpts,$cmpr.Kc,$cmpr.Gc; ndrange=$mpts.nmp);sync(CPU())
    end
    return suite
end

function run_bench(L,nel)
    simfile = ic_slump(L, nel; fid = "test/performance")
    ic = load(simfile, "ic")
    cfg = load(simfile, "cfg")
    suite = define_benchs(ic, cfg)
    if length(L) == 2
        dim = "2d" 
    elseif length(L) == 3
        dim = "3d"
    end

    @info "Run $dim benchmarks..."
    current = Dict()
    benchmark = mean(run(suite))
    for (group, results) ∈ benchmark
        #println("Group: kernel(s) in $group")
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

    L,nel   = [64.1584,64.1584/4.0],[40,10];
    current["2d"] = run_bench(L,nel)

    #L,nel   = [64.1584,64.1584/4.0,64.1584/4.0],[40,10,10];
    #current["3d"] = run_bench(L,nel)

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
        pct(cur, base) = base == 0 ? 0.0 : round((cur - base) / base * 100, digits=1)
        evidence = String[]
        tot_cur_time,   tot_base_time   = 0.0, 0.0
        tot_cur_memory, tot_base_memory = 0, 0
        tot_cur_allocs, tot_base_allocs = 0, 0
        for dim ∈ ("2d", "3d")
            for (group, res) ∈ baseline[cpu_name][dim]
                @testset "- $dim/$group" verbose = true begin
                    if !haskey(current[dim], group)
                        #@test_broken "Current performance results do not contain group $group"
                    else
                        cur = current[dim][group]
                        push!(evidence,"""
                        $dim/$group:
                          time   : $(cur.mean_time) ms ($(pct(cur.mean_time,res.mean_time))%)
                          memory : $(cur.memory) B ($(pct(cur.memory,res.memory))%)
                          allocs : $(cur.allocs) ($(pct(cur.allocs,res.allocs))%)                          
                        """)
                        tot_cur_time    += cur.mean_time;  tot_base_time   += res.mean_time
                        tot_cur_memory  += cur.memory;     tot_base_memory += res.memory
                        tot_cur_allocs  += cur.allocs;     tot_base_allocs += res.allocs
                    end
                end
            end
        end
        push!(evidence,"""
        overall:
          time   : $(round(tot_cur_time,digits=2)) ms ($(pct(tot_cur_time,tot_base_time))%)
          memory : $(tot_cur_memory) B ($(pct(tot_cur_memory,tot_base_memory))%)
          allocs : $(tot_cur_allocs) ($(pct(tot_cur_allocs,tot_base_allocs))%)
        """)

        # Offer to update the baseline when the *overall* totals (summed across all groups/dims)
        # for time, memory and allocs are all improved - never automatic, always opt-in.
        improved = tot_cur_time ≤ tot_base_time && tot_cur_memory ≤ tot_base_memory && tot_cur_allocs ≤ tot_base_allocs
        if improved && get(ENV, "GITHUB_ACTIONS", "false") != "true"
            println("Current overall performance (summed across all groups) is faster than the baseline:")
            foreach(print, evidence)
            try
                menu   = RadioMenu(["No", "Yes"], pagesize=2)
                choice = request("Update baseline with current performance results?", menu)
                if choice == 2
                    baseline[cpu_name] = current
                    save(path, baseline)
                    @info "Baseline updated with current performance results."
                end
            catch
                nothing
            end
        end
    end
end
