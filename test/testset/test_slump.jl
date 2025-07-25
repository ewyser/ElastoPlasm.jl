cases = [
    Dict(:locking => true,  :transfer => "musl", :deformation => "finite"       , :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => true,  :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => false, :transfer => "musl", :deformation => "finite"       , :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => false, :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "bsmpm", how = nothing, ghost = false)),
    #=
    Dict(:locking => true,  :transfer => "musl", :deformation => "finite"       , :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "musl", :deformation => "finite"       , :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "gimpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "gimpm", how = "undeformed", ghost = true)),

    Dict(:locking => true,  :transfer => "musl", :deformation => "finite"       , :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => true,  :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "musl", :deformation => "finite"       , :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "musl", :deformation => "infinitesimal", :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "finite"       , :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    Dict(:locking => false, :transfer => "tpic", :deformation => "infinitesimal", :basis => (; which = "smpm", how = "undeformed", ghost = true)),
    =#
]

@testset "+ $(basename(@__FILE__))" verbose = true begin
    function iter_slump(ic,cfg,msg)
        
        prog = Progress(length(cases)+1;dt=0.5,desc=msg,barlen=10);
        for case ∈ cases
            kwargs = Dict(
                :basis  => case[:basis],
                :fwrk   => (; deform = case[:deformation], trsfr = case[:transfer], locking = case[:locking]),
                :nonloc => (; status = false, ls = 0.5),
                :plot   => (; status = true, freq = 1.0, what = ["epII"], dims = (500.0,250.0),),
            )
            instr = kwargser(:instr,kwargs;dim=ic.mesh.dim)
            cfg = (;instr = instr, paths = cfg.paths)

            @testset "$(basename(@__FILE__)) executes safely with: $(instr[:fwrk].locking), $(instr[:fwrk].trsfr), $(instr[:fwrk].deform), $(instr[:basis].which)" begin
                status = false
                try
                    @suppress begin
                        status = slump(ic,cfg;)
                    end
                catch e
                    @warn "$(basename(@__FILE__)) failed with error" exception = e
                end
                @test status
            end
            next!(prog)
        end
        finish!(prog)   
    end
    
    @info "Launching ./$(basename(@__FILE__))"
    L,nel  = [64.1584,64.1584/4.0],[40,10];
    ic,cfg = ic_slump(L,nel; fid = "test/slump");
    @testset "- 2d geometry" verbose = true begin
        #iter_slump(ic,cfg,"Completion for 2d geometry:")
    end
    L,nel  = [64.1584,64.1584/4.0,64.1584/4.0],[40,10,10];
    ic,cfg = ic_slump(L,nel; fid = "test/slump");
    ic     = merge(ic, (; time = (; T = ic.time.te, te = ic.time.te, tg = ic.time.tg) ))
    @testset "- 3d geometry" verbose = true begin
        #iter_slump(ic,cfg,"Completion for 3d geometry:")
    end
end