@testset "+ $(basename(@__FILE__))" verbose = true begin
    function iter_slump(ic,cfg,basis,msg)
        cases = [
            Dict(:deformation => "finite",       :transfer => "musl", :locking => true ),
            Dict(:deformation => "finite",       :transfer => "musl", :locking => false),
            Dict(:deformation => "finite",       :transfer => "tpic", :locking => true ),
            Dict(:deformation => "finite",       :transfer => "tpic", :locking => false),
            Dict(:deformation => "infinitesimal",:transfer => "musl", :locking => true ),
            Dict(:deformation => "infinitesimal",:transfer => "musl", :locking => false),
            Dict(:deformation => "infinitesimal",:transfer => "tpic", :locking => true ),
            Dict(:deformation => "infinitesimal",:transfer => "tpic", :locking => false),
        ]
        prog = Progress(length(cases)+1;dt=0.5,desc=msg,barlen=10);
        for case ∈ cases
            kwargs = Dict(
                :basis  => basis,
                :fwrk   => (; deform = case[:deformation], trsfr = case[:transfer], locking = case[:locking]),
                :nonloc => (; status = false, ls = 0.5),
                :plot   => (; status = true, freq = 1.0, what = ["epII"], dims = (500.0,250.0),),
            )
            instr = kwargser(:instr,kwargs;dim=ic.mesh.dim)
            cfg = (;instr = instr, paths = cfg.paths)

            @testset "$(basename(@__FILE__)) executes with: $(instr[:fwrk].deform), $(instr[:fwrk].trsfr), $(instr[:fwrk].locking)" begin
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

    cases  = [
        (; which = "bsmpm", how = nothing     , ghost = false),
        (; which = "gimpm", how = "undeformed", ghost = true ),
        #(; which = "smpm" , how = nothing     , ghost = true ),
    ]
    for basis ∈ cases
        @info "Testing with $(basis.which) basis"
        @testset "- 2d geometry with $(basis.which) basis" verbose = true begin
            iter_slump(ic,cfg,basis,"Completion for 2d geometry:")
        end
    end

    #=
    L,nel  = [64.1584,64.1584/4.0,64.1584/4.0],[40,10,10];
    ic,cfg = ic_slump(L,nel; fid = "test/slump");
    ic     = merge(ic, (; time = (; T = ic.time.te, te = ic.time.te, tg = ic.time.tg) ))
    @testset "- 3d geometry" verbose = true begin
        iter_slump(ic,cfg,"Completion for 3d geometry:")
    end
    =#
end