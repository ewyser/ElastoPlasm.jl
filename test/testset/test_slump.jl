@testset "+ $(basename(@__FILE__))" verbose = true begin
    function iter_slump(ic,cfg,basis,plot,grf,msg)
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
                :plot   => plot,
                :grf    => grf,
            )
            instr = kwargser(:instr,kwargs;dim=ic.mesh.dim)
            cfg = (;instr = instr, paths = cfg.paths)

            @testset "$(basename(@__FILE__)) executes with: $(instr[:fwrk].deform), $(instr[:fwrk].trsfr), $(instr[:fwrk].locking)" begin
                status = false
                try
                    @suppress begin
                        status = slump(ic,cfg; workflow = "all-in-one").success
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
    
    cases  = [
        (; which = "bsmpm", how = nothing     , ghost = false),
        #(; which = "gimpm", how = "undeformed", ghost = true ),
        #(; which = "smpm" , how = nothing     , ghost = true ),
    ]
    plot = (; status = true, freq = 1.0, what = ["epII"], dims = (500.0,250.0),)
    grf  = (; status = true, covariance = "gaussian", param = (; Iₓ= [2.5,2.5,2.5], Nₕ = 5000, kₘ = 100,),)

    for basis ∈ cases
        @info "Testing with $(basis.which) basis"
        # 2d slump tests
        L,nel  = [64.1584,64.1584/4.0],[40,10];
        ic,cfg = ic_slump(L,nel; fid = "test/slump", grf);
        @testset "- 2d geometry with $(basis.which) basis" verbose = true begin
            iter_slump(ic,cfg,basis,plot,grf,"Completion for 2d geometry:")
        end
        # 3d slump tests
        #=
        L,nel  = [64.1584,64.1584/4.0,64.1584/4.0],[40,10,10];
        ic,cfg = ic_slump(L,nel; fid = "test/slump", grf);
        @testset "- 3d geometry with $(basis.which) basis" verbose = true begin
            iter_slump(ic,cfg,basis,plot,grf,"Completion for 3d geometry:")
        end
        =#
    end
end