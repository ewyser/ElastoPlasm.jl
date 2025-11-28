@testset "+ $(basename(@__FILE__))" verbose = true begin
    function iter_slump(sim,basis,plot,grf,msg)
        file = jldopen(sim,"r")        
        mesh,mpts,cmpr   = file["ic"]["mesh"]  , file["ic"]["mpts"]  , (file["ic"]["cmpr"])
        instr,paths,misc = file["cfg"]["instr"], file["cfg"]["paths"], file["cfg"]["misc"]
        time             = file["ic"]["time"]
        close(file)

        cases = [
            Dict(:deformation => "finite",       :transfer => "std" , :locking => true ),
            Dict(:deformation => "finite",       :transfer => "std" , :locking => false),
            Dict(:deformation => "finite",       :transfer => "tpic", :locking => true ),
            Dict(:deformation => "finite",       :transfer => "tpic", :locking => false),
            Dict(:deformation => "finite",       :transfer => "apic", :locking => true ),
            Dict(:deformation => "finite",       :transfer => "apic", :locking => false),
            Dict(:deformation => "infinitesimal",:transfer => "std" , :locking => true ),
            Dict(:deformation => "infinitesimal",:transfer => "std" , :locking => false),
            Dict(:deformation => "infinitesimal",:transfer => "tpic", :locking => true ),
            Dict(:deformation => "infinitesimal",:transfer => "tpic", :locking => false),
            Dict(:deformation => "infinitesimal",:transfer => "apic", :locking => true ),
            Dict(:deformation => "infinitesimal",:transfer => "apic", :locking => false),
        ]
        prog = Progress(length(cases)+1;dt=0.5,desc=msg,barlen=10);
        for (k,case) ∈ enumerate(cases)
            kwargs = Dict(
                :basis  => basis,
                :fwrk   => (; deform = case[:deformation], trsfr = case[:transfer], musl = true, C_pf = 0.99, locking = case[:locking], damping = 0.1,),
                :nonloc => (; status = false, ls = 0.5),
                :plot   => plot,
                :grf    => grf,
                :backend=> instr.backend,
            )
            instr = kwargser(kwargs, Instruction; dim=mesh.prprt.dim)
            misc = (;
                prefix = "case_$(k)_$(mesh.prprt.dim)d_$(instr.fwrk.trsfr)"
            )

            sim = export_setup(mesh,mpts,cmpr,time,instr,paths,misc; path = paths[:dat], file = "slump_simulation")

            @testset "$(basename(@__FILE__)) executes with: $(instr.fwrk.deform), $(instr.fwrk.trsfr), $(instr.fwrk.locking)" begin
                status = false
                try
                    @suppress begin
                        status = elastoplasm(sim; workflow = [elastodynamic!,elastoplastic!]).success
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
    plot = (; status = true, dpi = 500, freq = 1.0, what = [("mpts","epII")], dims = (500.0,250.0),cblim  = [(0.0, 1.5),],)
    grf  = (; status = true, covariance = "gaussian", param = (; Iₓ= [2.5,2.5,2.5], Nₕ = 5000, kₘ = 100,),)
    for basis ∈ cases
        @info "Testing with $(basis.which) basis"
        # 2d slump tests
        L,nel  = [64.1584,64.1584/4.0],[40,10];
        sim    = ic_slump(L,nel; fid = "test/slump", grf);
        @testset "- 2d geometry with $(basis.which) basis" verbose = true begin
            iter_slump(sim,basis,plot,grf,"Completion for 2d geometry:")
        end
        # 3d slump tests
        L,nel  = [64.1584,64.1584/4.0,64.1584/4.0],[40,10,10];
        sim    = ic_slump(L,nel; fid = "test/slump", grf);
        @testset "- 3d geometry with $(basis.which) basis" verbose = true begin
            iter_slump(sim,basis,plot,grf,"Completion for 3d geometry:")
        end
    end
end