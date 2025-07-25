push!(LOAD_PATH, "../src")

using Test,ProgressMeter,Suppressor,ElastoPlasm

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
    function iter_slump(L,nel,msg)
        count,n = 0,length(cases)
        prog = Progress(n;dt=0.5,desc=msg,barlen=10)
        for case ∈ cases
            locking = case[:locking]
            transfer = case[:transfer]
            deformation = case[:deformation]

            basis  = case[:basis]
            fwrk   = (; deform = deformation, trsfr = transfer, locking = locking)
            nonloc = (; status = false, ls = 0.5)
            plot   = (;status=true,freq=1.0,what=["epII"],dims=(500.0,250.0),)

            @testset "$(basename(@__FILE__)) executes safely with: $(locking), $(transfer), $(deformation), $(basis.which)" begin
                try
                    @suppress begin
                        slump(L,nel; basis,fwrk,nonloc,plot,fid = "test/slump")
                    end
                    @test true
                catch e
                    @warn "$(basename(@__FILE__)) failed with error" exception = e
                    @test false
                end
            end
            count+=1
            next!(prog)
        end
        finish!(prog)   
    end
    
    @info "Testing $(basename(@__FILE__)):"
    @testset "2D geometry" verbose = true begin
        @info "Considering 2D geometry"
         L,nel = [64.1584,64.1584],[40,40]
        iter_slump(L,nel,"Completion:")
    end
    @testset "3D geometry" verbose = true begin
        @info "Considering 3D geometry"
        L,nel = [64.1584,64.1584,12.80],20
        #iter_slump(L,nel,"Completion:")
    end
end