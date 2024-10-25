function slumpCheck(scheme,fwrk,shp)
    L,nel = [64.1584,12.80],40
    if shp == "gimpm"
        instr = slump(L,nel;fwrk=fwrk,trsfr=scheme,basis=(;which=shp,how="undeformed",ghost=true,),plot=(;cond=false,freq=1.0,what=["epII"],dims=(500.0,250.0),),);
    else
        instr = slump(L,nel;fwrk=fwrk,trsfr=scheme,basis=(;which=shp,how=nothing,ghost=false,),plot=(;cond=false,freq=1.0,what=["epII"],dims=(500.0,250.0),),);
    end
    return instr
end
function slumpTest()
    @testset "slump()" verbose = true begin
        for shp ∈ ["bsmpm","gimpm"]#["bsmpm","smpm","gimpm"]
            for fwrk ∈ ["finite","infinitesimal"]
                for scheme ∈ ["mUSL","tpicUSL"]
                    @testset "$(scheme), $(fwrk), $(shp)" verbose = true begin
                        @test isa(slumpCheck(scheme,fwrk,shp),Dict) 
                    end
                end
            end    
        end
    end
    return nothing
end
export slumpTest