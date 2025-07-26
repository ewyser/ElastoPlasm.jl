#= =#
@views @kernel inbounds = true function nonlocal(W,w,mp,mesh,ls,type)
    p = @index(Global)

    if type == "tplgy" && p ≤ mp.nmp
        for el ∈ findall(!iszero,mesh.e2e[:,mp.p2e[p]])
            mp.e2p[p,el] = p       
        end
    elseif type == "p->q" && p ≤ mp.nmp && mp.Δλ[p] != 0.0
        for (it,q) ∈ enumerate(findall(!iszero,mp.e2p[:,mp.p2e[p]]))
            ξ,η = (mp.x[1,p]-mp.x[1,q]),(mp.x[2,p]-mp.x[2,q])
            d   = sqrt(ξ^2+η^2)
            if w[p,q] == 0.0
                ω₀     = d/ls*exp(-(d/ls)^2)
                w[p,q] = ω₀
                w[q,p] = ω₀
                W[p]  += ω₀
                W[q]  += ω₀
                mp.p2p[q,p] = q
            end
        end
    elseif type == "p<-q" && p ≤ mp.nmp && mp.Δλ[p] != 0.0
        if isapprox(W[p]>1e-16,0.0,atol=1e-16)
            mp.ϵpII[2,p] = mp.ϵpII[1,p]
        else
            for (k,q) ∈ enumerate(findall(!iszero,mp.p2p[:,p]))
                mp.ϵpII[2,p]+= (w[p,q]/W[p])*mp.ϵpII[1,q]
            end
        end
    end
end
 
#= 
@views @kernel inbounds = true function nonlocal(W,w,mp,mesh,ls,type)
    p = @index(Global)

    if type == "p->q" && p ≤ mp.nmp
        for q ∈ mp.nmp
            ξ = (mp.x[p,1]-mp.x[q,1])
            η = (mp.x[p,2]-mp.x[q,2])
            d = sqrt(ξ^2+η^2)
        
            ω₀     = d/ls*exp(-(d/ls)^2)
            w[p,q] = ω₀
            w[q,p] = ω₀
            W[p]  += ω₀
            W[q]  += ω₀
      
        end
    elseif type == "p<-q" && p ≤ mp.nmp
        mp.ϵpII[p,2] = mp.ϵpII[p,1]
        for q ∈ mp.nmp
            if W[p] != 0.0
                mp.ϵpII[p,2]+= (w[p,q]/W[p])*mp.ϵpII[q,1]
            end
        end
    end
end
=#