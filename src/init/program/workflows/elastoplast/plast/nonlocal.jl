#= =#
@views @kernel inbounds = true function nonlocal(W,w,mpD,meD,ls,type)
    p = @index(Global)

    if type == "tplgy" && p ≤ mpD.nmp
        for el ∈ findall(!iszero,meD.e2e[:,mpD.p2e[p]])
            mpD.e2p[p,el] = p       
        end
    elseif type == "p->q" && p ≤ mpD.nmp && mpD.Δλ[p] != 0.0
        for (it,q) ∈ enumerate(findall(!iszero,mpD.e2p[:,mpD.p2e[p]]))
            ξ,η = (mpD.x[p,1]-mpD.x[q,1]),(mpD.x[p,2]-mpD.x[q,2])
            d   = sqrt(ξ^2+η^2)
            if w[p,q] == 0.0
                ω₀     = d/ls*exp(-(d/ls)^2)
                w[p,q] = ω₀
                w[q,p] = ω₀
                W[p]  += ω₀
                W[q]  += ω₀
                mpD.p2p[q,p] = q
            end
        end
    elseif type == "p<-q" && p ≤ mpD.nmp && mpD.Δλ[p] != 0.0
        if isapprox(W[p]>1e-16,0.0,atol=1e-16)
            mpD.ϵpII[p,2] = mpD.ϵpII[p,1]
        else
            for (k,q) ∈ enumerate(findall(!iszero,mpD.p2p[:,p]))
                mpD.ϵpII[p,2]+= (w[p,q]/W[p])*mpD.ϵpII[q,1]
            end
        end
    end
end
 
#= 
@views @kernel inbounds = true function nonlocal(W,w,mpD,meD,ls,type)
    p = @index(Global)

    if type == "p->q" && p ≤ mpD.nmp
        for q ∈ mpD.nmp
            ξ = (mpD.x[p,1]-mpD.x[q,1])
            η = (mpD.x[p,2]-mpD.x[q,2])
            d = sqrt(ξ^2+η^2)
        
            ω₀     = d/ls*exp(-(d/ls)^2)
            w[p,q] = ω₀
            w[q,p] = ω₀
            W[p]  += ω₀
            W[q]  += ω₀
      
        end
    elseif type == "p<-q" && p ≤ mpD.nmp
        mpD.ϵpII[p,2] = mpD.ϵpII[p,1]
        for q ∈ mpD.nmp
            if W[p] != 0.0
                mpD.ϵpII[p,2]+= (w[p,q]/W[p])*mpD.ϵpII[q,1]
            end
        end
    end
end
=#