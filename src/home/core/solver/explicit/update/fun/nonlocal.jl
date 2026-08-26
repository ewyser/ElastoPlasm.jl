#= =#
"""
    nonlocal(W, w, mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, ls::T2, type::String) where {T1,T2}

Apply nonlocal averaging for regularization of plastic strain at material points.

# Arguments
- `W`: Vector of normalization weights.
- `w`: Matrix of pairwise weights.
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `ls::T2`: Nonlocal length scale.
- `type::String`: Operation type ("tplgy", "p->q", or "p<-q").

# Returns
- Updates weights and plastic strain fields in-place.
"""
@views @kernel inbounds = true function nonlocal(W,w,mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},ls::T2,type::String) where {T1,T2}
    p = @index(Global)

    if type == "tplgy" && p ≤ mpts.nmp
        for k ∈ nzrange(basis.e2e, basis.p2e[p])
            el = rowvals(basis.e2e)[k]
            basis.e2p[p,el] = p
        end
    elseif type == "p->q" && p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        el = basis.p2e[p]
        for q ∈ axes(basis.e2p,1)
            if !iszero(basis.e2p[q,el])
                d   = norm(mpts.x[p] - mpts.x[q])
                if w[p,q] == T2(0.0)
                    ω₀     = d/ls*exp(-(d/ls)^2)
                    w[p,q] = ω₀
                    w[q,p] = ω₀
                    W[p]  += ω₀
                    W[q]  += ω₀
                    basis.p2p[q,p] = q
                end
            end
        end
    elseif type == "p<-q" && p ≤ mpts.nmp && mpts.s.Δλ[p] != T2(0.0)
        if isapprox(W[p]>T2(1e-16),T2(0.0),atol=T2(1e-16))
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], mpts.s.ϵpII[p][1])
        else
            acc = mpts.s.ϵpII[p][2]
            for q ∈ axes(basis.p2p,1)
                if !iszero(basis.p2p[q,p])
                    acc += (w[p,q]/W[p])*mpts.s.ϵpII[q][1]
                end
            end
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], acc)
        end
    end
end
 
#= 
@views @kernel inbounds = true function nonlocal(W,w,mpts,mesh,ls,type)
    p = @index(Global)

    if type == "p->q" && p ≤ mpts.nmp
        for q ∈ mpts.nmp
            ξ = (mpts.x[p,1]-mpts.x[q,1])
            η = (mpts.x[p,2]-mpts.x[q,2])
            d = sqrt(ξ^2+η^2)
        
            ω₀     = d/ls*exp(-(d/ls)^2)
            w[p,q] = ω₀
            w[q,p] = ω₀
            W[p]  += ω₀
            W[q]  += ω₀
      
        end
    elseif type == "p<-q" && p ≤ mpts.nmp
        mpts.ϵpII[p,2] = mpts.ϵpII[p,1]
        for q ∈ mpts.nmp
            if W[p] != 0.0
                mpts.ϵpII[p,2]+= (w[p,q]/W[p])*mpts.ϵpII[q,1]
            end
        end
    end
end
=#