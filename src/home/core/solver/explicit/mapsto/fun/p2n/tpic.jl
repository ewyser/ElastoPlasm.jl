# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    p2n!(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, basis::Basis{...,TR}, g::Vector{T2}) where {T1,T2,TR<:TpicTransfer}

Project material point data to mesh nodes (TPIC scheme). Dispatched via `Basis`'s `TR`
type parameter — see `Basis`'s docstring.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.

# Returns 
- Updates mesh fields in-place.
"""
@kernel inbounds = true function p2n!(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},basis::Basis{T1,T2,1,NN,K,TR},g::Vector{T2}) where {T1,T2,NN,K,TR<:TpicTransfer}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms, Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        xp     = mpts.x[p]        
        vp, ∇v = mpts.s.v[p]          , mpts.s.∇vᵢⱼ[p]
        σ      = get_voigt(mpts.s.σᵢⱼ[p]) 
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(1))
            δx    = mesh.x[no] - xp
            mv    = N * ms * (vp + ∇v * δx)
            # accumulation
            @atom mesh.s.m[no]   += N * ms
            @atom mesh.s.mv[no]  += mv[1]
            @atom mesh.s.oobf[no]-= Ω * (∂N[1] * σ[1]) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function p2n!(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2,NN,K,TR},g::Vector{T2}) where {T1,T2,NN,K,TR<:TpicTransfer}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms, Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        xp     = mpts.x[p]        
        vp, ∇v = mpts.s.v[p]          , mpts.s.∇vᵢⱼ[p]
        σ      = get_voigt(mpts.s.σᵢⱼ[p])      
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(2))
            δx    = mesh.x[no] - xp
            mv    = N * ms * (vp + ∇v * δx)
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += mv[1]
            @atom mesh.s.mv[2,no]  += mv[2]
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1] * σ[1] + ∂N[2] * σ[3])
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1] * σ[3] + ∂N[2] * σ[2]) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function p2n!(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},basis::Basis{T1,T2,3,NN,K,TR},g::Vector{T2}) where {T1,T2,NN,K,TR<:TpicTransfer}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms, Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        xp     = mpts.x[p]        
        vp, ∇v = mpts.s.v[p]          , mpts.s.∇vᵢⱼ[p]
        σ      = get_voigt(mpts.s.σᵢⱼ[p])
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(3))
            δx    = mesh.x[no] - xp
            mv    = N * ms * (vp + ∇v * δx)
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += mv[1]
            @atom mesh.s.mv[2,no]  += mv[2]
            @atom mesh.s.mv[3,no]  += mv[3]
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1]*σ[1]+∂N[2]*σ[6]+∂N[3]*σ[5])
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1]*σ[6]+∂N[2]*σ[2]+∂N[3]*σ[4])
            @atom mesh.s.oobf[3,no]-= Ω * (∂N[1]*σ[5]+∂N[2]*σ[4]+∂N[3]*σ[3]) - N*(ms*g[3])
        end
    end
end