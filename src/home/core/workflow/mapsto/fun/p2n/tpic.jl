# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    tpic_1d_p2n(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, g::Vector{T2}) where {T1,T2}

Project 1D material point data to mesh nodes (TPIC scheme).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.

# Returns 
- Updates mesh fields in-place.
"""
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx    = mpts.s.v[p][1]
        σxx   = mpts.s.σᵢ[p][1]
        ∇vxx  = mpts.s.∇vᵢⱼ[p][1,1]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            δx        = mesh.x[no][1]-mpts.x[1,p]
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no]   += N * ms
            @atom mesh.s.mv[no]  += N * ms * (vx + ∇vxx * δx)
            @atom mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        v_p   = mpts.s.v[p]
        σ     = mpts.s.σᵢ[p]
        ∇v    = mpts.s.∇vᵢⱼ[p]          # SMatrix: (∇v * δx) gives TPIC correction
        x_p   = SVector{2,T2}(mpts.x[1,p], mpts.x[2,p])  # zero-alloc, hoisted
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δx = mesh.x[no] - x_p
            mv  = N * ms * (v_p + ∇v * δx)
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += mv[1]
            @atom mesh.s.mv[2,no]  += mv[2]
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1] * σ[1] + ∂N[2] * σ[3])
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1] * σ[3] + ∂N[2] * σ[2]) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        v_p   = mpts.s.v[p]
        σ     = mpts.s.σᵢ[p]
        ∇v    = mpts.s.∇vᵢⱼ[p]
        x_p   = SVector{3,T2}(mpts.x[1,p], mpts.x[2,p], mpts.x[3,p])
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δx = mesh.x[no] - x_p
            mv  = N * ms * (v_p + ∇v * δx)
            @atom mesh.s.m[no]   += N * ms
            @atom mesh.s.mv[1,no]+= mv[1]
            @atom mesh.s.mv[2,no]+= mv[2]
            @atom mesh.s.mv[3,no]+= mv[3]
            @atom mesh.oobf[1,no]-= Ω * (∂N[1]*σ[1]+∂N[2]*σ[6]+∂N[3]*σ[5])
            @atom mesh.oobf[2,no]-= Ω * (∂N[1]*σ[6]+∂N[2]*σ[2]+∂N[3]*σ[4])
            @atom mesh.oobf[3,no]-= Ω * (∂N[1]*σ[5]+∂N[2]*σ[4]+∂N[3]*σ[3]) - N*(ms*g[3])
        end
    end
end