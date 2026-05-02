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
        vx    = mpts.s.v[1,p]        
        σxx   = mpts.s.σᵢ[1,p] 
        ∇vxx  = mpts.s.∇vᵢⱼ[1,1,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            δx        = mesh.x[1,no]-mpts.x[1,p]
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
        # buffering 
        ms ,Ω       = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx ,vy      = mpts.s.v[1,p]        ,mpts.s.v[2,p]
        σxx,σyy,σxy = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p]   ,mpts.s.σᵢ[3,p]
        ∇vxx,∇vxy   = mpts.s.∇vᵢⱼ[1,1,p]   ,mpts.s.∇vᵢⱼ[1,2,p]
        ∇vyx,∇vyy   = mpts.s.∇vᵢⱼ[2,1,p]   ,mpts.s.∇vᵢⱼ[2,2,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            δx, δy    = mesh.x[1,no]-mpts.x[1,p], mesh.x[2,no]-mpts.x[2,p]
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += N * ms * (vx + ∇vxx * δx + ∇vxy * δy)
            @atom mesh.s.mv[2,no]  += N * ms * (vy + ∇vyx * δx + ∇vyy * δy)
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])

        end
    end
end
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms  ,Ω         = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx  ,vy  ,vz   = mpts.s.v[1,p]        ,mpts.s.v[2,p]     ,mpts.s.v[3,p]
        σxx ,σyy ,σzz  = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p]    ,mpts.s.σᵢ[3,p]
        σyx ,σzy ,σzx  = mpts.s.σᵢ[6,p]       ,mpts.s.σᵢ[4,p]    ,mpts.s.σᵢ[5,p]
        ∇vxx,∇vxy,∇vxz = mpts.s.∇vᵢⱼ[1,1,p]  ,mpts.s.∇vᵢⱼ[1,2,p],mpts.s.∇vᵢⱼ[1,3,p]
        ∇vyx,∇vyy,∇vyz = mpts.s.∇vᵢⱼ[2,1,p]  ,mpts.s.∇vᵢⱼ[2,2,p],mpts.s.∇vᵢⱼ[2,3,p]
        ∇vzx,∇vzy,∇vzz = mpts.s.∇vᵢⱼ[3,1,p]  ,mpts.s.∇vᵢⱼ[3,2,p],mpts.s.∇vᵢⱼ[3,3,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N  = basis(mpts, mesh, p, nn)
            δx, δy, δz = mesh.x[1,no]-mpts.x[1,p], mesh.x[2,no]-mpts.x[2,p], mesh.x[3,no]-mpts.x[3,p]
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no]   += N * ms
            @atom mesh.s.mv[1,no]+= N * ms * (vx + ∇vxx * δx + ∇vxy * δy + ∇vxz * δz)
            @atom mesh.s.mv[2,no]+= N * ms * (vy + ∇vyx * δx + ∇vyy * δy + ∇vyz * δz)
            @atom mesh.s.mv[3,no]+= N * ms * (vz + ∇vzx * δx + ∇vzy * δy + ∇vzz * δz)
            @atom mesh.oobf[1,no]-= Ω * ( ∂N[1] * σxx + ∂N[2] * σyx + ∂N[3] * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂N[1] * σyx + ∂N[2] * σyy + ∂N[3] * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂N[1] * σzx + ∂N[2] * σzy + ∂N[3] * σzz) - N * (ms * g[3])
        end
    end
end