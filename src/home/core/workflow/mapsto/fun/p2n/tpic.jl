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
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx    = mpts.s.v[1,p]        
        σxx   = mpts.s.σᵢ[1,p] 
        ∇vxx  = mpts.s.∇vᵢⱼ[1,1,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no    = mpts.p2n[nn,p]
            N,∂Nx = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2]
            δx    = mpts.Δnp[nn,1,p]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * ms
            @atom mesh.mv[1,no]  += N * ms * (vx + ∇vxx * δx)
            @atom mesh.oobf[1,no]-= Ω * (∂Nx * σxy) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω       = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx ,vy      = mpts.s.v[1,p]        ,mpts.s.v[2,p]
        σxx,σyy,σxy = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p]   ,mpts.s.σᵢ[3,p]
        ∇vxx,∇vxy   = mpts.s.∇vᵢⱼ[1,1,p]  ,mpts.s.∇vᵢⱼ[1,2,p]
        ∇vyx,∇vyy   = mpts.s.∇vᵢⱼ[2,1,p]  ,mpts.s.∇vᵢⱼ[2,2,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no        = mpts.p2n[nn,p]
            N,∂Nx,∂Ny = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3]
            δx,δy     = mpts.Δnp[nn,1,p],mpts.Δnp[nn,2,p]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * ms
            @atom mesh.mv[1,no]  += N * ms * (vx + ∇vxx * δx + ∇vxy * δy)
            @atom mesh.mv[2,no]  += N * ms * (vy + ∇vyx * δx + ∇vyy * δy)
            @atom mesh.oobf[1,no]-= Ω * (∂Nx * σxx + ∂Ny * σxy)
            @atom mesh.oobf[2,no]-= Ω * (∂Nx * σxy + ∂Ny * σyy) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function tpic_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms  ,Ω          = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        vx  ,vy  ,vz    = mpts.s.v[1,p]        ,mpts.s.v[2,p]  ,mpts.s.v[3,p]
        σxx ,σyy ,σzz   = mpts.s.σᵢ[1,p]       ,mpts.s.σᵢ[2,p] ,mpts.s.σᵢ[3,p]
        σyx ,σzy ,σzx   = mpts.s.σᵢ[6,p]       ,mpts.s.σᵢ[4,p] ,mpts.s.σᵢ[5,p]
        ∇vxx,∇vxy,∇vxz = mpts.s.∇vᵢⱼ[1,1,p]   ,mpts.s.∇vᵢⱼ[1,2,p],mpts.s.∇vᵢⱼ[1,3,p]
        ∇vyx,∇vyy,∇vyz = mpts.s.∇vᵢⱼ[2,1,p]   ,mpts.s.∇vᵢⱼ[2,2,p],mpts.s.∇vᵢⱼ[2,3,p]
        ∇vzx,∇vzy,∇vzz = mpts.s.∇vᵢⱼ[3,1,p]   ,mpts.s.∇vᵢⱼ[3,2,p],mpts.s.∇vᵢⱼ[3,3,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no            = mpts.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3],mpts.ϕ∂ϕ[nn,p,4]
            δx,δy,δz      = mpts.Δnp[nn,1,p],mpts.Δnp[nn,2,p],mpts.Δnp[nn,3,p]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * ms
            @atom mesh.mv[1,no]  += N * ms * (vx + ∇vxx * δx + ∇vxy * δy + ∇vxz * δz)
            @atom mesh.mv[2,no]  += N * ms * (vy + ∇vyx * δx + ∇vyy * δy + ∇vyz * δz)
            @atom mesh.mv[3,no]  += N * ms * (vz + ∇vzx * δx + ∇vzy * δy + ∇vzz * δz)
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (ms * g[3])
        end
    end
end