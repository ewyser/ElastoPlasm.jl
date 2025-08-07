@views @kernel inbounds = true function transform(mp::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    # deformation framework dispatcher
    if p ≤ mp.nmp 
        J = mp.s.J[p]
        for c ∈ 1:size(mp.s.σᵢ,1)
            mp.s.σᵢ[c,p] = mp.s.τᵢ[c,p]/J
        end
    end   
end
@kernel inbounds = true function flip_1d_p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # buffering 
        m,Ω = mp.s.m[p]    ,mp.Ω[p]
        px  = m*mp.s.v[p]
        σxx = mp.s.σᵢ[1,p]
        for nn ∈ 1:mesh.nn
            # buffering 
            no        = mp.p2n[nn,p]
            N,∂Nx     = mp.ϕ∂ϕ[nn,p,1],mp.ϕ∂ϕ[nn,p,2]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]  += N * m
            @atom mesh.p[no]   += N * px
            @atom mesh.oobf[no]-= Ω * (∂Nx * σxx)
            @atom mesh.oobf[no]+= N * (m * g[1])
        end
    end
end
@kernel inbounds = true function flip_2d_p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp
        # buffering 
        m  ,Ω       = mp.s.m[p]    ,mp.Ω[p]
        px ,py      = m*mp.s.v[1,p],m*mp.s.v[2,p]
        σxx,σyy,σxy = mp.s.σᵢ[1,p] ,mp.s.σᵢ[2,p] ,mp.s.σᵢ[3,p]
        for nn ∈ 1:mesh.nn
            # buffering 
            no        = mp.p2n[nn,p]
            N,∂Nx,∂Ny = mp.ϕ∂ϕ[nn,p,1],mp.ϕ∂ϕ[nn,p,2],mp.ϕ∂ϕ[nn,p,3]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * m
            @atom mesh.p[1,no]   += N * px
            @atom mesh.p[2,no]   += N * py
            @atom mesh.oobf[1,no]-= Ω * (∂Nx * σxx + ∂Ny * σxy)
            @atom mesh.oobf[2,no]-= Ω * (∂Nx * σxy + ∂Ny * σyy) - N * (m * g[2])
        end
    end
end
@kernel inbounds = true function flip_3d_p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp
        # buffering 
        m  ,Ω       = mp.s.m[p]    ,mp.Ω[p]
        px ,py ,pz  = m*mp.s.v[1,p],m*mp.s.v[2,p],m*mp.s.v[3,p]
        σxx,σyy,σzz = mp.s.σᵢ[1,p] ,mp.s.σᵢ[2,p] ,mp.s.σᵢ[3,p]
        σyx,σzy,σzx = mp.s.σᵢ[6,p] ,mp.s.σᵢ[4,p] ,mp.s.σᵢ[5,p]
        for nn ∈ 1:mesh.nn
            # buffering
            no            = mp.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mp.ϕ∂ϕ[nn,p,1],mp.ϕ∂ϕ[nn,p,2],mp.ϕ∂ϕ[nn,p,3],mp.ϕ∂ϕ[nn,p,4]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * m
            @atom mesh.p[1,no]   += N * px
            @atom mesh.p[2,no]   += N * py
            @atom mesh.p[3,no]   += N * pz
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (m * g[3])
        end
    end
end
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
@kernel inbounds = true function tpic_2d_p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp
        # buffering 
        m  ,Ω       = mp.s.m[p]       ,mp.Ω[p]
        vx ,vy      = mp.s.v[1,p]     ,mp.s.v[2,p]
        σxx,σyy,σxy = mp.s.σᵢ[1,p]    ,mp.s.σᵢ[2,p]     ,mp.s.σᵢ[3,p]
        ∇vxx,∇vxy   = mp.s.∇vᵢⱼ[1,1,p],mp.s.∇vᵢⱼ[1,2,p]
        ∇vyx,∇vyy   = mp.s.∇vᵢⱼ[2,1,p],mp.s.∇vᵢⱼ[2,2,p]
        for nn ∈ 1:mesh.nn
            # buffering 
            no        = mp.p2n[nn,p]
            N,∂Nx,∂Ny = mp.ϕ∂ϕ[nn,p,1],mp.ϕ∂ϕ[nn,p,2],mp.ϕ∂ϕ[nn,p,3]
            δx,δy     = mp.δnp[nn,1,p],mp.δnp[nn,2,p]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * m
            @atom mesh.p[1,no]   += N * m * (vx + ∇vxx * δx + ∇vxy * δy)
            @atom mesh.p[2,no]   += N * m * (vy + ∇vyx * δx + ∇vyy * δy)
            @atom mesh.oobf[1,no]-= Ω * (∂Nx * σxx + ∂Ny * σxy)
            @atom mesh.oobf[2,no]-= Ω * (∂Nx * σxy + ∂Ny * σyy) - N * (m * g[2])
        end
    end
end
@kernel inbounds = true function tpic_3d_p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp
        # buffering 
        m   ,Ω          = mp.s.m[p]       ,mp.Ω[p]
        vx  ,vy  ,vz    = mp.s.v[1,p]     ,mp.s.v[2,p]  ,mp.s.v[3,p]
        σxx ,σyy ,σzz   = mp.s.σᵢ[1,p]    ,mp.s.σᵢ[2,p] ,mp.s.σᵢ[3,p]
        σyx ,σzy ,σzx   = mp.s.σᵢ[6,p]    ,mp.s.σᵢ[4,p] ,mp.s.σᵢ[5,p]
        σxx ,σyy ,σxy   = mp.s.σᵢ[1,p]    ,mp.s.σᵢ[2,p] ,mp.s.σᵢ[3,p]
        ∇vxx,∇vxy,∇vxz = mp.s.∇vᵢⱼ[1,1,p],mp.s.∇vᵢⱼ[1,2,p],mp.s.∇vᵢⱼ[1,3,p]
        ∇vyx,∇vyy,∇vyz = mp.s.∇vᵢⱼ[2,1,p],mp.s.∇vᵢⱼ[2,2,p],mp.s.∇vᵢⱼ[2,3,p]
        ∇vzx,∇vzy,∇vzz = mp.s.∇vᵢⱼ[3,1,p],mp.s.∇vᵢⱼ[3,2,p],mp.s.∇vᵢⱼ[3,3,p]
        for nn ∈ 1:mesh.nn
            # buffering
            no            = mp.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mp.ϕ∂ϕ[nn,p,1],mp.ϕ∂ϕ[nn,p,2],mp.ϕ∂ϕ[nn,p,3],mp.ϕ∂ϕ[nn,p,4]
            δx,δy,δz      = mp.δnp[nn,1,p],mp.δnp[nn,2,p],mp.δnp[nn,3,p]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * m
            @atom mesh.p[1,no]   += N * m * (vx + ∇vxx * δx + ∇vxy * δy + ∇vxz * δz)
            @atom mesh.p[2,no]   += N * m * (vy + ∇vyx * δx + ∇vyy * δy + ∇vyz * δz)
            @atom mesh.p[3,no]   += N * m * (vz + ∇vzx * δx + ∇vzy * δy + ∇vzz * δz)
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (m * g[3])
        end
    end
end
function p2n(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2},instr::Dict) where {T1,T2}
    # get cauchy stress 
    if instr[:fwrk][:deform] == "finite"
        instr[:cairn][:mapsto][:map].σᵢ!(ndrange=mp.nmp,mp);sync(CPU())
    end
    # initialize nodal quantities
    mesh.mᵢ  .= T2(0.0)
    mesh.p   .= T2(0.0)
    mesh.oobf.= T2(0.0)
    # mapping to mesh
    instr[:cairn][:mapsto][:map].p2n!(ndrange=mp.nmp,mp,mesh,g);sync(CPU())
    return nothing
end