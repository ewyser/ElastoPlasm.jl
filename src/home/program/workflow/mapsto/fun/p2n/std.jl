# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FLIP transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    std_1d_p2n(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, g::Vector{T2}) where {T1,T2}

Project 1D material point data to mesh nodes (FLIP scheme).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.

# Returns
- Updates mesh fields in-place.
"""
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # buffering 
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        px    = ms*mpts.s.v[p]
        σxx   = mpts.s.σᵢ[1,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no        = mpts.p2n[nn,p]
            N,∂Nx     = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]  += N * ms
            @atom mesh.mv[no]  += N * px
            @atom mesh.oobf[no]-= Ω * (∂Nx * σxx)
            @atom mesh.oobf[no]+= N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshThermalPhase{T1,T2,D}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms ,Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  ,T  = mpts.t.c[p]          , mpts.t.T[p]
        qx     = mpts.t.q[1,p]        
        γ      = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no     = mpts.p2n[nn,p]
            N, ∂Nx = mpts.ϕ∂ϕ[nn,p,1], mpts.ϕ∂ϕ[nn,p,2]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        px , py       = ms*mpts.s.v[1,p]     , ms*mpts.s.v[2,p]
        σxx, σyy, σxy = mpts.s.σᵢ[1,p]       , mpts.s.σᵢ[2,p] , mpts.s.σᵢ[3,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no        = mpts.p2n[nn,p]
            N,∂Nx,∂Ny = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * ms
            @atom mesh.mv[1,no]  += N * px
            @atom mesh.mv[2,no]  += N * py
            @atom mesh.oobf[1,no]-= Ω * (∂Nx * σxx + ∂Ny * σxy)
            @atom mesh.oobf[2,no]-= Ω * (∂Nx * σxy + ∂Ny * σyy) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshThermalPhase{T1,T2,D}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms , Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  , T  = mpts.t.c[p]          , mpts.t.T[p]
        qx , qy = mpts.t.q[1,p]        , mpts.t.q[2,p]
        γ       = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no        = mpts.p2n[nn,p]
            N,∂Nx,∂Ny = mpts.ϕ∂ϕ[nn,p,1],mpts.ϕ∂ϕ[nn,p,2],mpts.ϕ∂ϕ[nn,p,3]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx + ∂Ny * qy)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        px , py , pz  = ms*mpts.s.v[1,p]     , ms*mpts.s.v[2,p], ms*mpts.s.v[3,p]
        σxx, σyy, σzz = mpts.s.σᵢ[1,p]       , mpts.s.σᵢ[2,p]  , mpts.s.σᵢ[3,p]
        σyx, σzy, σzx = mpts.s.σᵢ[6,p]       , mpts.s.σᵢ[4,p]  , mpts.s.σᵢ[5,p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no            = mpts.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mpts.ϕ∂ϕ[nn,p,1], mpts.ϕ∂ϕ[nn,p,2], mpts.ϕ∂ϕ[nn,p,3], mpts.ϕ∂ϕ[nn,p,4]
            # accumulation
            if iszero(no) continue end
            @atom mesh.mᵢ[no]    += N * ms
            @atom mesh.mv[1,no]  += N * px
            @atom mesh.mv[2,no]  += N * py
            @atom mesh.mv[3,no]  += N * pz
            @atom mesh.oobf[1,no]-= Ω * ( ∂Nx * σxx + ∂Ny * σyx + ∂Nz * σzx)
            @atom mesh.oobf[2,no]-= Ω * ( ∂Nx * σyx + ∂Ny * σyy + ∂Nz * σzy)
            @atom mesh.oobf[3,no]-= Ω * ( ∂Nx * σzx + ∂Ny * σzy + ∂Nz * σzz) - N * (ms * g[3])
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,D},mesh::MeshThermalPhase{T1,T2,D}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms , Ω      = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  , T      = mpts.t.c[p]          , mpts.t.T[p]
        qx , qy, qz = mpts.t.q[1,p]        , mpts.t.q[2,p], mpts.t.q[3,p]
        γ           = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering 
            no            = mpts.p2n[nn,p]
            N,∂Nx,∂Ny,∂Nz = mpts.ϕ∂ϕ[nn,p,1], mpts.ϕ∂ϕ[nn,p,2], mpts.ϕ∂ϕ[nn,p,3], mpts.ϕ∂ϕ[nn,p,4]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx + ∂Ny * qy + ∂Nz * qz)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end
