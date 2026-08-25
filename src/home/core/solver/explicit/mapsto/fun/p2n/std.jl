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
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},basis::Basis{T1,T2,1},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        px    = ms*mpts.s.v[p][1]
        σxx   = get_vector(mpts.s.σᵢ[p])[1]
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(1))
            # accumulation
            @atom mesh.s.m[no]  += N * ms
            @atom mesh.s.mv[no]  += N * px
            @atom mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        mv            = ms*mpts.s.v[p]
        σ             = get_vector(mpts.s.σᵢ[p])
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(2))
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += N * mv[1]
            @atom mesh.s.mv[2,no]  += N * mv[2]
            @atom mesh.s.oobf[1,no]-= Ω * (∂N[1] * σ[1] + ∂N[2] * σ[3])
            @atom mesh.s.oobf[2,no]-= Ω * (∂N[1] * σ[3] + ∂N[2] * σ[2]) - N * (ms * g[2])
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},basis::Basis{T1,T2,3},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        mv            = ms*mpts.s.v[p]
        σ             = get_vector(mpts.s.σᵢ[p])
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N, ∂N = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(3))
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += N * mv[1]
            @atom mesh.s.mv[2,no]  += N * mv[2]
            @atom mesh.s.mv[3,no]  += N * mv[3]
            @atom mesh.s.oobf[1,no]-= Ω * ( ∂N[1] * σ[1] + ∂N[2] * σ[6] + ∂N[3] * σ[5])
            @atom mesh.s.oobf[2,no]-= Ω * ( ∂N[1] * σ[6] + ∂N[2] * σ[2] + ∂N[3] * σ[4])
            @atom mesh.s.oobf[3,no]-= Ω * ( ∂N[1] * σ[5] + ∂N[2] * σ[4] + ∂N[3] * σ[3]) - N * (ms * g[3])
        end
    end
end














@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,1},mesh::MeshThermalPhase{T1,T2,1},basis::Basis{T1,T2,1}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms ,Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  ,T  = mpts.t.c[p]          , mpts.t.T[p]
        qx     = mpts.t.q[1,p]
        γ      = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no     = basis.p2n[p][nn]
            N, ∂N  = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(1))
            ∂Nx    = ∂N[1]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end

@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,2},mesh::MeshThermalPhase{T1,T2,2},basis::Basis{T1,T2,2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms , Ω  = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  , T  = mpts.t.c[p]          , mpts.t.T[p]
        qx , qy = mpts.t.q[1,p]        , mpts.t.q[2,p]
        γ       = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no          = basis.p2n[p][nn]
            N, ∂N       = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(2))
            ∂Nx,∂Ny     = ∂N[1], ∂N[2]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx + ∂Ny * qy)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end
@kernel inbounds = true function std_p2n(mpts::Point{T1,T2,3},mesh::MeshThermalPhase{T1,T2,3},basis::Basis{T1,T2,3}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering
        ms , Ω      = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        c  , T      = mpts.t.c[p]          , mpts.t.T[p]
        qx , qy, qz = mpts.t.q[1,p]        , mpts.t.q[2,p], mpts.t.q[3,p]
        γ           = T2(0.0) # heat source
        for nn ∈ 1:mesh.prprt.nn
            # buffering
            no              = basis.p2n[p][nn]
            N, ∂N           = basis.N[nn,p], ∂Nrow(basis.∂N, nn, p, Val(3))
            ∂Nx,∂Ny,∂Nz     = ∂N[1], ∂N[2], ∂N[3]
            # accumulation
            if iszero(no) continue end
            @atom mesh.cᵢ[no]  += N * ms * c
            @atom mesh.mcT[no] += N * ms * c * T
            @atom mesh.oobq[no]+= Ω * (∂Nx * qx + ∂Ny * qy + ∂Nz * qz)
            #@atom mesh.oobq[no]+= Ω * γ * N
        end
    end
end

