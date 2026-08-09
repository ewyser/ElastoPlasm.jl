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
@kernel inbounds = true function assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # buffering 
        ms ,Ω = mpts.s.ρ[p]*mpts.Ω[p],mpts.Ω[p]
        px    = ms*mpts.s.v[p][1]
        σxx   = mpts.s.σᵢ[p][1]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no]  += N * ms
            @atom mesh.s.mv[no]  += N * px
            #@atom mesh.s.oobf[no]-= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end
@kernel inbounds = true function mass_assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms = mpts.s.ρ[p]*mpts.Ω[p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no] += N * ms
        end
    end
end
@kernel inbounds = true function assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        px , py , pz  = ms*mpts.s.v[p][1]     , ms*mpts.s.v[p][2], ms*mpts.s.v[p][3]
        σxx, σyy, σzz = mpts.s.σᵢ[p][1]       , mpts.s.σᵢ[p][2]  , mpts.s.σᵢ[p][3]
        σyx, σzy, σzx = mpts.s.σᵢ[p][6]       , mpts.s.σᵢ[p][4]  , mpts.s.σᵢ[p][5]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no]     += N * ms
            @atom mesh.s.mv[1,no]  += N * px
            @atom mesh.s.mv[2,no]  += N * py
            @atom mesh.s.mv[3,no]  += N * pz
            @atom mesh.s.oobf[1,no]-= Ω * ( ∂N[1] * σxx + ∂N[2] * σyx + ∂N[3] * σzx)
            @atom mesh.s.oobf[2,no]-= Ω * ( ∂N[1] * σyx + ∂N[2] * σyy + ∂N[3] * σzy)
            @atom mesh.s.oobf[3,no]-= Ω * ( ∂N[1] * σzx + ∂N[2] * σzy + ∂N[3] * σzz) - N * (ms * g[3])
        end
    end
end