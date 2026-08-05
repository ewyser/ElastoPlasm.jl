# Incremental continuity residual: ∫N(div v + Δp/K)dΩ = 0  where Δp = p - p_n (pressure increment this step)

@kernel inbounds = true function div_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},Kc::T2,p_n::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        Ω    = mpts.Ω[p]
        divv = mpts.s.∇vᵢⱼ[1,1,p] + mpts.s.∇vᵢⱼ[2,2,p]
        Δpp  = mpts.s.p[p] - p_n[p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, _ = basis(mpts, mesh, T1(p), T1(nn))
            if iszero(no) continue end
            @atom mesh.s.oobp[no] += N * Ω * (Δpp / Kc + divv)
        end
    end
end

@kernel inbounds = true function div_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},Kc::T2,p_n::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        Ω    = mpts.Ω[p]
        divv = mpts.s.∇vᵢⱼ[1,1,p] + mpts.s.∇vᵢⱼ[2,2,p] + mpts.s.∇vᵢⱼ[3,3,p]
        Δpp  = mpts.s.p[p] - p_n[p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, _ = basis(mpts, mesh, T1(p), T1(nn))
            if iszero(no) continue end
            @atom mesh.s.oobp[no] += N * Ω * (Δpp / Kc + divv)
        end
    end
end
