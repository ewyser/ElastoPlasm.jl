# accumulate only f_int and body force to nodes — mass and momentum are NOT touched
# used in PT sub-iterations where mass is frozen from the start of the time step

@kernel inbounds = true function fint_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms, Ω = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx   = mpts.s.σᵢ[1,p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end

@kernel inbounds = true function fint_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx, σyy, σxy = mpts.s.σᵢ[1,p], mpts.s.σᵢ[2,p], mpts.s.σᵢ[3,p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no] -= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end
    end
end

@kernel inbounds = true function fint_p2n(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx, σyy, σzz = mpts.s.σᵢ[1,p], mpts.s.σᵢ[2,p], mpts.s.σᵢ[3,p]
        σyx, σzy, σzx = mpts.s.σᵢ[6,p], mpts.s.σᵢ[4,p], mpts.s.σᵢ[5,p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * ( ∂N[1] * σxx + ∂N[2] * σyx + ∂N[3] * σzx)
            @atom mesh.s.oobf[2,no] -= Ω * ( ∂N[1] * σyx + ∂N[2] * σyy + ∂N[3] * σzy)
            @atom mesh.s.oobf[3,no] -= Ω * ( ∂N[1] * σzx + ∂N[2] * σzy + ∂N[3] * σzz) - N * (ms * g[3])
        end
    end
end
