@kernel inbounds = true function deform(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},dt::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        # accumulate ∇vᵢⱼ into a local mutable buffer then convert
        TM = eltype(mpts.s.∇vᵢⱼ)
        ∇v = zeros(MMatrix{size(TM,1),size(TM,2),T2})
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            ∂N = ∂Nrow(basis.∂N, nn, p, Val(D))
            for i ∈ 1:(length(mesh.prprt.nel)-1)
                for j ∈ 1:(length(mesh.prprt.nel)-1)
                    ∇v[i,j] += ∂N[j]*mesh.s.v[no][i]
                end
            end
        end
        ∇vᵢⱼ   = TM(∇v)
        ΔFᵢⱼ   = TM(I) + dt * ∇vᵢⱼ
        Fᵢⱼ    = ΔFᵢⱼ * mpts.s.Fᵢⱼ[p]
        mpts.s.∇vᵢⱼ[p] = ∇vᵢⱼ
        mpts.s.ΔFᵢⱼ[p] = ΔFᵢⱼ
        mpts.s.Fᵢⱼ[p]  = Fᵢⱼ
        ΔJ = det(ΔFᵢⱼ)
        J  = det(Fᵢⱼ)
        mpts.ΔJ[p]  = ΔJ
        mpts.J[p]   = J
        mpts.Ω[p]   = J * mpts.Ω₀[p]
        n           = T2(1.0) - T2(1.0)/J*(T2(1.0)-mpts.n₀[p])
        
        mpts.s.ρ[p] = (T2(1.0)-n)*mpts.s.ρ₀[p]
        mpts.n[p]   = min(max(n, T2(0.0)), T2(1.0))
    end
end
