# compute ∇vᵢⱼ and ΔFᵢⱼ from current nodal velocities WITHOUT side effects
# F, J, ΔJ, ρ, Ω, n are intentionally not updated — this is a trial/linearized computation
# used inside the PT implicit solve loop

@views @kernel inbounds = true function trial_deform(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        ∇v = zeros(MMatrix{2,2,T2})
        for nn ∈ 1:mesh.prprt.nn
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            for i ∈ 1:mesh.prprt.dim
                for j ∈ 1:mesh.prprt.dim
                    ∇v[i,j] += ∂N[j]*mesh.s.v[i,no]
                end
            end
        end
        ∇vᵢⱼ = SMatrix(∇v)
        mpts.s.∇vᵢⱼ[p] = ∇vᵢⱼ
        mpts.s.ΔFᵢⱼ[p] = SMatrix{2,2,T2}(I) + dt .* ∇vᵢⱼ
    end
end
