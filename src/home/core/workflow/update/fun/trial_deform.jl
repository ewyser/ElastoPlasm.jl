# compute ∇vᵢⱼ and ΔFᵢⱼ from current nodal velocities WITHOUT side effects
# F, J, ΔJ, ρ, Ω, n are intentionally not updated — this is a trial/linearized computation
# used inside the PT implicit solve loop

@views @kernel inbounds = true function trial_deform(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        mpts.s.∇vᵢⱼ[:,:,p] .= T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            for i ∈ 1:mesh.prprt.dim
                for j ∈ 1:mesh.prprt.dim
                    mpts.s.∇vᵢⱼ[i,j,p] += ∂N[j]*mesh.s.v[i,no]
                end
            end
        end
        mpts.s.ΔFᵢⱼ[:,:,p] .= mpts.δᵢⱼ .+ (dt .* mpts.s.∇vᵢⱼ[:,:,p])
    end
end
