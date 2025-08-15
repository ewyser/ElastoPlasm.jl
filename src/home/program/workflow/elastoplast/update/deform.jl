@views @kernel inbounds = true function finite_deform(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # compute velocity & displacement gradients
        mpts.s.∇vᵢⱼ[:,:,p].= T2(0.0)
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            for i ∈ 1:mesh.dim , j ∈ 1:mesh.dim
                mpts.s.∇vᵢⱼ[i,j,p]+= mpts.ϕ∂ϕ[nn,p,j+1]*mesh.v[i,no]
            end
        end
        # compute incremental deformation and update
        mpts.s.ΔFᵢⱼ[:,:,p].= mpts.δᵢⱼ+(dt.*mpts.s.∇vᵢⱼ[:,:,p])
        mpts.s.Fᵢⱼ[:,:,p] .= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpts.ΔJ[p]       = det(mpts.s.ΔFᵢⱼ[:,:,p])
        mpts.J[p]        = det(mpts.s.Fᵢⱼ[:,:,p])
        mpts.Ω[p]        = mpts.J[p]*mpts.Ω₀[p]
        # update material point's positivity-preserving porosity
        mpts.s.ρ[p]      = mpts.s.ρ₀[p]/mpts.J[p]     
        mpts.n[p]        = T2(1.0)-T2(1.0)/mpts.J[p]*(T2(1.0)-mpts.n₀[p])
    end
end
@views @kernel inbounds = true function infinitesimal_deform(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # compute velocity & displacement gradients
        mpts.s.∇vᵢⱼ[:,:,p].= T2(0.0)
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            for i ∈ 1:mesh.dim , j ∈ 1:mesh.dim
                mpts.s.∇vᵢⱼ[i,j,p]+= mpts.ϕ∂ϕ[nn,p,j+1]*mesh.v[i,no]
            end
        end
        # compute incremental deformation and update
        mpts.s.ΔFᵢⱼ[:,:,p].= mpts.δᵢⱼ+(dt.*mpts.s.∇vᵢⱼ[:,:,p])
        mpts.s.Fᵢⱼ[:,:,p] .= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpts.ΔJ[p]       = det(mpts.s.ΔFᵢⱼ[:,:,p])
        mpts.J[p]        = det(mpts.s.Fᵢⱼ[:,:,p])
        mpts.Ω[p]        = mpts.J[p]*mpts.Ω₀[p]
    end
end