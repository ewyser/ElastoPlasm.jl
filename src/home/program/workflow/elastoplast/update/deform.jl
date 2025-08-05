@views @kernel inbounds = true function finite_deform(mp::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # compute velocity & displacement gradients
        mp.s.∇vᵢⱼ[:,:,p].= T2(0.0)
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if iszero(no) continue end
            for i ∈ 1:mesh.dim , j ∈ 1:mesh.dim
                mp.s.∇vᵢⱼ[i,j,p]+= mp.ϕ∂ϕ[nn,p,j+1]*mesh.v[i,no]
            end
        end
        # compute incremental deformation and update
        mp.s.ΔFᵢⱼ[:,:,p].= mp.s.δᵢⱼ+(dt.*mp.s.∇vᵢⱼ[:,:,p])
        mp.s.Fᵢⱼ[:,:,p] .= mp.s.ΔFᵢⱼ[:,:,p]*mp.s.Fᵢⱼ[:,:,p]
        # update material point's volume
        mp.s.ΔJ[p]       = det(mp.s.ΔFᵢⱼ[:,:,p])
        mp.s.J[p]        = det(mp.s.Fᵢⱼ[:,:,p])
        mp.Ω[p]          = mp.s.J[p]*mp.Ω₀[p]
    end
end
@views @kernel inbounds = true function infinitesimal_deform(mp::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # compute velocity & displacement gradients
        mp.s.∇vᵢⱼ[:,:,p].= T2(0.0)
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if iszero(no) continue end
            for i ∈ 1:mesh.dim , j ∈ 1:mesh.dim
                mp.s.∇vᵢⱼ[i,j,p]+= mp.ϕ∂ϕ[nn,p,j+1]*mesh.v[i,no]
            end
        end
        # compute incremental deformation and update
        mp.s.ΔFᵢⱼ[:,:,p].= mp.s.δᵢⱼ+(dt.*mp.s.∇vᵢⱼ[:,:,p])
        mp.s.Fᵢⱼ[:,:,p] .= mp.s.ΔFᵢⱼ[:,:,p]*mp.s.Fᵢⱼ[:,:,p]
        # update material point's volume
        mp.s.ΔJ[p]       = det(mp.s.ΔFᵢⱼ[:,:,p])
        mp.s.J[p]        = det(mp.s.Fᵢⱼ[:,:,p])
        mp.Ω[p]          = mp.s.J[p]*mp.Ω₀[p]
    end
end
































