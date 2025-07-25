@views @kernel inbounds = true function deform(mp,mesh,dt)
    p = @index(Global)
    if p≤mp.nmp 
        # compute velocity & displacement gradients
        mp.∇vᵢⱼ[:,:,p].= 0.0
        for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
            for i ∈ 1:mesh.dim , j ∈ 1:mesh.dim
                mp.∇vᵢⱼ[i,j,p]+= mp.ϕ∂ϕ[nn,p,j+1]*mesh.v[i,no]
            end
        end
        # compute incremental deformation and update
        mp.ΔFᵢⱼ[:,:,p].= mp.δᵢⱼ+(dt.*mp.∇vᵢⱼ[:,:,p])
        mp.Fᵢⱼ[:,:,p] .= mp.ΔFᵢⱼ[:,:,p]*mp.Fᵢⱼ[:,:,p]
        # update material point's volume
        mp.ΔJ[p]       = det(mp.ΔFᵢⱼ[:,:,p])
        mp.J[p]        = det(mp.Fᵢⱼ[:,:,p])
        mp.Ω[p]        = mp.J[p]*mp.Ω₀[p]
    end
end

































