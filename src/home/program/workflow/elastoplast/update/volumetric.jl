@kernel inbounds = true function ΔJn(mp,mesh)
    p = @index(Global)
    if p≤mp.nmp 
        # accumulation
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            @atom mesh.ΔJ[no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*mp.ΔJ[p])  
        end
    end
end
@kernel inbounds = true function ΔJs(mesh)
    no = @index(Global)
    if no≤mesh.nno[end] 
        # solve
        if iszero(mesh.mᵢ[no])
            mesh.ΔJ[no] = 0.0
        else
            mesh.ΔJ[no]/= mesh.mᵢ[no]
        end
    end
end
@views @kernel inbounds = true function ΔJp(mp,mesh,dim)
    p = @index(Global)
    if p≤mp.nmp 
        # mapping back to mp's
        ΔJ = 0.0
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            ΔJ += mp.ϕ∂ϕ[nn,p,1]*mesh.ΔJ[no]/mp.ΔJ[p]
        end
        # update
        mp.ΔFᵢⱼ[:,:,p].*= ΔJ^dim
    end
end