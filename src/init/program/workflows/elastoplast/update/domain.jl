@views @kernel inbounds = true function undeformed(mp)
    p = @index(Global)
    if p≤mp.nmp 
        mp.ℓ[:,p].= mp.ℓ₀[:,p]
    end
end
@views @kernel inbounds = true function ΔUᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mp.ΔFᵢⱼ[:,:,p]'*mp.ΔFᵢⱼ[:,:,p],sortby=nothing)
        ΔU         = (n*diagm(sqrt.(λ))*n')
        mp.ℓ[:,p].= diag(ΔU).*mp.ℓ[:,p]
    end
end
@views @kernel inbounds = true function Uᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mp.Fᵢⱼ[:,:,p]'*mp.Fᵢⱼ[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mp.ℓ[:,p].= diag(U).*mp.ℓ₀[:,p]
    end
end
@views @kernel inbounds = true function detΔFᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mp.ℓ[:,p].= mp.ℓ[:,p].*det(mp.ΔFᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function detFᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mp.ℓ[:,p].= mp.ℓ₀[:,p].*det(mp.Fᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function ΔFᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mp.ℓ[:,p].= mp.ℓ[:,p].*diag(mp.ΔFᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function Fᵢᵢ(mp)
    p = @index(Global)
    if p≤mp.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mp.ℓ[:,p].= mp.ℓ₀[:,p].*diag(mp.Fᵢⱼ[:,:,p])
    end
end
