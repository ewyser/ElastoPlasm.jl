@views @kernel inbounds = true function undeformed(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.ℓ[p,:].= mpD.ℓ₀[p,:]
    end
end
@views @kernel inbounds = true function Uᵢᵢ(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpD.Fᵢⱼ[:,:,p]'*mpD.Fᵢⱼ[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mpD.ℓ[p,:].= diag(U).*mpD.ℓ₀[p,:]
    end
end
@views @kernel inbounds = true function detFᵢᵢ(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpD.ℓ[p,:].= mpD.ℓ₀[p,:].*det(mpD.Fᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function Fᵢᵢ(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpD.ℓ[p,:].= mpD.ℓ₀[p,:].*diag(mpD.Fᵢⱼ[:,:,p])
    end
end
function init_domain(dim::Number,basis::NamedTuple)
    if basis[:which] == "gimpm"
        if basis[:how] == "undeformed"
            update! = undeformed(CPU())
        elseif basis[:how] == "detFii"
            update! = detFᵢᵢ(CPU())
        elseif basis[:how] == "Fii"
            update! = Fᵢᵢ(CPU())
        elseif basis[:how] == "Uii"
            update! = Uᵢᵢ(CPU())
        end
    else
        update! = nothing
    end
    return update! 
end
function domain!(mpD,instr)
    if instr[:basis][:which] == "gimpm"
        instr[:cairn].update!(mpD; ndrange=mpD.nmp);sync(CPU())
    end
    return nothing
end 