@views @kernel inbounds = true function undeformed(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        mpts.ℓ[:,p].= mpts.ℓ₀[:,p]
    end
end
@views @kernel inbounds = true function ΔUᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpts.ΔFᵢⱼ[:,:,p]'*mpts.ΔFᵢⱼ[:,:,p],sortby=nothing)
        ΔU         = (n*diagm(sqrt.(λ))*n')
        mpts.ℓ[:,p].= diag(ΔU).*mpts.ℓ[:,p]
    end
end
@views @kernel inbounds = true function Uᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpts.Fᵢⱼ[:,:,p]'*mpts.Fᵢⱼ[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mpts.ℓ[:,p].= diag(U).*mpts.ℓ₀[:,p]
    end
end
@views @kernel inbounds = true function detΔFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[:,p].= mpts.ℓ[:,p].*det(mpts.ΔFᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function detFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[:,p].= mpts.ℓ₀[:,p].*det(mpts.Fᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function ΔFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[:,p].= mpts.ℓ[:,p].*diag(mpts.ΔFᵢⱼ[:,:,p])
    end
end
@views @kernel inbounds = true function Fᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[:,p].= mpts.ℓ₀[:,p].*diag(mpts.Fᵢⱼ[:,:,p])
    end
end
