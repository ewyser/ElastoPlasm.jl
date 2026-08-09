@views @kernel inbounds = true function undeformed(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        mpts.ℓ[p] = mpts.ℓ₀[p]
    end
end
@views @kernel inbounds = true function ΔUᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        ΔF = mpts.s.ΔFᵢⱼ[p]
        λ,n = eigen(ΔF'ΔF)
        ΔU = (n*diagm(sqrt.(λ))*n')
        mpts.ℓ[p] = eltype(mpts.ℓ)(diag(ΔU)) .* mpts.ℓ[p]
    end
end
@views @kernel inbounds = true function Uᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        F = mpts.s.Fᵢⱼ[p]
        λ,n = eigen(F'F)
        U = (n*diagm(sqrt.(λ))*n')
        mpts.ℓ[p] = eltype(mpts.ℓ)(diag(U)) .* mpts.ℓ₀[p]
    end
end
@views @kernel inbounds = true function detΔFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[p] = mpts.ℓ[p] .* det(mpts.s.ΔFᵢⱼ[p])
    end
end
@views @kernel inbounds = true function detFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[p] = mpts.ℓ₀[p] .* det(mpts.s.Fᵢⱼ[p])
    end
end
@views @kernel inbounds = true function ΔFᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[p] = mpts.ℓ[p] .* diag(mpts.s.ΔFᵢⱼ[p])
    end
end
@views @kernel inbounds = true function Fᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using diagonal terms of deformation gradient
        mpts.ℓ[p] = mpts.ℓ₀[p] .* diag(mpts.s.Fᵢⱼ[p])
    end
end
