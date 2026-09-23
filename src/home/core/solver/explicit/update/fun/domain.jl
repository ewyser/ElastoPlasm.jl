@views @kernel inbounds = true function Uᵢᵢ(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        F   = mpts.s.Fᵢⱼ[p]
        λ,n = eigen(Symmetric(F'*F))
        U   = (n*diagm(sqrt.(λ))*n')
        mpts.ℓ[p] = diag(U) .* mpts.ℓ₀[p]
    end
end
