@inline function get_J2(σ0::SVector{3,T}) where {T}
    P  = (σ0[1]+σ0[2])/T(2.0)
    ξ  = σ0 .- SVector{3,T}(P,P,zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+T(2.0)*ξ[3]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function get_J2(σ0::SVector{6,T}) where {T}
    P  = (σ0[1]+σ0[2]+σ0[3])/T(3.0)
    ξ  = σ0 .- SVector{6,T}(P,P,P,zero(T),zero(T),zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+ξ[3]^2+T(2.0)*ξ[4]^2+T(2.0)*ξ[5]^2+T(2.0)*ξ[6]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function yield_J2(σ,κ::T) where {T}
    ξn,n = get_J2(σ)
    f    = ξn-κ
    return f,n
end

# Rewritten to match the current Point/PointSolidPhase/AbstractConstitutiveModel layout
# (mpts.s.τᵢ/mpts.s.σᵢ are Vector{SVector}, mpts.s.ϵpII is a Matrix, and constitutive
# constants live on mpts.s.cmp[p]::DruckerPrager — see DP.jl for the sibling kernel this
# mirrors). The pre-refactor version of this file indexed flat, non-existent fields
# (mpts.cᵣ[p]/mpts.c₀[p] instead of mpts.s.*, mpts.bᵢⱼ which never existed) and had a
# kernel arg-count mismatch with its call site — this rewrite fixes both as a natural
# side effect of rerouting to mpts.s.cmp[p].
@views @kernel inbounds = true function finite_J2(mpts::Point{T1,T2}; ftol::Real=1e-9,ηmax::Int=20) where {T1,T2} # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mpts.nmp
        cmp = mpts.s.cmp[p]
        # calculate yield function
        κ   = max(cmp.cᵣ,cmp.c₀+cmp.Hp*mpts.s.ϵpII[2,p])
        f,n = yield_J2(mpts.s.τᵢ[p],κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>T2(0.0)
            γ0,σ0,ηit = mpts.s.ϵpII[2,p],mpts.s.τᵢ[p],1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmp.Del*∂f∂σ)
                Δσ   = Δλ .* (cmp.Del*∂f∂σ)
                σ0   = σ0 .- Δσ
                γ0  += Δλ
                κ    = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit += 1
            end
            mpts.s.ϵpII[1,p] = γ0
            mpts.s.τᵢ[p]     = σ0
            # update strain tensor
            mpts.s.ϵᵢⱼ[p]    = mutate(cmp.Del\σ0,T2(0.5),:tensor)
        end
    end
end
@views @kernel inbounds = true function infinitesimal_J2(mpts::Point{T1,T2}; ftol::Real=1e-9,ηmax::Int=20) where {T1,T2} # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mpts.nmp
        cmp = mpts.s.cmp[p]
        # calculate yield function
        κ   = max(cmp.cᵣ,cmp.c₀+cmp.Hp*mpts.s.ϵpII[2,p])
        f,n = yield_J2(mpts.s.σᵢ[p],κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>T2(0.0)
            γ0,σ0,ηit = mpts.s.ϵpII[2,p],mpts.s.σᵢ[p],1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmp.Del*∂f∂σ)
                Δσ   = Δλ .* (cmp.Del*∂f∂σ)
                σ0   = σ0 .- Δσ
                γ0  += Δλ
                κ    = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit += 1
            end
            mpts.s.ϵpII[1,p] = γ0
            mpts.s.σᵢ[p]     = σ0
        end
    end
end
