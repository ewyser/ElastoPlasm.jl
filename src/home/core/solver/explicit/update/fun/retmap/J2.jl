# NOTE ON NAMING: these two `_yield_normal` methods used to be called `get_J2`, which
# collided semantically with `get_J2(::AbstractStress)` in `concrete/tensor.jl` — same
# name, different return value (this one returns `(‖ξ‖, n̂)`, the yield-surface normal
# and the norm of the deviator; that one returns the scalar invariant J₂). They were
# deliberately NOT merged during the tensor port: unifying them would have meant either
# silently changing what one call site gets back, or reformulating the return mapping in
# terms of the stored (p,dev) split, which is not canonically trace-free (see the
# `AbstractTensor` docstring). Renamed instead, so the collision is gone and each name
# means one thing.
@inline function _yield_normal(σ0::SVector{3,T}) where {T}
    P  = (σ0[1]+σ0[2])/T(2.0)
    ξ  = σ0 .- SVector{3,T}(P,P,zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+T(2.0)*ξ[3]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function _yield_normal(σ0::SVector{6,T}) where {T}
    P  = (σ0[1]+σ0[2]+σ0[3])/T(3.0)
    ξ  = σ0 .- SVector{6,T}(P,P,P,zero(T),zero(T),zero(T))
    J2 = T(0.5)*(ξ[1]^2+ξ[2]^2+ξ[3]^2+T(2.0)*ξ[4]^2+T(2.0)*ξ[5]^2+T(2.0)*ξ[6]^2) # Borja (2013), p.33
    ξn = sqrt(T(2.0)*J2)
    n  = ξ./ξn
    return ξn,n
end
@inline function yield_J2(σ,κ::T) where {T}
    ξn,n = _yield_normal(σ)
    f    = ξn-κ
    return f,n
end

# Rewritten to match the current Point/PointSolidPhase/AbstractConstitutiveModel layout
# (mpts.s.τᵢⱼ/mpts.s.σᵢⱼ/mpts.s.ϵpII are Vector{SVector}, and constitutive
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
        κ   = max(cmp.cᵣ,cmp.c₀+cmp.Hp*mpts.s.ϵpII[p][2])
        f,n = yield_J2(get_voigt(mpts.s.τᵢⱼ[p]),κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>T2(0.0)
            γ0,σ0,ηit = mpts.s.ϵpII[p][2],get_voigt(mpts.s.τᵢⱼ[p]),1
            Δλacc     = T2(0.0)
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmp.Del*∂f∂σ)
                Δσ   = Δλ .* (cmp.Del*∂f∂σ)
                σ0   = σ0 .- Δσ
                γ0  += Δλ
                Δλacc += Δλ
                κ    = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit += 1
            end
            mpts.s.ϵpII[p]   = SVector{2,T2}(γ0, mpts.s.ϵpII[p][2])
            # Δλ was previously a pure loop-local: it fed Δσ/γ0 and was then thrown away,
            # so mpts.s.Δλ[p] stayed at its zero-initialized value forever under
            # plast.constitutive="J2" — permanently disabling nonlocal.jl's activity gate
            # (`mpts.s.Δλ[p] != 0`). Accumulate it over the CPA iterations and store it.
            mpts.s.Δλ[p]     = Δλacc
            mpts.s.τᵢⱼ[p]     = KirchhoffStress(σ0)
            # update strain tensor: Del\σ0 is already an engineering-Voigt strain
            # vector, so LogarithmicStrain(::SVector) (tensor.jl) does the
            # engineering→tensor conversion and vol/dev split in one step
            mpts.s.ϵᵢⱼ[p]    = LogarithmicStrain(cmp.Del\σ0)
        else
            mpts.s.Δλ[p]     = T2(0.0)
        end
    end
end
@views @kernel inbounds = true function infinitesimal_J2(mpts::Point{T1,T2}; ftol::Real=1e-9,ηmax::Int=20) where {T1,T2} # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mpts.nmp
        cmp = mpts.s.cmp[p]
        # calculate yield function
        κ   = max(cmp.cᵣ,cmp.c₀+cmp.Hp*mpts.s.ϵpII[p][2])
        f,n = yield_J2(get_voigt(mpts.s.σᵢⱼ[p]),κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>T2(0.0)
            γ0,σ0,ηit = mpts.s.ϵpII[p][2],get_voigt(mpts.s.σᵢⱼ[p]),1
            Δλacc     = T2(0.0)
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmp.Del*∂f∂σ)
                Δσ   = Δλ .* (cmp.Del*∂f∂σ)
                σ0   = σ0 .- Δσ
                γ0  += Δλ
                Δλacc += Δλ
                κ    = max(cmp.cᵣ,cmp.c₀+cmp.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit += 1
            end
            mpts.s.ϵpII[p]   = SVector{2,T2}(γ0, mpts.s.ϵpII[p][2])
            # see finite_J2 above — Δλ was never written, silently disabling nonlocal.jl
            mpts.s.Δλ[p]     = Δλacc
            mpts.s.σᵢⱼ[p]     = CauchyStress(σ0)
        else
            mpts.s.Δλ[p]     = T2(0.0)
        end
    end
end
