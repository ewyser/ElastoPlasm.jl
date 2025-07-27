@views function get_J2(σ0; χ::Real=3.0/2.0)
    if length(σ0) == 3
        P  = (σ0[1]+σ0[2])/2.0
        ξ  = σ0.-[P,P,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+2.0*ξ[3]^2) # Borja (2013), p.33
    elseif length(σ0) == 6
        P  = (σ0[1]+σ0[2]+σ0[3])/3.0
        ξ  = σ0.-[P,P,P,0.0,0.0,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2+2.0*ξ[5]^2+2.0*ξ[6]^2) # Borja (2013), p.33
    end
        ξn = sqrt(2.0*J2) 
        n  = ξ./ξn
        #q  = sqrt(χ)*ξn
    return ξn,n
end
@views function yield_J2(σ,κ)
    ξn,n = get_J2(σ)
    f    = ξn-κ
    return f,n
end
@views @kernel inbounds = true function finite_J2(mp,cmpr,instr; ftol::Real= 1e-9,ηmax::Int= 20) # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mp.nmp 
        # calculate yield function
        κ   = max(mp.cᵣ[p],mp.c₀[p]+cmpr.Hp*mp.ϵpII[2,p])
        f,n = yield_J2(mp.τᵢ[:,p],κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>0.0 
            γ0,σ0,ηit = copy(mp.ϵpII[2,p]),copy(mp.τᵢ[:,p]),1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmpr.Del*∂f∂σ)
                Δσ   = (Δλ*cmpr.Del*∂f∂σ)        
                σ0 .-= Δσ 
                γ0  += Δλ
                κ    = max(mp.cᵣ[p],mp.c₀[p]+cmpr.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit +=1
            end
            mp.ϵpII[p,1] = γ0
            mp.τᵢ[:,p]  .= σ0
            # update strain tensor
            mp.ϵᵢⱼ[:,:,p].= mutate(cmpr.Del\mp.τᵢ[:,p],0.5,:tensor)
            # update left cauchy green deformation tensor
            λ,n           = eigen(mp.ϵᵢⱼ[:,:,p],sortby=nothing)
            mp.Bᵢⱼ[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            # update plastic corrector increment
            ηmax = max(ηit,ηmax)
        end        
    end
end
@views @kernel inbounds = true function infinitesimal_J2(mp,cmpr,instr; ftol::Real= 1e-9,ηmax::Int= 20) # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mp.nmp 
        # calculate yield function
        κ   = max(mp.cᵣ[p],mp.c₀[p]+cmpr.Hp*mp.ϵpII[2,p])
        f,n = yield_J2(mp.σᵢ[:,p],κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>0.0 
            γ0,σ0,ηit = copy(mp.ϵpII[2,p]),copy(mp.σᵢ[:,p]),1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmpr.Del*∂f∂σ)
                Δσ   = (Δλ*cmpr.Del*∂f∂σ)        
                σ0 .-= Δσ 
                γ0  += Δλ
                κ    = max(mp.cᵣ[p],mp.c₀[p]+cmpr.Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit +=1
            end
            mp.ϵpII[p,1] = γ0
            mp.σᵢ[:,p]  .= σ0
            # update plastic corrector increment
            ηmax = max(ηit,ηmax)
        end        
    end
end