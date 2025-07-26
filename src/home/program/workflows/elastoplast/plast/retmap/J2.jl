@views function getJ2(σ0,χ,nstr)
    if nstr == 3
        P  = (σ0[1]+σ0[2])/2.0
        ξ  = σ0.-[P,P,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+2.0*ξ[3]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
        n  = ξ./ξn
        q  = sqrt(χ)*ξn
    elseif nstr == 6
        P  = (σ0[1]+σ0[2]+σ0[3])/3.0
        ξ  = σ0.-[P,P,P,0.0,0.0,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2+2.0*ξ[5]^2+2.0*ξ[6]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
        n  = ξ./ξn
        q  = sqrt(χ)*ξn
    end
    return P,q,n,ξn
end
@views function hasyield(σ,κ; χ::Real=3.0/2.0)
    nstr     = length(σ)
    P,q,n,ξn = getJ2(σ,χ,nstr)
    return ξn-κ,n
end
@views @kernel inbounds = true function J2!(mp,ϵIIp,cmp,instr; ftol::Real= 1e-9,ηmax::Int= 20) # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mp.nmp 
        Hp = 0.35*cmp.Hp
        # create aliases
        if instr[:fwrk] == "finite"
            σ = mp.τᵢ
        elseif instr[:fwrk] == "infinitesimal"
            σ = mp.σᵢ
        end
        if first(instr[:nonloc])
            ϵII0 = mp.ϵpII[:,2]
        else
            ϵII0 = mp.ϵpII[:,1]
        end
        # calculate yield function
        κ   = max(mp.cr[p],2.5*mp.c0[p]+cmp.Hp*ϵpII[p,1])
        f,n = hasyield(σ[:,p],κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>0.0 
            γ0,σ0,ηit = copy(mp.ϵpII[p]),copy(σ[:,p]),1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmp.Del*∂f∂σ)
                Δσ   = (Δλ*cmp.Del*∂f∂σ)        
                σ0 .-= Δσ 
                γ0  += Δλ
                κ    = max(mp.cr[p],2.5*mp.c0[p]+cmp.Hp*ϵpII[p,1])
                f    = hasyield(σ[:,p],κ)
                ηit +=1
            end
            mp.ϵpII[p,1] = γ0
            σ[:,p]      .= σ0
            if instr[:fwrk] == "finite"
                # update strain tensor
                mp.ϵᵢⱼ[:,:,p].= mutate(cmp.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green deformation tensor
                λ,n            = eigen(mp.ϵᵢⱼ[:,:,p],sortby=nothing)
                mp.Bᵢⱼ[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
            ηmax = max(ηit,ηmax)
        end        
    end
end