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
@views @kernel inbounds = true function J2!(mp,cmpr,instr; ftol::Real= 1e-9,ηmax::Int= 20) # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mp.nmp 
        Hp = cmpr.Hp
        cᵣ = mp.cᵣ[p]
        c₀ = mp.c₀[p]
        # create aliases
        if instr[:fwrk][:deform] == "finite"
            σ = mp.τᵢ[:,p]
        elseif instr[:fwrk][:deform] == "infinitesimal"
            σ = mp.σᵢ[:,p]
        end
        if instr[:nonloc][:status]
            ϵpII = mp.ϵpII[2,p]
        else
            ϵpII = mp.ϵpII[1,p]
        end
        # calculate yield function
        κ   = max(cᵣ,c₀+Hp*ϵpII)
        f,n = yield_J2(σ,κ)
        # return mapping using CPA (non-quadratic convergence)
        if f>0.0 
            γ0,σ0,ηit = copy(ϵpII),copy(σ),1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*cmpr.Del*∂f∂σ)
                Δσ   = (Δλ*cmpr.Del*∂f∂σ)        
                σ0 .-= Δσ 
                γ0  += Δλ
                κ    = max(cᵣ,c₀+Hp*γ0)
                f,n  = yield_J2(σ0,κ)
                ηit +=1
            end
            mp.ϵpII[p,1] = γ0
            if instr[:fwrk][:deform] == "finite"
                mp.τᵢ[:,p] .= σ0
                # update strain tensor
                mp.ϵᵢⱼ[:,:,p].= mutate(cmpr.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green deformation tensor
                λ,n            = eigen(mp.ϵᵢⱼ[:,:,p],sortby=nothing)
                mp.Bᵢⱼ[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            elseif instr[:fwrk][:deform] == "infinitesimal"
                mp.σᵢ[:,p] .= σ0
            end   
            ηmax = max(ηit,ηmax)
            #==#
        end        
    end
end