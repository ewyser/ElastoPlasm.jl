@views function σTr(σ0,nstr)
    if nstr == 3
        P   = (σ0[1]+σ0[2])/2.0
        τ0  = σ0.-[P,P,0.0]
        τII = sqrt(0.5*(τ0[1]^2+τ0[2]^2)+τ0[3]^2)
    elseif nstr == 6
        P   = (σ0[1]+σ0[2]+σ0[3])/3.0
        τ0  = σ0.-[P,P,P,0.0,0.0,0.0]
        τII = sqrt(0.5*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    end
    return P,τ0,τII
end
@views function materialParam(ϕ,ψ,c,nstr)
    if nstr == 3
        η   = 3.0*tan(ϕ)/(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
        ηB  = 3.0*tan(ψ)/(sqrt(9.0+12.0*tan(ψ)*tan(ψ)))
        ξ   = 3.0*c     /(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
    elseif nstr == 6
        η   = 6.0  *sin(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ)))
        ηB  = 6.0  *sin(ψ)/(sqrt(3.0)*(3.0+sin(ψ)))
        ξ   = 6.0*c*cos(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ))) 
    end
    return η,ηB,ξ
end
@views function σn(Pn,τ0,τn,τII,nstr)
    if nstr == 3
        σn = τ0.*(τn/τII).+[Pn,Pn,0.0]
    elseif nstr == 6
        σn = τ0.*(τn/τII).+[Pn,Pn,Pn,0.0,0.0,0.0]
    end
    return σn 
end
@views @kernel inbounds = true function DP!(mpD,cmParam,instr)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.Δλ[p] = 0.0
        ψ,nstr   = 0.0*π/180.0,size(mpD.σᵢ,1)
        # create an alias for stress tensor
        if instr[:fwrk] == "finite"
            σ = mpD.τᵢ
        elseif instr[:fwrk] == "infinitesimal"
            σ = mpD.σᵢ
        end
        if first(instr[:nonloc])
            ϵII0 = mpD.ϵpII[:,2]
        else
            ϵII0 = mpD.ϵpII[:,1]
        end
        # closed-form solution return-mapping for D-P
        c   = mpD.c₀[p]+cmParam.Hp*ϵII0[p]
        if c<mpD.cᵣ[p] c = mpD.cᵣ[p] end
        P,τ0,τII = σTr(σ[:,p],nstr)
        η,ηB,ξ   = materialParam(mpD.ϕ[p],ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(1.0+η^2)-η,τII-τP-(sqrt(1.0+η^2))*(P-σm)  
        if fs>0.0 && P<σm || h>0.0 && P>=σm
            Δλ        = fs/(cmParam.Gc+cmParam.Kc*η*ηB)
            mpD.Δλ[p] = Δλ
            Pn,τn     = P-cmParam.Kc*ηB*Δλ,ξ-η*(P-cmParam.Kc*ηB*Δλ)
            σ[:,p]   .= σn(Pn,τ0,τn,τII,nstr)
            if instr[:fwrk] == "finite"
                mpD.ϵᵢⱼ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n            = eigen(mpD.ϵᵢⱼ[:,:,p],sortby=nothing)
                mpD.Bᵢⱼ[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
            mpD.ϵpII[p,1]+= Δλ*sqrt(1/3+2/9*ηB^2)
        end
        if h<=0.0 && P>=σm
            Δλ        = (P-σm)/cmParam.Kc
            mpD.Δλ[p] = Δλ
            Pn        = σm-P
            σ[:,p]   .= σn(Pn,τ0,0.0,τII,nstr)
            if instr[:fwrk] == "finite"
                mpD.ϵᵢⱼ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n            = eigen(mpD.ϵᵢⱼ[:,:,p],sortby=nothing)
                mpD.Bᵢⱼ[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
            mpD.ϵpII[p,1]+= sqrt(2.0)*Δλ/3.0
        end
    end
end