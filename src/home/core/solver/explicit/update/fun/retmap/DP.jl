@views function σTr(σ0::AbstractVector{T2},nstr::T1) where {T1,T2}
    if nstr == T1(3)
        P   = (σ0[1]+σ0[2])/T2(2.0)
        τ0  = σ0.-[P,P,T2(0.0)]
        τII = sqrt(T2(0.5)*(τ0[1]^2+τ0[2]^2)+τ0[3]^2)
    elseif nstr == T1(6)
        P   = (σ0[1]+σ0[2]+σ0[3])/T2(3.0)
        τ0  = σ0.-[P,P,P,T2(0.0),T2(0.0),T2(0.0)]
        τII = sqrt(T2(0.5)*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    end
    return P,τ0,τII
end
@views function materialParam(ϕ::T2,ψ::T2,c::T2,nstr::T1) where {T1,T2}
    if nstr == T1(3)
        η   = T2(3.0)*tan(ϕ)/(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
        ηB  = T2(3.0)*tan(ψ)/(sqrt(T2(9.0)+T2(12.0)*tan(ψ)*tan(ψ)))
        ξ   = T2(3.0)*c     /(sqrt(T2(9.0)+T2(12.0)*tan(ϕ)*tan(ϕ)))
    elseif nstr == T1(6)
        η   = T2(6.0)*sin(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ)))
        ηB  = T2(6.0)*sin(ψ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ψ)))
        ξ   = T2(6.0)*c*cos(ϕ)/(sqrt(T2(3.0))*(T2(3.0)+sin(ϕ))) 
    end
    return η,ηB,ξ
end
@views function σn(Pn::T2,τ0,τn::T2,τII::T2,nstr::T1) where {T1,T2}
    if nstr == T1(3)
        σn = τ0.*(τn/τII).+[Pn,Pn,T2(0.0)]
    elseif nstr == T1(6)
        σn = τ0.*(τn/τII).+[Pn,Pn,Pn,T2(0.0),T2(0.0),T2(0.0)]
    end
    return σn 
end
@views @kernel inbounds = true function finite_DP(mpts::Point{T1,T2},cmp::NamedTuple) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        mpts.s.Δλ[p] = T2(0.0)
        ψ,nstr   = T2(0.0*π/180.0),length(mpts.s.τᵢ[1])

        # closed-form solution return-mapping for D-P
        c = mpts.s.c₀[p]+cmp.Hp*mpts.s.ϵpII[2,p]
        if c<mpts.s.cᵣ[p] 
            c = mpts.s.cᵣ[p] 
        end
        P,τ0,τII = σTr(mpts.s.τᵢ[p],nstr)
        η,ηB,ξ   = materialParam(mpts.s.ϕ₀[p],ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T2(1.0)+η^2)-η,τII-τP-(sqrt(T2(1.0)+η^2))*(P-σm)  
        if fs>T2(0.0) && P<σm || h>T2(0.0) && P≥σm 
            Δλ             = fs/(cmp.Gc+cmp.Kc*η*ηB)
            Pn,τn          = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            mpts.s.Δλ[p]     = Δλ
            mpts.s.τᵢ[p]   = eltype(mpts.s.τᵢ)(σn(Pn,τ0,τn,τII,nstr))
            mpts.s.ϵpII[1,p]+= Δλ*sqrt(T2(1/3)+T2(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ             = (P-σm)/cmp.Kc
            Pn             = σm-P
            mpts.s.Δλ[p]     = Δλ
            mpts.s.τᵢ[p]   = eltype(mpts.s.τᵢ)(σn(Pn,τ0,T2(0.0),τII,nstr))
            mpts.s.ϵpII[1,p]+= sqrt(T2(2.0))*Δλ/T2(3.0)
        end
        # update logarithmicstrain tensor
        ϵ_mat = mutate(cmp.Del\mpts.s.τᵢ[p],T2(0.5),:tensor)
        mpts.s.ϵᵢⱼ[p] = eltype(mpts.s.ϵᵢⱼ)(ϵ_mat)
    end
end
@views @kernel inbounds = true function infinitesimal_DP(mpts::Point{T1,T2},cmp::NamedTuple) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        mpts.s.Δλ[p] = T2(0.0)
        ψ,nstr   = T2(0.0*π/180.0),length(mpts.s.σᵢ[1])

        # closed-form solution return-mapping for D-P
        c   = mpts.s.c₀[p]+cmp.Hp*mpts.s.ϵpII[2,p]
        if c<mpts.s.cᵣ[p] 
            c = mpts.s.cᵣ[p] 
        end
        P,τ0,τII = σTr(mpts.s.σᵢ[p],nstr)
        η,ηB,ξ   = materialParam(mpts.s.ϕ₀[p],ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T2(1.0)+η^2)-η,τII-τP-(sqrt(T2(1.0)+η^2))*(P-σm)  
        if fs>T2(0.0) && P<σm || h>T2(0.0) && P≥σm 
            Δλ             = fs/(cmp.Gc+cmp.Kc*η*ηB)
            mpts.s.Δλ[p]     = Δλ
            Pn,τn          = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            mpts.s.σᵢ[p]   = eltype(mpts.s.σᵢ)(σn(Pn,τ0,τn,τII,nstr))
            mpts.s.ϵpII[1,p]+= Δλ*sqrt(T2(1/3)+T2(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ             = (P-σm)/cmp.Kc
            mpts.s.Δλ[p]     = Δλ
            Pn             = σm-P
            mpts.s.σᵢ[p]   = eltype(mpts.s.σᵢ)(σn(Pn,τ0,T2(0.0),τII,nstr))
            mpts.s.ϵpII[1,p]+= sqrt(T2(2.0))*Δλ/T2(3.0)
        end
    end
end