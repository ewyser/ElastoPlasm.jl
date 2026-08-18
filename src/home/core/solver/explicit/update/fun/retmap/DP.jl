@inline function σTr(σ0::SVector{3,T}) where {T}
    P   = (σ0[1]+σ0[2])/2.0
    τ0  = σ0 .- SVector{3,T}(P,P,zero(T))
    τII = sqrt(0.5 * (τ0[1]^2 + τ0[2]^2) + τ0[3]^2)
    return P,τ0,τII
end
@inline function σTr(σ0::SVector{6,T}) where {T}
    P   = (σ0[1]+σ0[2]+σ0[3])/3.0
    τ0  = σ0 .- SVector{6,T}(P,P,P,zero(T),zero(T),zero(T))
    τII = sqrt(0.5*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    return P,τ0,τII
end
@inline function materialParam(ϕ::T2,ψ::T2,c::T2,nstr::T1) where {T1,T2}
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

@inline function σn(Pn::T,τ0::SVector{3,T},τn::T,τII::T) where {T}
    return τ0 .* (τn / τII) .+ SVector{3,T}(Pn,Pn,zero(T))
end
@inline function σn(Pn::T,τ0::SVector{6,T},τn::T,τII::T) where {T}
    return τ0 .* (τn / τII) .+ SVector{6,T}(Pn,Pn,Pn,zero(T),zero(T),zero(T))
end

@inline function drucker_prager(τᵢ::SVector{D,T},ϵᵢⱼ::SMatrix{S,S,T,L},ϵpII::MVector{2,T},ϕ₀::T,c₀::T,cᵣ::T,cmp::NamedTuple) where {D,S,T,L}
        Δλ     = T(0.0)
        ψ,nstr = T(0.0*π/180.0),length(τᵢ)

        # closed-form solution return-mapping for D-P
        c = c₀+cmp.Hp*ϵpII[2]
        if c<cᵣ 
            c = cᵣ
        end
        P,τ0,τII = σTr(τᵢ)
        η,ηB,ξ   = materialParam(ϕ₀,ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T(1.0)+η^2)-η,τII-τP-(sqrt(T(1.0)+η^2))*(P-σm)
        if fs>T(0.0) && P<σm || h>T(0.0) && P≥σm 
            Δλ       = fs/(cmp.Gc+cmp.Kc*η*ηB)
            Pn,τn    = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            τᵢ       = σn(Pn,τ0,τn,τII)
            ϵpII[1] += Δλ*sqrt(T(1/3)+T(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ      = (P-σm)/cmp.Kc
            Pn      = σm-P
            τᵢ      = (σn(Pn,τ0,T(0.0),τII))
            ϵpII[1]+= sqrt(T(2.0))*Δλ/T(3.0)
        end
        # update logarithmic strain tensor
        if Δλ > T(0.0)
            #=
            τᵢⱼ = SMatrix{2,2}(τᵢ[1], τᵢ[3], τᵢ[3], τᵢ[2])
            P   = tr(τᵢⱼ) / T(3.0)
            dev = τᵢⱼ .- SMatrix{2,2,T,4}(I) * P
            ϵᵢⱼ = P / (3.0 * cmp.Kc) * SMatrix{2,2,T,4}(I) + dev / (2.0 * cmp.Gc)
            =#
            ϵᵢⱼ = mutate(cmp.Del\τᵢ,T(0.5),:tensor)

        end
        return ϵᵢⱼ,τᵢ,Δλ,ϵpII
end

@views @kernel inbounds = true function finite_DP(mpts::Point{T1,T2},cmp::NamedTuple) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        ϵᵢⱼ,τᵢ,Δλ,ϵpII = drucker_prager(mpts.s.τᵢ[p],mpts.s.ϵᵢⱼ[p],MVector{2,T2}(mpts.s.ϵpII[1,p],mpts.s.ϵpII[2,p]),mpts.s.ϕ₀[p],mpts.s.c₀[p],mpts.s.cᵣ[p],cmp)
        if Δλ > T2(0.0)
            mpts.s.ϵᵢⱼ[p]     = ϵᵢⱼ
            mpts.s.τᵢ[p]      = τᵢ
            mpts.s.Δλ[p]      = Δλ
            mpts.s.ϵpII[:,p] .= ϵpII
        end
        #=
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
        =#
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
        P,τ0,τII = σTr(mpts.s.σᵢ[p])
        η,ηB,ξ   = materialParam(mpts.s.ϕ₀[p],ψ,c,nstr)
        σm,τP    = ξ/η,ξ-η*(ξ/η)
        fs,ft    = τII+η*P-ξ,P-σm         
        αP,h     = sqrt(T2(1.0)+η^2)-η,τII-τP-(sqrt(T2(1.0)+η^2))*(P-σm)  
        if fs>T2(0.0) && P<σm || h>T2(0.0) && P≥σm 
            Δλ             = fs/(cmp.Gc+cmp.Kc*η*ηB)
            mpts.s.Δλ[p]     = Δλ
            Pn,τn          = P-cmp.Kc*ηB*Δλ,ξ-η*(P-cmp.Kc*ηB*Δλ)
            mpts.s.σᵢ[p]   = eltype(mpts.s.σᵢ)(σn(Pn,τ0,τn,τII))
            mpts.s.ϵpII[1,p]+= Δλ*sqrt(T2(1/3)+T2(2/9)*ηB^2)
        end
        if h≤0.0 && P≥σm
            Δλ             = (P-σm)/cmp.Kc
            mpts.s.Δλ[p]     = Δλ
            Pn             = σm-P
            mpts.s.σᵢ[p]   = eltype(mpts.s.σᵢ)(σn(Pn,τ0,T2(0.0),τII))
            mpts.s.ϵpII[1,p]+= sqrt(T2(2.0))*Δλ/T2(3.0)
        end
    end
end