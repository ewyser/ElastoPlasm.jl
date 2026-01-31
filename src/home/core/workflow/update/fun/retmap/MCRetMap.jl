@views function MCRetMap!(mpts,ϵIIp,cmp,fwrkDeform)
    ftol,ηtol,ηmax = 1e-6,1e4,0
    ψ              = 0.5*π/180.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpts.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpts.σ
    end
    for p ∈ 1:mpts.nmp
        ϕ,H,ϵII0 = mpts.ϕ[p],cos(mpts.ϕ[p])*cmp.Hp,ϵIIp[p]
        c0,cr    = mpts.c₀[p]+cmp.Hp*ϵII0,mpts.cᵣ[p]
        if c0<cr c0 = cr end
        σm,τII   = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
        f        = τII+σm*sin(ϕ)-c0*cos(ϕ)    
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,3)
            ηit = 0
            while abs(f)>ftol
                ηit+= 1
                ∂σf = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                         σ[3,p]/τII                     ]
                ∂σg = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                         σ[3,p]/τII                     ] 

                Δγ  = f/(H+∂σf'*cmp.Del*∂σg)
                Δσ  = Δγ*cmp.Del*∂σg
                Δϵ.+= cmp.Del\Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+2*Δϵ[3]^2))
                c0  = mpts.c0[p]+cmp.Hp*ϵII
                if c0<cr c0 = cr end
                σ[:,p].-= Δσ
                σm,τII  = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
                f       = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if ηit>ηtol
                    err_msg = "CPA: η_it>$(ηit): program killed..."
                    throw(error(err_msg))
                end
                ηmax = max(ηit,ηmax)
            end
            mpts.ϵ[:,:,p].= mutate(cmp.Del\σ[:,p],0.5,:tensor)
            mpts.ϵpII[p]  = ϵII 
            if fwrkDeform == :finite
                # update left cauchy green tensor
                λ,n           = eigen(mpts.ϵ[:,:,p],sortby=nothing)
                mpts.b[:,:,p] .= n*diagm(exp.(2.0.*λ))*n'
            end
        end        
    end
    return ηmax::Int64
end