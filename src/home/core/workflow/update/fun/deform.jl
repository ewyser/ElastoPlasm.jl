@views @kernel inbounds = true function deform(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # compute velocity & displacement gradients
        mpts.s.∇vᵢⱼ[:,:,p].= T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end                                     
            for i ∈ 1:mesh.prprt.dim  
                for j ∈ 1:mesh.prprt.dim
                    mpts.s.∇vᵢⱼ[i,j,p]+= ∂N[j]*mesh.s.v[i,no]
                end
            end
        end
        # compute incremental deformation and update
        mpts.s.ΔFᵢⱼ[:,:,p].= mpts.δᵢⱼ+(dt.*mpts.s.∇vᵢⱼ[:,:,p])
        mpts.ΔJ[p]         = det(mpts.s.ΔFᵢⱼ[:,:,p])
        mpts.s.Fᵢⱼ[:,:,p] .= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.Fᵢⱼ[:,:,p]
        mpts.J[p]          = det(mpts.s.Fᵢⱼ[:,:,p])
        # update material point's volume        
        mpts.Ω[p]          = mpts.J[p]*mpts.Ω₀[p]
        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p]        = mpts.s.ρ[p]/mpts.ΔJ[p]     
        mpts.n[p]          = T2(1.0)-T2(1.0)/mpts.J[p]*(T2(1.0)-mpts.n[p])
    end
end

@kernel inbounds = true function deform_fast(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D<:TwoDimension}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # compute velocity & displacement gradients
        ∇vxx,∇vxy,∇vyx,∇vyy = T2(0.0),T2(0.0),T2(0.0),T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end   
            ∇vxx += ∂N[1]*mesh.s.v[1,no]
            ∇vxy += ∂N[1]*mesh.s.v[2,no]
            ∇vyx += ∂N[2]*mesh.s.v[1,no]
            ∇vyy += ∂N[2]*mesh.s.v[2,no]
        end
        # compute incremental deformation tensor
        ΔFxx,ΔFxy = T2(1.0)+(dt*∇vxx),        (dt*∇vxy) 
        ΔFyx,ΔFyy =         (dt*∇vyx),T2(1.0)+(dt*∇vyy) 
        mpts.s.ΔFᵢⱼ[1,1,p], mpts.s.ΔFᵢⱼ[1,2,p] = ΔFxx, ΔFxy
        mpts.s.ΔFᵢⱼ[2,1,p], mpts.s.ΔFᵢⱼ[2,2,p] = ΔFyx, ΔFyy        
        # compute incremental deformation jacobian
        ΔJ   = ΔFxx*ΔFyy-ΔFxy*ΔFyx
        # update deformation gradient
        fxx,fxy = mpts.s.Fᵢⱼ[1,1,p], mpts.s.Fᵢⱼ[1,2,p]
        fyx,fyy = mpts.s.Fᵢⱼ[2,1,p], mpts.s.Fᵢⱼ[2,2,p]
        Fxx,Fxy = ΔFxx*fxx+ΔFxy*fyx, ΔFxx*fxy+ΔFxy*fyy
        Fyx,Fyy = ΔFyx*fxx+ΔFyy*fyx, ΔFyx*fxy+ΔFyy*fyy
        mpts.s.Fᵢⱼ[1,1,p], mpts.s.Fᵢⱼ[1,2,p] = Fxx, Fxy
        mpts.s.Fᵢⱼ[2,1,p], mpts.s.Fᵢⱼ[2,2,p] = Fyx, Fyy
        # compute deformation jacobian
        J           = Fxx*Fyy-Fxy*Fyx
        # update deformation gradient and jacobian
        mpts.ΔJ[p]  = ΔJ
        mpts.J[p]   = J
        # update material point's volume        
        mpts.Ω[p]   = J*mpts.Ω₀[p]
        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p] = mpts.s.ρ[p]/ΔJ
        mpts.n[p]   = T2(1.0)-T2(1.0)/J*(T2(1.0)-mpts.n[p])
    end
end

@kernel inbounds = true function deform_fast(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # compute velocity & displacement gradients
        ∇vxx,∇vxy,∇vxz = T2(0.0),T2(0.0),T2(0.0)
        ∇vyx,∇vyy,∇vyz = T2(0.0),T2(0.0),T2(0.0)
        ∇vzx,∇vzy,∇vzz = T2(0.0),T2(0.0),T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            ∇vxx += ∂N[1]*mesh.s.v[1,no]
            ∇vxy += ∂N[1]*mesh.s.v[2,no]
            ∇vxz += ∂N[1]*mesh.s.v[3,no]
            ∇vyx += ∂N[2]*mesh.s.v[1,no]
            ∇vyy += ∂N[2]*mesh.s.v[2,no]
            ∇vyz += ∂N[2]*mesh.s.v[3,no]
            ∇vzx += ∂N[3]*mesh.s.v[1,no]
            ∇vzy += ∂N[3]*mesh.s.v[2,no]
            ∇vzz += ∂N[3]*mesh.s.v[3,no]
        end
        # compute incremental deformation tensor
        ΔFxx,ΔFxy,ΔFxz = T2(1.0)+(dt*∇vxx),        (dt*∇vxy),        (dt*∇vxz) 
        ΔFyx,ΔFyy,ΔFyz =         (dt*∇vyx),T2(1.0)+(dt*∇vyy),        (dt*∇vyz) 
        ΔFzx,ΔFzy,ΔFzz =         (dt*∇vzx),        (dt*∇vzy),T2(1.0)+(dt*∇vzz)
        
        # compute incremental deformation jacobian
        ΔJ  = ΔFxx*(ΔFyy*ΔFzz-ΔFyz*ΔFzy)
        ΔJ += ΔFxy*(ΔFyx*ΔFzz-ΔFzx*ΔFyz)
        ΔJ += ΔFxz*(ΔFyx*ΔFzy-ΔFzx*ΔFyy)

        # update deformation gradient
        fxx,fxy,fxz = mpts.s.Fᵢⱼ[1,1,p], mpts.s.Fᵢⱼ[1,2,p], mpts.s.Fᵢⱼ[1,3,p]
        fyx,fyy,fyz = mpts.s.Fᵢⱼ[2,1,p], mpts.s.Fᵢⱼ[2,2,p], mpts.s.Fᵢⱼ[2,3,p]
        fzx,fzy,fzz = mpts.s.Fᵢⱼ[3,1,p], mpts.s.Fᵢⱼ[3,2,p], mpts.s.Fᵢⱼ[3,3,p]

        Fxx,Fxy,Fxz = ΔFxx*fxx+ΔFxy*fyx+ΔFxz*fzx, ΔFxx*fxy+ΔFxy*fyy+ΔFxz*fzy, ΔFxx*fxz+ΔFxy*fyz+ΔFxz*fzz
        Fyx,Fyy,Fyz = ΔFyx*fxx+ΔFyy*fyx+ΔFyz*fzx, ΔFyx*fxy+ΔFyy*fyy+ΔFyz*fzy, ΔFyx*fxz+ΔFyy*fyz+ΔFyz*fzz
        Fzx,Fzy,Fzz = ΔFzx*fxx+ΔFzy*fyx+ΔFzz*fzx, ΔFzx*fxy+ΔFzy*fyy+ΔFzz*fzy, ΔFzx*fxz+ΔFzy*fyz+ΔFzz*fzz

        mpts.s.Fᵢⱼ[1,1,p], mpts.s.Fᵢⱼ[1,2,p], mpts.s.Fᵢⱼ[1,3,p] = Fxx, Fxy, Fxz
        mpts.s.Fᵢⱼ[2,1,p], mpts.s.Fᵢⱼ[2,2,p], mpts.s.Fᵢⱼ[2,3,p] = Fyx, Fyy, Fyz
        mpts.s.Fᵢⱼ[3,1,p], mpts.s.Fᵢⱼ[3,2,p], mpts.s.Fᵢⱼ[3,3,p] = Fzx, Fzy, Fzz 

        # compute deformation jacobian
        J  = Fxx*(Fyy*Fzz-Fyz*Fzy)
        J += Fxy*(Fyx*Fzz-Fzx*Fyz)
        J += Fxz*(Fyx*Fzy-Fzx*Fyy)

        # update deformation gradient and jacobian
        mpts.ΔJ[p]  = ΔJ
        mpts.J[p]   = J
        # update material point's volume        
        mpts.Ω[p]   = J*mpts.Ω₀[p]
        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p] = mpts.s.ρ[p]/ΔJ
        mpts.n[p]   = T2(1.0)-T2(1.0)/J*(T2(1.0)-mpts.n[p])
    end
end