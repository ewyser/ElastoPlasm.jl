@kernel inbounds = true function deform(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # accumulate ∇vᵢⱼ into a local mutable buffer then convert
        TM = eltype(mpts.s.∇vᵢⱼ)
        ∇v = zeros(MMatrix{size(TM,1),size(TM,2),T2})
        for nn ∈ 1:mesh.prprt.nn
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            for i ∈ 1:(length(mesh.prprt.nel)-1)
                for j ∈ 1:(length(mesh.prprt.nel)-1)
                    ∇v[i,j] += ∂N[j]*mesh.s.v[no][i]
                end
            end
        end
        ∇vᵢⱼ   = TM(∇v)
        ΔFᵢⱼ   = TM(I) + dt * ∇vᵢⱼ
        Fᵢⱼ    = ΔFᵢⱼ * mpts.s.Fᵢⱼ[p]
        mpts.s.∇vᵢⱼ[p] = ∇vᵢⱼ
        mpts.s.ΔFᵢⱼ[p] = ΔFᵢⱼ
        mpts.s.Fᵢⱼ[p]  = Fᵢⱼ
        ΔJ = det(ΔFᵢⱼ)
        J  = det(Fᵢⱼ)
        mpts.ΔJ[p] = ΔJ
        mpts.J[p]  = J
        mpts.Ω[p]  = J * mpts.Ω₀[p]
        mpts.s.ρ[p] = mpts.s.ρ[p] / ΔJ
        mpts.n[p]   = T2(1.0) - T2(1.0)/J*(T2(1.0)-mpts.n[p])
    end
end

@kernel inbounds = true function deform_fast(mpts::Point{T1,T2,2},mesh::Mesh{T1,T2,2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # compute velocity & displacement gradients
        ∇vxx,∇vxy,∇vyx,∇vyy = T2(0.0),T2(0.0),T2(0.0),T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end   
            ∇vxx += ∂N[1]*mesh.s.v[no][1]
            ∇vxy += ∂N[1]*mesh.s.v[no][2]
            ∇vyx += ∂N[2]*mesh.s.v[no][1]
            ∇vyy += ∂N[2]*mesh.s.v[no][2]
        end
        # compute incremental deformation tensor (SMatrix, column-major: [col1..., col2...])
        ΔFxx,ΔFxy = T2(1.0)+(dt*∇vxx),        (dt*∇vxy)
        ΔFyx,ΔFyy =         (dt*∇vyx),T2(1.0)+(dt*∇vyy)
        ΔFᵢⱼ = SMatrix{2,2,T2,4}(ΔFxx, ΔFyx, ΔFxy, ΔFyy)
        mpts.s.ΔFᵢⱼ[p] = ΔFᵢⱼ
        ΔJ = ΔFxx*ΔFyy-ΔFxy*ΔFyx
        F  = mpts.s.Fᵢⱼ[p]
        Fxx,Fxy = ΔFxx*F[1,1]+ΔFxy*F[2,1], ΔFxx*F[1,2]+ΔFxy*F[2,2]
        Fyx,Fyy = ΔFyx*F[1,1]+ΔFyy*F[2,1], ΔFyx*F[1,2]+ΔFyy*F[2,2]
        mpts.s.Fᵢⱼ[p] = SMatrix{2,2,T2,4}(Fxx, Fyx, Fxy, Fyy)
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

@kernel inbounds = true function deform_fast(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},dt::T2) where {T1,T2}
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
            ∇vxx += ∂N[1]*mesh.s.v[no][1]
            ∇vxy += ∂N[1]*mesh.s.v[no][2]
            ∇vxz += ∂N[1]*mesh.s.v[no][3]
            ∇vyx += ∂N[2]*mesh.s.v[no][1]
            ∇vyy += ∂N[2]*mesh.s.v[no][2]
            ∇vyz += ∂N[2]*mesh.s.v[no][3]
            ∇vzx += ∂N[3]*mesh.s.v[no][1]
            ∇vzy += ∂N[3]*mesh.s.v[no][2]
            ∇vzz += ∂N[3]*mesh.s.v[no][3]
        end
        ΔFxx,ΔFxy,ΔFxz = T2(1.0)+(dt*∇vxx),        (dt*∇vxy),        (dt*∇vxz)
        ΔFyx,ΔFyy,ΔFyz =         (dt*∇vyx),T2(1.0)+(dt*∇vyy),        (dt*∇vyz)
        ΔFzx,ΔFzy,ΔFzz =         (dt*∇vzx),        (dt*∇vzy),T2(1.0)+(dt*∇vzz)
        # SMatrix is column-major: SMatrix{3,3}(col1..., col2..., col3...)
        ΔFᵢⱼ = SMatrix{3,3,T2,9}(ΔFxx,ΔFyx,ΔFzx, ΔFxy,ΔFyy,ΔFzy, ΔFxz,ΔFyz,ΔFzz)
        mpts.s.ΔFᵢⱼ[p] = ΔFᵢⱼ
        ΔJ  = ΔFxx*(ΔFyy*ΔFzz-ΔFyz*ΔFzy)
        ΔJ += ΔFxy*(ΔFyx*ΔFzz-ΔFzx*ΔFyz)
        ΔJ += ΔFxz*(ΔFyx*ΔFzy-ΔFzx*ΔFyy)
        F  = mpts.s.Fᵢⱼ[p]
        Fxx,Fxy,Fxz = ΔFxx*F[1,1]+ΔFxy*F[2,1]+ΔFxz*F[3,1], ΔFxx*F[1,2]+ΔFxy*F[2,2]+ΔFxz*F[3,2], ΔFxx*F[1,3]+ΔFxy*F[2,3]+ΔFxz*F[3,3]
        Fyx,Fyy,Fyz = ΔFyx*F[1,1]+ΔFyy*F[2,1]+ΔFyz*F[3,1], ΔFyx*F[1,2]+ΔFyy*F[2,2]+ΔFyz*F[3,2], ΔFyx*F[1,3]+ΔFyy*F[2,3]+ΔFyz*F[3,3]
        Fzx,Fzy,Fzz = ΔFzx*F[1,1]+ΔFzy*F[2,1]+ΔFzz*F[3,1], ΔFzx*F[1,2]+ΔFzy*F[2,2]+ΔFzz*F[3,2], ΔFzx*F[1,3]+ΔFzy*F[2,3]+ΔFzz*F[3,3]
        mpts.s.Fᵢⱼ[p] = SMatrix{3,3,T2,9}(Fxx,Fyx,Fzx, Fxy,Fyy,Fzy, Fxz,Fyz,Fzz)

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