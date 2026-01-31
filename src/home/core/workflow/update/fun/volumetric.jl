@kernel inbounds = true function ΔJn(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # accumulation
        for nn ∈ 1:mesh.prprt.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            @atom mesh.ΔJ[no]+= mpts.ϕ∂ϕ[nn,p,1]*((mpts.s.ρ[p]*mpts.Ω[p])*mpts.ΔJ[p])  
        end
    end
end
@kernel inbounds = true function ΔJs(mesh::Mesh{T1,T2}) where {T1,T2}
    no = @index(Global)
    if no ≤ mesh.prprt.nno[end] 
        # solve
        if iszero(mesh.s.mᵢ[no])
            mesh.ΔJ[no] = T2(0.0)
        else
            mesh.ΔJ[no]/= mesh.s.mᵢ[no]
        end
    end
end
@views @kernel inbounds = true function ΔJp(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dim::T2) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # mapping back to mpts's
        ΔJ = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            ΔJ += mpts.ϕ∂ϕ[nn,p,1]*mesh.ΔJ[no]/mpts.ΔJ[p]
        end
        # update
        mpts.s.ΔFᵢⱼ[:,:,p].*= ΔJ^dim
    end
end