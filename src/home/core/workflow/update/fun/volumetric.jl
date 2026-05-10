@kernel inbounds = true function ΔJn(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # buffering
        ms,ΔJ = mpts.s.ρ[p]*mpts.Ω[p], mpts.ΔJ[p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end 
            # accumulation
            @atom mesh.ΔJ[no]+= N*ms*ΔJ  
        end
    end
end
@kernel inbounds = true function ΔJp(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},dim::T2) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # mapping back to mpts's
        ΔJ,ΔJ₀⁻¹ = T2(0.0),T2(1.0)/mpts.ΔJ[p] 
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            ΔJ += N*(mesh.ΔJ[no]/mesh.s.m[no])*ΔJ₀⁻¹
        end
        # update
        for i ∈ 1:mesh.prprt.dim  
            for j ∈ 1:mesh.prprt.dim
                mpts.s.ΔFᵢⱼ[i,j,p]*= ΔJ^dim
            end
        end
    end
end